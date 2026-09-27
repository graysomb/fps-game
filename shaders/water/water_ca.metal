#include <metal_stdlib>
using namespace metal;

struct WaterUniforms { int mode; int cell_count; int padding0; int padding1; };
constant uint WATER_MAX = 65535u;
constant uint LATERAL_SUPPORT_MASS = 60000u;
constant int VELOCITY_MAX = 127, VELOCITY_GRAVITY = 12;
constant int IMPACT_THRESHOLD = 16, HORIZONTAL_KICK = 8;
constant int SIZE_X = 80, SIZE_Y = 80, SIZE_Z = 80;
constant int LAYER = SIZE_X * SIZE_Z, TILE_PAIRS = 256;

inline uint load_mass(device const uint *packed, int i) {
    return (packed[i >> 1] >> ((i & 1) * 16)) & WATER_MAX;
}
inline bool mask_bit(device const uint *mask, int i) {
    return ((mask[i >> 5] >> (i & 31)) & 1u) != 0u;
}
inline bool blocked(int i, device const uint *sm, device const uint *dm) {
    return mask_bit(sm, i) || mask_bit(dm, i);
}
inline int tile_for_index(int i) {
    int y = i / LAYER, rem = i - y * LAYER;
    int z = rem / SIZE_X, x = rem - z * SIZE_X;
    return (x / 8) + 10 * ((z / 8) + 10 * (y / 8));
}
inline bool tile_active(device const uint *active, int tile) {
    return ((active[tile >> 5] >> (tile & 31)) & 1u) != 0u;
}
inline int unpack_component(uint packed, int shift) {
    int value = int((packed >> shift) & 255u);
    return value > 127 ? value - 256 : value;
}
inline int3 unpack_velocity(uint packed) {
    return int3(unpack_component(packed, 0), unpack_component(packed, 8),
                unpack_component(packed, 16));
}
inline uint pack_velocity(int3 velocity) {
    velocity = clamp(velocity, int3(-VELOCITY_MAX), int3(VELOCITY_MAX));
    return uint(velocity.x & 255) | (uint(velocity.y & 255) << 8) |
           (uint(velocity.z & 255) << 16);
}
inline int3 accelerated_velocity(uint packed, uint mass) {
    if (mass == 0u) return int3(0);
    int3 velocity = unpack_velocity(packed);
    velocity.y = max(-VELOCITY_MAX, velocity.y - VELOCITY_GRAVITY);
    return velocity;
}
inline uint down_flow(int i, int y, device const uint *mass, device const uint *velocity,
                      device const uint *sm, device const uint *dm,
                      device const uint *active) {
    if (blocked(i, sm, dm) || y == 0) return 0u;
    if (accelerated_velocity(velocity[i], load_mass(mass, i)).y > 0) return 0u;
    int below = i - LAYER;
    if (blocked(below, sm, dm) || !tile_active(active, tile_for_index(i)) ||
        !tile_active(active, tile_for_index(below))) return 0u;
    return min(load_mass(mass, i), WATER_MAX - load_mass(mass, below));
}
inline uint gravity_mass(int i, int y, device const uint *mass, device const uint *velocity,
                         device const uint *sm, device const uint *dm,
                         device const uint *active) {
    if (blocked(i, sm, dm)) return load_mass(mass, i);
    uint outgoing = down_flow(i, y, mass, velocity, sm, dm, active);
    uint incoming = y + 1 < SIZE_Y ? down_flow(i + LAYER, y + 1, mass, velocity, sm, dm, active) : 0u;
    return load_mass(mass, i) - outgoing + incoming;
}
inline uint gravity_velocity(int i, int y, device const uint *mass,
                             device const uint *velocity,
                             device const uint *sm, device const uint *dm,
                             device const uint *active) {
    if (blocked(i, sm, dm)) return 0u;
    uint source_mass = load_mass(mass, i);
    uint outgoing = down_flow(i, y, mass, velocity, sm, dm, active);
    uint incoming = y + 1 < SIZE_Y ? down_flow(i + LAYER, y + 1, mass, velocity, sm, dm, active) : 0u;
    uint retained = source_mass - outgoing, total = retained + incoming;
    if (total == 0u) return 0u;
    int3 source_v = accelerated_velocity(velocity[i], source_mass);
    int3 incoming_v = incoming > 0u
        ? accelerated_velocity(velocity[i + LAYER], load_mass(mass, i + LAYER)) : int3(0);
    int3 mixed = (int(retained) * source_v + int(incoming) * incoming_v) / int(total);
    return pack_velocity(mixed);
}
inline int3 effective_velocity(int i, device const uint *scratch,
                               device const uint *velocity_scratch,
                               device const uint *sm, device const uint *dm) {
    if (load_mass(scratch, i) == 0u || blocked(i, sm, dm)) return int3(0);
    int3 velocity = unpack_velocity(velocity_scratch[i]);
    int y = i / LAYER;
    bool impact = y == 0;
    if (y > 0) {
        int below = i - LAYER;
        if (blocked(below, sm, dm)) impact = true;
        else if (load_mass(scratch, below) >= LATERAL_SUPPORT_MASS)
            impact = unpack_velocity(velocity_scratch[below]).y >= -IMPACT_THRESHOLD;
    }
    if (impact && velocity.y < 0)
        velocity.y = velocity.y < -IMPACT_THRESHOLD ? ((-velocity.y) * 3) / 4 : 0;
    return velocity;
}
inline bool laterally_supported(int i, device const uint *scratch,
                                device const uint *sm, device const uint *dm) {
    int y = i / LAYER;
    if (y == 0) return true;
    int below = i - LAYER;
    return blocked(below, sm, dm) || load_mass(scratch, below) >= LATERAL_SUPPORT_MASS;
}
inline int directional_speed(int3 velocity, int direction) {
    if (direction == 0) return velocity.x;
    if (direction == 1) return -velocity.x;
    if (direction == 2) return velocity.z;
    if (direction == 3) return -velocity.z;
    return velocity.y;
}
inline int surface_splash_speed(int index, device const uint *scratch,
                                device const uint *velocity_scratch,
                                device const uint *sm, device const uint *dm) {
    int fall_speed = 0;
    bool impact = false;
    for (int depth = 0; depth < 8; ++depth) {
        if (load_mass(scratch, index) == 0u || blocked(index, sm, dm)) return 0;
        int3 velocity = unpack_velocity(velocity_scratch[index]);
        fall_speed = max(fall_speed, -velocity.y);
        int y = index / LAYER;
        if (y == 0) { impact = true; break; }
        int below = index - LAYER;
        if (blocked(below, sm, dm)) { impact = true; break; }
        if (load_mass(scratch, below) < LATERAL_SUPPORT_MASS) return 0;
        if (unpack_velocity(velocity_scratch[below]).y >= -IMPACT_THRESHOLD) {
            impact = true;
            break;
        }
        index = below;
    }
    if (!impact || fall_speed <= IMPACT_THRESHOLD) return 0;
    uint hash = uint(index) * 747796405u + 2891336453u;
    if ((hash >> 30) != 0u) return 0;
    return min(VELOCITY_MAX, fall_speed * 2);
}
inline uint directional_flow(int from, int to, int direction,
                             device const uint *scratch,
                             device const uint *velocity_scratch,
                             device const uint *sm, device const uint *dm,
                             device const uint *active) {
    if (from < 0 || to < 0 || from >= SIZE_X * SIZE_Y * SIZE_Z ||
        to >= SIZE_X * SIZE_Y * SIZE_Z || blocked(from, sm, dm) || blocked(to, sm, dm) ||
        !tile_active(active, tile_for_index(from)) || !tile_active(active, tile_for_index(to))) return 0u;
    uint a = load_mass(scratch, from), b = load_mass(scratch, to);
    if (a == 0u || b >= WATER_MAX) return 0u;
    int to_y = to / LAYER;
    if (direction != 4 && to_y > 0) {
        int below_to = to - LAYER;
        uint below_mass = load_mass(scratch, below_to);
        if (surface_splash_speed(below_to, scratch, velocity_scratch, sm, dm) > 0 ||
            (below_mass < LATERAL_SUPPORT_MASS &&
             effective_velocity(below_to, scratch, velocity_scratch, sm, dm).y > 0)) return 0u;
    }
    uint base = direction != 4 && laterally_supported(from, scratch, sm, dm) && a > b
        ? (a - b) >> 3 : 0u;
    int3 effective = effective_velocity(from, scratch, velocity_scratch, sm, dm);
    int speed = directional_speed(effective, direction);
    int surface_splash = surface_splash_speed(from, scratch, velocity_scratch, sm, dm);
    bool rising_droplet = direction == 4 && surface_splash == 0 &&
                          a < LATERAL_SUPPORT_MASS && effective.y > 0;
    if (direction != 4 && (surface_splash > 0 ||
        (a < LATERAL_SUPPORT_MASS && effective.y > 0))) return 0u;
    if (direction == 4 && surface_splash == 0 && effective.y > 0 && !rising_droplet) return 0u;
    if (direction == 4) speed = max(speed, surface_splash);
    else speed += max(effective.y, surface_splash) / 2;
    speed = min(speed, VELOCITY_MAX);
    bool splash_launch = direction == 4 && surface_splash > 0;
    uint momentum = rising_droplet ? a :
        (speed > 0 ? (a * uint(speed)) / (127u * 4u) : 0u);
    uint source_cap = rising_droplet ? a : (splash_launch ? a / 4u : a / 5u);
    uint destination_cap = (rising_droplet || splash_launch)
        ? WATER_MAX - b : (WATER_MAX - b) / 5u;
    return min(base + momentum, min(source_cap, destination_cap));
}
inline int3 carried_velocity(int from, int direction, device const uint *scratch,
                             device const uint *velocity_scratch,
                             device const uint *sm, device const uint *dm) {
    int3 velocity = effective_velocity(from, scratch, velocity_scratch, sm, dm);
    if (direction == 4) {
        int splash = surface_splash_speed(from, scratch, velocity_scratch, sm, dm);
        if (splash > 0) {
            velocity.y = max(velocity.y, splash);
            uint hash = uint(from) * 747796405u + 2891336453u;
            velocity.x += (hash & 1u) != 0u ? 64 : -64;
            velocity.z += (hash & 2u) != 0u ? 64 : -64;
        }
    }
    if (direction == 0) velocity.x += HORIZONTAL_KICK;
    else if (direction == 1) velocity.x -= HORIZONTAL_KICK;
    else if (direction == 2) velocity.z += HORIZONTAL_KICK;
    else if (direction == 3) velocity.z -= HORIZONTAL_KICK;
    return clamp(velocity, int3(-VELOCITY_MAX), int3(VELOCITY_MAX));
}
inline uint directional_mass(int i, int x, int y, int z, device const uint *scratch,
                             device const uint *velocity_scratch,
                             device const uint *sm, device const uint *dm,
                             device const uint *active) {
    uint value = load_mass(scratch, i);
    if (blocked(i, sm, dm)) return value;
    int neighbors[5] = { x + 1 < SIZE_X ? i + 1 : -1, x > 0 ? i - 1 : -1,
        z + 1 < SIZE_Z ? i + SIZE_X : -1, z > 0 ? i - SIZE_X : -1,
        y + 1 < SIZE_Y ? i + LAYER : -1 };
    int reverse[4] = { 1, 0, 3, 2 };
    uint outgoing = 0u, incoming = 0u;
    for (int n = 0; n < 5; ++n) if (neighbors[n] >= 0)
        outgoing += directional_flow(i, neighbors[n], n, scratch, velocity_scratch, sm, dm, active);
    for (int n = 0; n < 4; ++n) if (neighbors[n] >= 0)
        incoming += directional_flow(neighbors[n], i, reverse[n], scratch, velocity_scratch, sm, dm, active);
    if (y > 0) incoming += directional_flow(i - LAYER, i, 4, scratch, velocity_scratch, sm, dm, active);
    return value - outgoing + incoming;
}
inline uint directional_velocity(int i, int x, int y, int z, device const uint *scratch,
                                 device const uint *velocity_scratch,
                                 device const uint *sm, device const uint *dm,
                                 device const uint *active) {
    if (blocked(i, sm, dm)) return 0u;
    uint value = load_mass(scratch, i);
    int neighbors[5] = { x + 1 < SIZE_X ? i + 1 : -1, x > 0 ? i - 1 : -1,
        z + 1 < SIZE_Z ? i + SIZE_X : -1, z > 0 ? i - SIZE_X : -1,
        y + 1 < SIZE_Y ? i + LAYER : -1 };
    int reverse[4] = { 1, 0, 3, 2 };
    uint outgoing = 0u, incoming[5] = { 0u, 0u, 0u, 0u, 0u }, incoming_sum = 0u;
    for (int n = 0; n < 5; ++n) if (neighbors[n] >= 0)
        outgoing += directional_flow(i, neighbors[n], n, scratch, velocity_scratch, sm, dm, active);
    for (int n = 0; n < 4; ++n) if (neighbors[n] >= 0) {
        incoming[n] = directional_flow(neighbors[n], i, reverse[n], scratch, velocity_scratch, sm, dm, active);
        incoming_sum += incoming[n];
    }
    int below = y > 0 ? i - LAYER : -1;
    if (below >= 0) {
        incoming[4] = directional_flow(below, i, 4, scratch, velocity_scratch, sm, dm, active);
        incoming_sum += incoming[4];
    }
    uint retained = value - outgoing, total = retained + incoming_sum;
    if (total == 0u) return 0u;
    int3 momentum = int(retained) * effective_velocity(i, scratch, velocity_scratch, sm, dm);
    for (int n = 0; n < 4; ++n) if (incoming[n] > 0u)
        momentum += int(incoming[n]) * carried_velocity(neighbors[n], reverse[n], scratch, velocity_scratch, sm, dm);
    if (incoming[4] > 0u)
        momentum += int(incoming[4]) * carried_velocity(below, 4, scratch, velocity_scratch, sm, dm);
    int3 velocity = (momentum / int(total)) * 7 / 8;
    return pack_velocity(velocity);
}

kernel void water_ca(device uint *mass [[buffer(0)]], device uint *scratch [[buffer(1)]],
                     device const uint *sm [[buffer(2)]], device const uint *dm [[buffer(3)]],
                     device const uint *active [[buffer(4)]],
                     device const uint *active_ids [[buffer(5)]],
                     device uint *velocity [[buffer(6)]],
                     device uint *velocity_scratch [[buffer(7)]],
                     constant WaterUniforms &uniforms [[buffer(8)]],
                     uint gid [[thread_position_in_grid]]) {
    int invocation = int(gid);
    if (invocation >= uniforms.cell_count) return;
    int tile_job = invocation / TILE_PAIRS;
    int pair = invocation - tile_job * TILE_PAIRS;
    int tile = int(active_ids[tile_job]);
    int ty = tile / 100, rem = tile - ty * 100;
    int tz = rem / 10, tx = rem - tz * 10;
    int x = tx * 8 + ((pair & 3) << 1);
    int z = tz * 8 + ((pair >> 2) & 7);
    int y = ty * 8 + (pair >> 5);
    int i0 = x + SIZE_X * (z + SIZE_Z * y), i1 = i0 + 1;
    uint lo, hi;
    if (uniforms.mode == 0) {
        lo = gravity_mass(i0, y, mass, velocity, sm, dm, active);
        hi = gravity_mass(i1, y, mass, velocity, sm, dm, active);
        scratch[i0 >> 1] = lo | (hi << 16);
        velocity_scratch[i0] = gravity_velocity(i0, y, mass, velocity, sm, dm, active);
        velocity_scratch[i1] = gravity_velocity(i1, y, mass, velocity, sm, dm, active);
    } else {
        lo = directional_mass(i0, x, y, z, scratch, velocity_scratch, sm, dm, active);
        hi = directional_mass(i1, x + 1, y, z, scratch, velocity_scratch, sm, dm, active);
        mass[i0 >> 1] = lo | (hi << 16);
        velocity[i0] = directional_velocity(i0, x, y, z, scratch, velocity_scratch, sm, dm, active);
        velocity[i1] = directional_velocity(i1, x + 1, y, z, scratch, velocity_scratch, sm, dm, active);
    }
}
