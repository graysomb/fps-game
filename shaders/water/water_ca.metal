#include <metal_stdlib>
using namespace metal;

struct WaterUniforms { int mode; int cell_count; int padding0; int padding1; };
constant uint WATER_MAX = 65535u;
constant uint LATERAL_SUPPORT_MASS = 60000u;
constant int SIZE_X = 80, SIZE_Y = 80, SIZE_Z = 80;
constant int LAYER = SIZE_X * SIZE_Z;
constant int TILE_PAIRS = 256;

inline uint load_mass(device const uint *packed, int i) {
    return (packed[i >> 1] >> ((i & 1) * 16)) & WATER_MAX;
}
inline bool mask_bit(device const uint *mask, int i) {
    return ((mask[i >> 5] >> (i & 31)) & 1u) != 0u;
}
inline bool water_blocked(int i, device const uint *sm, device const uint *dm) {
    return mask_bit(sm, i) || mask_bit(dm, i);
}
inline int water_tile_for_index(int i) {
    int y = i / LAYER, rem = i - y * LAYER;
    int z = rem / SIZE_X, x = rem - z * SIZE_X;
    return (x / 8) + 10 * ((z / 8) + 10 * (y / 8));
}
inline bool tile_active(device const uint *active, int tile) {
    return ((active[tile >> 5] >> (tile & 31)) & 1u) != 0u;
}
inline uint down_flow(int i, int y, device const uint *mass,
                      device const uint *sm, device const uint *dm,
                      device const uint *active) {
    if (water_blocked(i, sm, dm) || y == 0) return 0u;
    int below = i - LAYER;
    if (water_blocked(below, sm, dm) || !tile_active(active, water_tile_for_index(i)) ||
        !tile_active(active, water_tile_for_index(below))) return 0u;
    return min(load_mass(mass, i), WATER_MAX - load_mass(mass, below));
}
inline uint horizontal_flow(int from, int to, device const uint *scratch,
                            device const uint *sm, device const uint *dm,
                            device const uint *active) {
    if (to < 0 || to >= SIZE_X * SIZE_Y * SIZE_Z || water_blocked(from, sm, dm) ||
        water_blocked(to, sm, dm) || !tile_active(active, water_tile_for_index(from)) ||
        !tile_active(active, water_tile_for_index(to))) return 0u;
    int from_y = from / LAYER;
    if (from_y > 0) {
        int below = from - LAYER;
        if (!water_blocked(below, sm, dm) && load_mass(scratch, below) < LATERAL_SUPPORT_MASS)
            return 0u;
    }
    uint a = load_mass(scratch, from), b = load_mass(scratch, to);
    return a > b ? (a - b) >> 3 : 0u;
}
inline uint gravity_value(int i, int y, device const uint *mass,
                          device const uint *sm, device const uint *dm,
                          device const uint *active) {
    if (water_blocked(i, sm, dm)) return load_mass(mass, i);
    uint outgoing = down_flow(i, y, mass, sm, dm, active);
    uint incoming = y + 1 < SIZE_Y ? down_flow(i + LAYER, y + 1, mass, sm, dm, active) : 0u;
    return load_mass(mass, i) - outgoing + incoming;
}
inline uint horizontal_value(int i, int x, int z, device const uint *scratch,
                             device const uint *sm, device const uint *dm,
                             device const uint *active) {
    uint value = load_mass(scratch, i);
    if (water_blocked(i, sm, dm)) return value;
    int neighbors[4] = { x + 1 < SIZE_X ? i + 1 : -1, x > 0 ? i - 1 : -1,
                         z + 1 < SIZE_Z ? i + SIZE_X : -1, z > 0 ? i - SIZE_X : -1 };
    uint outgoing = 0u, incoming = 0u;
    for (int n = 0; n < 4; ++n) {
        if (neighbors[n] < 0) continue;
        outgoing += horizontal_flow(i, neighbors[n], scratch, sm, dm, active);
        incoming += horizontal_flow(neighbors[n], i, scratch, sm, dm, active);
    }
    return value - outgoing + incoming;
}

kernel void water_ca(device uint *mass [[buffer(0)]], device uint *scratch [[buffer(1)]],
                     device const uint *sm [[buffer(2)]], device const uint *dm [[buffer(3)]],
                     device const uint *active [[buffer(4)]],
                     device const uint *active_ids [[buffer(5)]],
                     constant WaterUniforms &uniforms [[buffer(6)]],
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
        lo = gravity_value(i0, y, mass, sm, dm, active);
        hi = gravity_value(i1, y, mass, sm, dm, active);
        scratch[i0 >> 1] = lo | (hi << 16);
    } else {
        lo = horizontal_value(i0, x, z, scratch, sm, dm, active);
        hi = horizontal_value(i1, x + 1, z, scratch, sm, dm, active);
        mass[i0 >> 1] = lo | (hi << 16);
    }
}
