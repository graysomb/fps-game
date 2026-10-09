#include <metal_stdlib>
using namespace metal;

struct WaterUniforms {
    int mode;
    int cell_count;
    int padding0;
    int padding1;
};

constant uint WATER_MAX = 65535u;
constant int SIZE_X = 80;
constant int SIZE_Y = 80;
constant int SIZE_Z = 80;
constant int LAYER = SIZE_X * SIZE_Z;

inline bool water_blocked(int i, device const uint *static_mask,
                          device const uint *dynamic_mask)
{
    return static_mask[i] != 0u || dynamic_mask[i] != 0u;
}

inline uint water_down_flow(int i, int y, device const uint *mass,
                            device const uint *static_mask,
                            device const uint *dynamic_mask)
{
    if (water_blocked(i, static_mask, dynamic_mask) || y == 0) return 0u;
    int below = i - LAYER;
    if (water_blocked(below, static_mask, dynamic_mask)) return 0u;
    return min(mass[i], WATER_MAX - mass[below]);
}

inline uint water_horizontal_flow(int from, int to, int count,
                                  device const uint *scratch,
                                  device const uint *static_mask,
                                  device const uint *dynamic_mask)
{
    if (to < 0 || to >= count ||
        water_blocked(from, static_mask, dynamic_mask) ||
        water_blocked(to, static_mask, dynamic_mask)) return 0u;
    return scratch[from] > scratch[to] ? (scratch[from] - scratch[to]) >> 3 : 0u;
}

kernel void water_ca(device uint *mass [[buffer(0)]],
                     device uint *scratch [[buffer(1)]],
                     device const uint *static_mask [[buffer(2)]],
                     device const uint *dynamic_mask [[buffer(3)]],
                     device const uint *active_tile [[buffer(4)]],
                     constant WaterUniforms &uniforms [[buffer(5)]],
                     uint gid [[thread_position_in_grid]])
{
    int i = int(gid);
    if (i >= uniforms.cell_count) return;
    int y = i / LAYER;
    int rem = i - y * LAYER;
    int z = rem / SIZE_X;
    int x = rem - z * SIZE_X;
    int tile = (x / 8) + 10 * ((z / 8) + 10 * (y / 8));
    if (uniforms.mode == 0) {
        if (active_tile[tile] == 0u || water_blocked(i, static_mask, dynamic_mask)) {
            scratch[i] = mass[i];
            return;
        }
        uint outgoing = water_down_flow(i, y, mass, static_mask, dynamic_mask);
        uint incoming = y + 1 < SIZE_Y
            ? water_down_flow(i + LAYER, y + 1, mass, static_mask, dynamic_mask) : 0u;
        scratch[i] = mass[i] - outgoing + incoming;
        return;
    }
    if (active_tile[tile] == 0u || water_blocked(i, static_mask, dynamic_mask)) {
        mass[i] = scratch[i];
        return;
    }
    int neighbors[4] = {
        x + 1 < SIZE_X ? i + 1 : -1,
        x > 0 ? i - 1 : -1,
        z + 1 < SIZE_Z ? i + SIZE_X : -1,
        z > 0 ? i - SIZE_X : -1
    };
    uint outgoing = 0u;
    uint incoming = 0u;
    for (int n = 0; n < 4; ++n) {
        if (neighbors[n] < 0) continue;
        outgoing += water_horizontal_flow(i, neighbors[n], uniforms.cell_count,
                                           scratch, static_mask, dynamic_mask);
        incoming += water_horizontal_flow(neighbors[n], i, uniforms.cell_count,
                                           scratch, static_mask, dynamic_mask);
    }
    mass[i] = scratch[i] - outgoing + incoming;
}
