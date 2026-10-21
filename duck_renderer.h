#ifndef DUCK_RENDERER_H
#define DUCK_RENDERER_H
#include "raylib.h"

/* Cosmetic mallard. Position is the floor/waterline anchor; forward is -Z.
 * Caller owns movement and samples time once per frame for all viewports. */
typedef struct DuckPose {
    Vector3 position;
    float yaw;              /* degrees */
    float scale;            /* 1 = 1.55 units bill-to-tail */
    float time;
    float speed;            /* world units/second; zero gives idle breathing */
    float phase;            /* radians, optional per-duck variation */
} DuckPose;

/* NULL selects production shader files. A missing directory tests fallback. */
bool duck_renderer_init(const char *shaderDirectory);
void duck_renderer_shutdown(void);
/* Opaque color pass. No simulation, allocation, or time queries. */
void duck_draw(DuckPose pose);
/* Optional contact shadow and restrained wake. Draw on an exposed flat floor
 * after opaque geometry, before duck_draw; does not write depth. */
void duck_draw_floor(DuckPose pose, bool wake);
int duck_renderer_triangle_count(void);
bool duck_renderer_has_shader(void);
#endif
