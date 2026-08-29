#ifndef FPS_PHYSICS_BACKEND_H
#define FPS_PHYSICS_BACKEND_H

#include <stdbool.h>

typedef enum PhysicsBackendKind {
    PHYSICS_BACKEND_AUTO = 0,
    PHYSICS_BACKEND_GPU_GL43,
    PHYSICS_BACKEND_GPU_METAL,
    PHYSICS_BACKEND_CPU_MT,
    PHYSICS_BACKEND_CPU_ST
} PhysicsBackendKind;

typedef struct PhysicsBackendStatus {
    PhysicsBackendKind requested;
    PhysicsBackendKind active;
    bool initialized;
    bool gpu_compiled;
    bool gpu_available;
    bool sticky_fallback;
    int cpu_workers;
    double last_step_ms;
    char fallback_reason[192];
} PhysicsBackendStatus;

const char *physics_backend_name(PhysicsBackendKind kind);
bool physics_backend_parse(const char *value, PhysicsBackendKind *out);

#endif
