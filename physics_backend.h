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

typedef enum PhysicsGpuFailureKind {
    PHYSICS_GPU_FAILURE_NONE = 0,
    PHYSICS_GPU_FAILURE_TRANSIENT,
    PHYSICS_GPU_FAILURE_VALIDATION,
    PHYSICS_GPU_FAILURE_PERMANENT
} PhysicsGpuFailureKind;

typedef struct PhysicsBackendStatus {
    PhysicsBackendKind requested;
    PhysicsBackendKind active;
    bool initialized;
    bool gpu_compiled;
    bool gpu_available;
    bool sticky_fallback;
    bool gpu_recovery_pending;
    bool gpu_recovery_probe;
    int gpu_retry_attempts;
    int gpu_retry_limit;
    double gpu_retry_at;
    PhysicsGpuFailureKind last_gpu_failure;
    int cpu_workers;
    double last_step_ms;
    char fallback_reason[192];
} PhysicsBackendStatus;

const char *physics_backend_name(PhysicsBackendKind kind);
bool physics_backend_parse(const char *value, PhysicsBackendKind *out);

#endif
