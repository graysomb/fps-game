typedef struct WaterGpuState {
    bool ready;
    bool failed;
    int retry_attempts;
    double retry_at;
    uint32_t *static_stage;
    uint32_t *dynamic_stage;
    uint32_t *active_stage;
#if defined(__APPLE__) && defined(FPS_GPU_METAL)
    void *buffer[5];
#elif defined(GRAPHICS_API_OPENGL_43)
    unsigned int buffer[5];
    unsigned int program;
    int mode_location;
    int count_location;
#endif
} WaterGpuState;

static WaterGpuState waterGpu;

static void water_gpu_release(void) {
#if defined(__APPLE__) && defined(FPS_GPU_METAL)
    for (int i = 0; i < 5; ++i) {
        if (waterGpu.buffer[i]) fps_metal_buffer_destroy(waterGpu.buffer[i]);
        waterGpu.buffer[i] = NULL;
    }
#elif defined(GRAPHICS_API_OPENGL_43)
    for (int i = 0; i < 5; ++i) {
        if (waterGpu.buffer[i]) rlUnloadShaderBuffer(waterGpu.buffer[i]);
        waterGpu.buffer[i] = 0;
    }
    if (waterGpu.program) glDeleteProgram(waterGpu.program);
    waterGpu.program = 0;
#endif
    waterGpu.ready = false;
}

static void water_gpu_fail(const char *reason) {
    if (physicsReportRequested) {
        fprintf(stderr, "water GPU fallback: %s\n", reason ? reason : "unknown failure");
    }
    water_gpu_release();
    waterGpu.failed = true;
    static const double retry_delay[] = { 5.0, 15.0, 60.0 };
    int index = waterGpu.retry_attempts;
    if (index > 2) index = 2;
    waterGpu.retry_at = GetTime() + retry_delay[index];
    waterGpu.retry_attempts++;
}

static bool water_gpu_allocate_staging(void) {
    if (!waterGpu.static_stage)
        waterGpu.static_stage = (uint32_t *)malloc((size_t)WATER_CELL_COUNT * sizeof(uint32_t));
    if (!waterGpu.dynamic_stage)
        waterGpu.dynamic_stage = (uint32_t *)malloc((size_t)WATER_CELL_COUNT * sizeof(uint32_t));
    if (!waterGpu.active_stage)
        waterGpu.active_stage = (uint32_t *)malloc((size_t)WATER_TILE_COUNT * sizeof(uint32_t));
    return waterGpu.static_stage && waterGpu.dynamic_stage && waterGpu.active_stage;
}

static bool water_gpu_initialize_backend(void) {
    if (waterGpu.ready) return true;
    if (waterGpu.failed && GetTime() < waterGpu.retry_at) return false;
    waterGpu.failed = false;
    if (!water_gpu_allocate_staging()) {
        water_gpu_fail("staging allocation failed");
        return false;
    }
#if defined(__APPLE__) && defined(FPS_GPU_METAL)
    if (!gpuPhysics.ready || !fps_metal_water_available()) {
        water_gpu_fail("water_ca is unavailable in the Metal library");
        return false;
    }
    size_t cell_bytes = (size_t)WATER_CELL_COUNT * sizeof(uint32_t);
    size_t tile_bytes = (size_t)WATER_TILE_COUNT * sizeof(uint32_t);
    waterGpu.buffer[0] = fps_metal_buffer_create(cell_bytes);
    waterGpu.buffer[1] = fps_metal_buffer_create(cell_bytes);
    waterGpu.buffer[2] = fps_metal_buffer_create(cell_bytes);
    waterGpu.buffer[3] = fps_metal_buffer_create(cell_bytes);
    waterGpu.buffer[4] = fps_metal_buffer_create(tile_bytes);
    for (int i = 0; i < 5; ++i) {
        if (!waterGpu.buffer[i]) {
            water_gpu_fail("Metal water buffer allocation failed");
            return false;
        }
    }
#elif defined(GRAPHICS_API_OPENGL_43)
    if (!gpuPhysics.ready || rlGetVersion() != RL_OPENGL_43) {
        water_gpu_fail("OpenGL 4.3 compute is unavailable");
        return false;
    }
    char *source = LoadFileText("shaders/water/water_ca.comp");
    if (!source) {
        water_gpu_fail("missing shaders/water/water_ca.comp");
        return false;
    }
    unsigned int shader = glCreateShader(GL_COMPUTE_SHADER);
    const char *shader_source = source;
    glShaderSource(shader, 1, &shader_source, NULL);
    glCompileShader(shader);
    UnloadFileText(source);
    GLint compiled = GL_FALSE;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
    if (compiled == GL_FALSE) {
        char log[1024] = { 0 };
        glGetShaderInfoLog(shader, (GLsizei)sizeof(log) - 1, NULL, log);
        glDeleteShader(shader);
        water_gpu_fail(log);
        return false;
    }
    waterGpu.program = glCreateProgram();
    glAttachShader(waterGpu.program, shader);
    glLinkProgram(waterGpu.program);
    glDeleteShader(shader);
    GLint linked = GL_FALSE;
    glGetProgramiv(waterGpu.program, GL_LINK_STATUS, &linked);
    if (linked == GL_FALSE) {
        char log[1024] = { 0 };
        glGetProgramInfoLog(waterGpu.program, (GLsizei)sizeof(log) - 1, NULL, log);
        water_gpu_fail(log);
        return false;
    }
    waterGpu.mode_location = glGetUniformLocation(waterGpu.program, "uMode");
    waterGpu.count_location = glGetUniformLocation(waterGpu.program, "uCellCount");
    if (waterGpu.mode_location < 0 || waterGpu.count_location < 0) {
        water_gpu_fail("water compute uniforms are missing");
        return false;
    }
    size_t cell_bytes = (size_t)WATER_CELL_COUNT * sizeof(uint32_t);
    size_t tile_bytes = (size_t)WATER_TILE_COUNT * sizeof(uint32_t);
    for (int i = 0; i < 5; ++i) {
        size_t bytes = i == 4 ? tile_bytes : cell_bytes;
        waterGpu.buffer[i] = rlLoadShaderBuffer((unsigned int)bytes, NULL, RL_DYNAMIC_COPY);
        if (!waterGpu.buffer[i]) {
            water_gpu_fail("OpenGL water buffer allocation failed");
            return false;
        }
    }
#else
    water_gpu_fail("this build has no native water compute backend");
    return false;
#endif
    waterGpu.ready = true;
    return true;
}

static bool water_gpu_run_batch(int fixed_steps) {
    if (fixed_steps <= 0 || !water_gpu_initialize_backend()) return false;
    for (int i = 0; i < WATER_CELL_COUNT; ++i) {
        waterGpu.static_stage[i] = waterSystem.static_solid[i] ? 1u : 0u;
        waterGpu.dynamic_stage[i] = waterSystem.dynamic_solid[i] ? 1u : 0u;
    }
    for (int i = 0; i < WATER_TILE_COUNT; ++i)
        waterGpu.active_stage[i] = waterSystem.tile_active[i] ? 1u : 0u;
    bool ok = true;
#if defined(__APPLE__) && defined(FPS_GPU_METAL)
    size_t cell_bytes = (size_t)WATER_CELL_COUNT * sizeof(uint32_t);
    size_t tile_bytes = (size_t)WATER_TILE_COUNT * sizeof(uint32_t);
    ok = fps_metal_buffer_update(waterGpu.buffer[0], waterSystem.mass, cell_bytes, 0) &&
         fps_metal_buffer_update(waterGpu.buffer[1], waterSystem.scratch, cell_bytes, 0) &&
         fps_metal_buffer_update(waterGpu.buffer[2], waterGpu.static_stage, cell_bytes, 0) &&
         fps_metal_buffer_update(waterGpu.buffer[3], waterGpu.dynamic_stage, cell_bytes, 0) &&
         fps_metal_buffer_update(waterGpu.buffer[4], waterGpu.active_stage, tile_bytes, 0);
    if (ok) ok = fps_metal_water_dispatch(waterGpu.buffer[0], waterGpu.buffer[1],
                                          waterGpu.buffer[2], waterGpu.buffer[3],
                                          waterGpu.buffer[4], fixed_steps, WATER_CELL_COUNT);
    if (ok) ok = fps_metal_buffer_read(waterGpu.buffer[0], waterSystem.mass, cell_bytes, 0);
    if (ok) waterSystem.diagnostics.backend = "gpu-metal";
#elif defined(GRAPHICS_API_OPENGL_43)
    size_t cell_bytes = (size_t)WATER_CELL_COUNT * sizeof(uint32_t);
    size_t tile_bytes = (size_t)WATER_TILE_COUNT * sizeof(uint32_t);
    rlUpdateShaderBuffer(waterGpu.buffer[0], waterSystem.mass, (unsigned int)cell_bytes, 0);
    rlUpdateShaderBuffer(waterGpu.buffer[1], waterSystem.scratch, (unsigned int)cell_bytes, 0);
    rlUpdateShaderBuffer(waterGpu.buffer[2], waterGpu.static_stage, (unsigned int)cell_bytes, 0);
    rlUpdateShaderBuffer(waterGpu.buffer[3], waterGpu.dynamic_stage, (unsigned int)cell_bytes, 0);
    rlUpdateShaderBuffer(waterGpu.buffer[4], waterGpu.active_stage, (unsigned int)tile_bytes, 0);
    while (glGetError() != GL_NO_ERROR) { }
    glUseProgram(waterGpu.program);
    glUniform1i(waterGpu.count_location, WATER_CELL_COUNT);
    for (int i = 0; i < 5; ++i) rlBindShaderBuffer(waterGpu.buffer[i], (unsigned int)i);
    GLuint groups = (GLuint)((WATER_CELL_COUNT + 127) / 128);
    for (int step = 0; step < fixed_steps; ++step) {
        glUniform1i(waterGpu.mode_location, 0);
        glDispatchCompute(groups, 1, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
        glUniform1i(waterGpu.mode_location, 1);
        glDispatchCompute(groups, 1, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    }
    glMemoryBarrier(GL_BUFFER_UPDATE_BARRIER_BIT | GL_SHADER_STORAGE_BARRIER_BIT);
    ok = glGetError() == GL_NO_ERROR;
    if (ok) rlReadShaderBuffer(waterGpu.buffer[0], waterSystem.mass, (unsigned int)cell_bytes, 0);
    if (ok) waterSystem.diagnostics.backend = "gpu-gl43";
#else
    ok = false;
#endif
    if (ok) {
        for (int i = 0; i < WATER_CELL_COUNT; ++i) {
            if (waterSystem.mass[i] > WATER_MAX_MASS) { ok = false; break; }
        }
    }
    if (!ok) {
        water_gpu_fail("dispatch, readback, or mass validation failed");
        return false;
    }
    waterGpu.failed = false;
    waterGpu.retry_attempts = 0;
    return true;
}

static void water_gpu_shutdown_backend(void) {
    water_gpu_release();
    free(waterGpu.static_stage);
    free(waterGpu.dynamic_stage);
    free(waterGpu.active_stage);
    memset(&waterGpu, 0, sizeof(waterGpu));
}
