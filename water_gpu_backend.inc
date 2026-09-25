typedef struct WaterGpuState {
    bool ready;
    bool failed;
    int retry_attempts;
    double retry_at;
    uint32_t *static_stage;
    uint32_t *dynamic_stage;
    uint32_t *dynamic_next_stage;
    uint32_t *active_stage;
    uint32_t *active_id_stage;
    uint64_t synced_mass_generation;
    uint64_t synced_static_generation;
    uint64_t synced_dynamic_generation;
#if defined(__APPLE__) && defined(FPS_GPU_METAL)
    void *buffer[6];
#elif defined(GRAPHICS_API_OPENGL_43)
    unsigned int buffer[6];
    unsigned int program;
    int mode_location;
    int count_location;
#endif
} WaterGpuState;

static WaterGpuState waterGpu;

#define WATER_MASK_WORD_COUNT ((WATER_CELL_COUNT + 31) / 32)
#define WATER_ACTIVE_WORD_COUNT ((WATER_TILE_COUNT + 31) / 32)

static void water_gpu_release(void) {
#if defined(__APPLE__) && defined(FPS_GPU_METAL)
    for (int i = 0; i < 6; ++i) {
        if (waterGpu.buffer[i]) fps_metal_buffer_destroy(waterGpu.buffer[i]);
        waterGpu.buffer[i] = NULL;
    }
#elif defined(GRAPHICS_API_OPENGL_43)
    for (int i = 0; i < 6; ++i) {
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
        waterGpu.static_stage = (uint32_t *)malloc((size_t)WATER_MASK_WORD_COUNT * sizeof(uint32_t));
    if (!waterGpu.dynamic_stage)
        waterGpu.dynamic_stage = (uint32_t *)calloc((size_t)WATER_MASK_WORD_COUNT, sizeof(uint32_t));
    if (!waterGpu.dynamic_next_stage)
        waterGpu.dynamic_next_stage = (uint32_t *)malloc((size_t)WATER_MASK_WORD_COUNT * sizeof(uint32_t));
    if (!waterGpu.active_stage)
        waterGpu.active_stage = (uint32_t *)malloc((size_t)WATER_ACTIVE_WORD_COUNT * sizeof(uint32_t));
    if (!waterGpu.active_id_stage)
        waterGpu.active_id_stage = (uint32_t *)malloc((size_t)WATER_TILE_COUNT * sizeof(uint32_t));
    return waterGpu.static_stage && waterGpu.dynamic_stage && waterGpu.dynamic_next_stage && waterGpu.active_stage &&
           waterGpu.active_id_stage;
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
    size_t mass_bytes = (size_t)WATER_CELL_COUNT * sizeof(uint16_t);
    size_t mask_bytes = (size_t)WATER_MASK_WORD_COUNT * sizeof(uint32_t);
    waterGpu.buffer[0] = fps_metal_buffer_create(mass_bytes);
    waterGpu.buffer[1] = fps_metal_buffer_create(mass_bytes);
    waterGpu.buffer[2] = fps_metal_buffer_create(mask_bytes);
    waterGpu.buffer[3] = fps_metal_buffer_create(mask_bytes);
    waterGpu.buffer[4] = fps_metal_buffer_create((size_t)WATER_ACTIVE_WORD_COUNT * sizeof(uint32_t));
    waterGpu.buffer[5] = fps_metal_buffer_create((size_t)WATER_TILE_COUNT * sizeof(uint32_t));
    for (int i = 0; i < 6; ++i) {
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
    size_t mass_bytes = (size_t)WATER_CELL_COUNT * sizeof(uint16_t);
    size_t mask_bytes = (size_t)WATER_MASK_WORD_COUNT * sizeof(uint32_t);
    for (int i = 0; i < 6; ++i) {
        size_t bytes = i < 2 ? mass_bytes : (i < 4 ? mask_bytes :
                       (i == 4 ? (size_t)WATER_ACTIVE_WORD_COUNT * sizeof(uint32_t) :
                                 (size_t)WATER_TILE_COUNT * sizeof(uint32_t)));
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
    waterGpu.synced_mass_generation = UINT64_MAX;
    waterGpu.synced_static_generation = UINT64_MAX;
    waterGpu.synced_dynamic_generation = UINT64_MAX;
    return true;
}

static bool water_gpu_run_batch(int fixed_steps) {
    if (fixed_steps <= 0 || !water_gpu_initialize_backend()) return false;
    bool upload_mass = waterGpu.synced_mass_generation != waterSystem.mass_generation;
    bool upload_static = waterGpu.synced_static_generation != waterSystem.static_generation;
    bool upload_dynamic = waterGpu.synced_dynamic_generation != waterSystem.dynamic_generation;
    if (upload_static) {
        memset(waterGpu.static_stage, 0, (size_t)WATER_MASK_WORD_COUNT * sizeof(uint32_t));
        for (int i = 0; i < WATER_CELL_COUNT; ++i)
            if (waterSystem.static_solid[i]) waterGpu.static_stage[i >> 5] |= 1u << (i & 31);
    }
    int dynamic_changed_words = 0;
    bool dynamic_full_upload = waterGpu.synced_dynamic_generation == UINT64_MAX;
    if (upload_dynamic) {
        memset(waterGpu.dynamic_next_stage, 0, (size_t)WATER_MASK_WORD_COUNT * sizeof(uint32_t));
        for (uint32_t n = 0; n < waterSystem.dynamic_cell_count; ++n) {
            int i = (int)waterSystem.dynamic_cells[n];
            waterGpu.dynamic_next_stage[i >> 5] |= 1u << (i & 31);
        }
        for (int w = 0; w < WATER_MASK_WORD_COUNT; ++w)
            dynamic_changed_words += waterGpu.dynamic_stage[w] != waterGpu.dynamic_next_stage[w];
        if (dynamic_changed_words > WATER_MASK_WORD_COUNT / 4) dynamic_full_upload = true;
    }
    (void)dynamic_changed_words;
    (void)dynamic_full_upload;
    memset(waterGpu.active_stage, 0, (size_t)WATER_ACTIVE_WORD_COUNT * sizeof(uint32_t));
    for (int i = 0; i < WATER_TILE_COUNT; ++i)
        if (waterSystem.batch_tile_active[i]) waterGpu.active_stage[i >> 5] |= 1u << (i & 31);
    for (int i = 0; i < waterSystem.active_tile_count; ++i)
        waterGpu.active_id_stage[i] = waterSystem.active_tile_ids[i];
    bool ok = true;
    double upload_started = pbd_time_now_ms();
#if defined(__APPLE__) && defined(FPS_GPU_METAL)
    size_t mass_bytes = (size_t)WATER_CELL_COUNT * sizeof(uint16_t);
    size_t mask_bytes = (size_t)WATER_MASK_WORD_COUNT * sizeof(uint32_t);
    if (upload_mass) ok = fps_metal_buffer_update(waterGpu.buffer[0], waterSystem.mass, mass_bytes, 0) &&
                          fps_metal_buffer_update(waterGpu.buffer[1], waterSystem.scratch, mass_bytes, 0);
    if (ok && upload_static) ok = fps_metal_buffer_update(waterGpu.buffer[2], waterGpu.static_stage, mask_bytes, 0);
    if (ok && upload_dynamic && dynamic_full_upload)
        ok = fps_metal_buffer_update(waterGpu.buffer[3], waterGpu.dynamic_next_stage, mask_bytes, 0);
    if (ok && upload_dynamic && !dynamic_full_upload) {
        for (int first = 0; first < WATER_MASK_WORD_COUNT && ok; ) {
            while (first < WATER_MASK_WORD_COUNT &&
                   waterGpu.dynamic_stage[first] == waterGpu.dynamic_next_stage[first]) first++;
            int last = first;
            while (last < WATER_MASK_WORD_COUNT &&
                   waterGpu.dynamic_stage[last] != waterGpu.dynamic_next_stage[last]) last++;
            if (last > first) ok = fps_metal_buffer_update(waterGpu.buffer[3],
                waterGpu.dynamic_next_stage + first, (size_t)(last - first) * sizeof(uint32_t),
                (size_t)first * sizeof(uint32_t));
            first = last;
        }
    }
    if (ok) ok = fps_metal_buffer_update(waterGpu.buffer[4], waterGpu.active_stage,
                                         (size_t)WATER_ACTIVE_WORD_COUNT * sizeof(uint32_t), 0);
    if (ok) ok = fps_metal_buffer_update(waterGpu.buffer[5], waterGpu.active_id_stage,
                                         (size_t)waterSystem.active_tile_count * sizeof(uint32_t), 0);
    waterSystem.diagnostics.last_gpu_upload_ms = pbd_time_now_ms() - upload_started;
    double dispatch_started = pbd_time_now_ms();
    if (ok) ok = fps_metal_water_dispatch(waterGpu.buffer[0], waterGpu.buffer[1],
                                          waterGpu.buffer[2], waterGpu.buffer[3],
                                          waterGpu.buffer[4], waterGpu.buffer[5], fixed_steps,
                                          waterSystem.active_tile_count * 256);
    waterSystem.diagnostics.last_gpu_dispatch_ms = pbd_time_now_ms() - dispatch_started;
    double readback_started = pbd_time_now_ms();
    if (ok) ok = fps_metal_buffer_read(waterGpu.buffer[0], waterSystem.gpu_readback, mass_bytes, 0);
    waterSystem.diagnostics.last_gpu_readback_ms = pbd_time_now_ms() - readback_started;
    if (ok) waterSystem.diagnostics.backend = "gpu-metal";
#elif defined(GRAPHICS_API_OPENGL_43)
    size_t mass_bytes = (size_t)WATER_CELL_COUNT * sizeof(uint16_t);
    size_t mask_bytes = (size_t)WATER_MASK_WORD_COUNT * sizeof(uint32_t);
    if (upload_mass) {
        rlUpdateShaderBuffer(waterGpu.buffer[0], waterSystem.mass, (unsigned int)mass_bytes, 0);
        rlUpdateShaderBuffer(waterGpu.buffer[1], waterSystem.scratch, (unsigned int)mass_bytes, 0);
    }
    if (upload_static) rlUpdateShaderBuffer(waterGpu.buffer[2], waterGpu.static_stage, (unsigned int)mask_bytes, 0);
    if (upload_dynamic && dynamic_full_upload)
        rlUpdateShaderBuffer(waterGpu.buffer[3], waterGpu.dynamic_next_stage, (unsigned int)mask_bytes, 0);
    if (upload_dynamic && !dynamic_full_upload) {
        for (int first = 0; first < WATER_MASK_WORD_COUNT; ) {
            while (first < WATER_MASK_WORD_COUNT &&
                   waterGpu.dynamic_stage[first] == waterGpu.dynamic_next_stage[first]) first++;
            int last = first;
            while (last < WATER_MASK_WORD_COUNT &&
                   waterGpu.dynamic_stage[last] != waterGpu.dynamic_next_stage[last]) last++;
            if (last > first) rlUpdateShaderBuffer(waterGpu.buffer[3], waterGpu.dynamic_next_stage + first,
                (unsigned int)((last - first) * sizeof(uint32_t)),
                (unsigned int)(first * sizeof(uint32_t)));
            first = last;
        }
    }
    rlUpdateShaderBuffer(waterGpu.buffer[4], waterGpu.active_stage,
                         (unsigned int)(WATER_ACTIVE_WORD_COUNT * sizeof(uint32_t)), 0);
    rlUpdateShaderBuffer(waterGpu.buffer[5], waterGpu.active_id_stage,
                         (unsigned int)(waterSystem.active_tile_count * sizeof(uint32_t)), 0);
    waterSystem.diagnostics.last_gpu_upload_ms = pbd_time_now_ms() - upload_started;
    while (glGetError() != GL_NO_ERROR) { }
    glUseProgram(waterGpu.program);
    glUniform1i(waterGpu.count_location, waterSystem.active_tile_count * 256);
    for (int i = 0; i < 6; ++i) rlBindShaderBuffer(waterGpu.buffer[i], (unsigned int)i);
    GLuint groups = (GLuint)(((int)waterSystem.active_tile_count * 256 + 127) / 128);
    double dispatch_started = pbd_time_now_ms();
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
    waterSystem.diagnostics.last_gpu_dispatch_ms = pbd_time_now_ms() - dispatch_started;
    double readback_started = pbd_time_now_ms();
    if (ok) rlReadShaderBuffer(waterGpu.buffer[0], waterSystem.gpu_readback, (unsigned int)mass_bytes, 0);
    waterSystem.diagnostics.last_gpu_readback_ms = pbd_time_now_ms() - readback_started;
    if (ok) waterSystem.diagnostics.backend = "gpu-gl43";
#else
    ok = false;
#endif
    int invalid_mass_index = -1;
    if (ok) {
        uint64_t before_mass = 0, after_mass = 0;
        for (int i = 0; i < WATER_CELL_COUNT; ++i) {
            if (waterSystem.gpu_readback[i] > WATER_MAX_MASS) {
                invalid_mass_index = i;
                ok = false;
                break;
            }
            before_mass += waterSystem.batch_snapshot[i];
            after_mass += waterSystem.gpu_readback[i];
        }
        if (ok && before_mass != after_mass) ok = false;
    }
    if (!ok) {
        char reason[160];
        if (invalid_mass_index >= 0) {
            snprintf(reason, sizeof(reason), "mass validation failed at %d: %u",
                     invalid_mass_index, (unsigned)waterSystem.gpu_readback[invalid_mass_index]);
            water_gpu_fail(reason);
        } else water_gpu_fail("dispatch or readback failed");
        return false;
    }
    waterGpu.failed = false;
    waterGpu.retry_attempts = 0;
    memcpy(waterSystem.mass, waterSystem.gpu_readback, sizeof(waterSystem.mass));
    if (++waterSystem.mass_generation == 0) waterSystem.mass_generation = 1;
    waterGpu.synced_mass_generation = waterSystem.mass_generation;
    waterGpu.synced_static_generation = waterSystem.static_generation;
    waterGpu.synced_dynamic_generation = waterSystem.dynamic_generation;
    if (upload_dynamic) memcpy(waterGpu.dynamic_stage, waterGpu.dynamic_next_stage,
                               (size_t)WATER_MASK_WORD_COUNT * sizeof(uint32_t));
    waterSystem.cpu_scratch_valid = false;
    return true;
}

static void water_gpu_shutdown_backend(void) {
    water_gpu_release();
    free(waterGpu.static_stage);
    free(waterGpu.dynamic_stage);
    free(waterGpu.dynamic_next_stage);
    free(waterGpu.active_stage);
    free(waterGpu.active_id_stage);
    memset(&waterGpu, 0, sizeof(waterGpu));
}
