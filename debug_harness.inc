#include <errno.h>

#define DEBUG_PATH_CAPACITY 2048
#define DEBUG_MAX_CAPTURE_STEPS 64
#define DEBUG_MAX_SAMPLES 96
#define DEBUG_MAX_CAPTURES 96
#define DEBUG_MAX_FAILURES 32
#define DEBUG_MAX_TRACKED_CELLS 1024
#define DEBUG_MESSAGE_CAPACITY 256
#define DEBUG_TEST_TAG 63
#define DEBUG_SUPPORT_TAG 62
#define DEBUG_THROW_TAG 61
#define DEBUG_EXPECTED_VOXELS 18
#define DEBUG_MIN_STEPS 60
#define DEBUG_DEFAULT_STEPS 60
#define DEBUG_RENDER_WIDTH 1280
#define DEBUG_RENDER_HEIGHT 720

typedef struct DebugRunConfig {
    bool requested;
    bool show_window;
    bool capture_steps_explicit;
    int steps;
    int capture_steps[DEBUG_MAX_CAPTURE_STEPS];
    int capture_step_count;
    char scenario[64];
    char output_dir[DEBUG_PATH_CAPACITY];
} DebugRunConfig;

typedef struct DebugSample {
    char phase[48];
    int step;
    int tagged_voxels;
    int tagged_dynamic;
    int tagged_static;
    int sim_particles;
    int zero_mass_corners;
    int invalid_particle_refs;
    int static_glued_faces;
    float centroid_y;
    float centroid_drop;
    float min_y;
    float max_y;
    float mean_velocity_y;
    uint64_t static_hash_generation;
    uint64_t gpu_static_uploaded_generation;
    double physics_step_ms;
} DebugSample;

typedef struct DebugTrackedCell {
    int x;
    int y;
    int z;
} DebugTrackedCell;

typedef struct DebugRunState {
    RenderTexture2D target;
    Camera3D camera;
    int expected_voxels;
    bool activation_succeeded;
    bool sleeping_conversion_succeeded;
    bool support_removal_succeeded;
    bool drop_sampled;
    bool delayed_activation;
    bool overhang_impact;
    bool overhang_anchor_observed;
    bool overhang_deflection_sampled;
    bool gpu_readback_performed;
    int activation_step;
    float initial_centroid_y;
    float drop_at_step_60;
    float overhang_initial_root_y;
    float overhang_initial_tip_y;
    float overhang_root_drop_at_step_30;
    float overhang_tip_drop_at_step_30;
    float overhang_relative_deflection_at_step_30;
    int gpu_static_ssbo_diff_entries;
    int gpu_static_ssbo_stale_cells;
    int gpu_sim_id_mismatches;
    int gpu_tagged_corners_missing_from_sim_ids;
    int gpu_zero_mass_corners;
    DebugTrackedCell tracked_cells[DEBUG_MAX_TRACKED_CELLS];
    int tracked_cell_count;
    DebugSample samples[DEBUG_MAX_SAMPLES];
    int sample_count;
    char captures[DEBUG_MAX_CAPTURES][DEBUG_PATH_CAPACITY];
    int capture_count;
    char failures[DEBUG_MAX_FAILURES][DEBUG_MESSAGE_CAPACITY];
    int failure_count;
    char output_dir[DEBUG_PATH_CAPACITY];
    char report_path[DEBUG_PATH_CAPACITY];
} DebugRunState;

typedef bool (*DebugScenarioSetupFn)(DebugRunState *state);

typedef struct DebugScenarioDefinition {
    const char *name;
    const char *description;
    DebugScenarioSetupFn setup;
} DebugScenarioDefinition;

static DebugRunConfig debugRunConfig = {
    .steps = DEBUG_DEFAULT_STEPS
};

static bool debug_setup_dynamic_freefall(DebugRunState *state);
static bool debug_setup_activation_floating(DebugRunState *state);
static bool debug_setup_sleep_wake_floating(DebugRunState *state);
static bool debug_setup_tether_throw_floating(DebugRunState *state);
static bool debug_setup_overhang_impact(DebugRunState *state);

static const DebugScenarioDefinition debugScenarios[] = {
    { "dynamic-freefall", "Dynamic cluster freefall control", debug_setup_dynamic_freefall },
    { "activation-floating", "Static cluster activated while unsupported", debug_setup_activation_floating },
    { "sleep-wake-floating", "Supported cluster sleeps, loses support, then activates", debug_setup_sleep_wake_floating },
    { "tether-throw-floating", "Released tether voxel wakes an unsupported sleeping cluster", debug_setup_tether_throw_floating },
    { "overhang-impact", "Falling dynamic voxel activates and bends a supported static overhang", debug_setup_overhang_impact }
};

static const int debugScenarioCount = (int)(sizeof(debugScenarios) / sizeof(debugScenarios[0]));

static bool debug_copy_option(char *out, size_t capacity, const char *value, const char *label)
{
    if (!out || capacity == 0 || !value || !value[0]) {
        fprintf(stderr, "Missing value for %s\n", label);
        return false;
    }
    size_t length = strlen(value);
    if (length >= capacity) {
        fprintf(stderr, "Value for %s is too long\n", label);
        return false;
    }
    memcpy(out, value, length + 1);
    return true;
}

static bool debug_parse_integer(const char *value, int minimum, int maximum,
                                int *out, const char *label)
{
    if (!value || !value[0] || !out) {
        fprintf(stderr, "Missing integer for %s\n", label);
        return false;
    }
    char *end = NULL;
    errno = 0;
    long parsed = strtol(value, &end, 10);
    if (errno != 0 || !end || *end != '\0' || parsed < minimum || parsed > maximum) {
        fprintf(stderr, "Invalid %s '%s' (expected %d..%d)\n",
                label, value, minimum, maximum);
        return false;
    }
    *out = (int)parsed;
    return true;
}

static int debug_compare_ints(const void *left, const void *right)
{
    int a = *(const int *)left;
    int b = *(const int *)right;
    return (a > b) - (a < b);
}

static bool debug_parse_capture_steps(const char *value)
{
    if (!value || !value[0]) {
        fprintf(stderr, "Missing value for --debug-capture-steps\n");
        return false;
    }
    char copy[1024];
    if (strlen(value) >= sizeof(copy)) {
        fprintf(stderr, "--debug-capture-steps is too long\n");
        return false;
    }
    strcpy(copy, value);
    debugRunConfig.capture_step_count = 0;
    char *token = strtok(copy, ",");
    while (token) {
        if (debugRunConfig.capture_step_count >= DEBUG_MAX_CAPTURE_STEPS) {
            fprintf(stderr, "Too many debug capture steps (maximum %d)\n", DEBUG_MAX_CAPTURE_STEPS);
            return false;
        }
        int step = 0;
        if (!debug_parse_integer(token, 0, 10000, &step, "debug capture step")) {
            return false;
        }
        debugRunConfig.capture_steps[debugRunConfig.capture_step_count++] = step;
        token = strtok(NULL, ",");
    }
    if (debugRunConfig.capture_step_count <= 0) {
        fprintf(stderr, "--debug-capture-steps requires at least one step\n");
        return false;
    }
    qsort(debugRunConfig.capture_steps, (size_t)debugRunConfig.capture_step_count,
          sizeof(debugRunConfig.capture_steps[0]), debug_compare_ints);
    int write = 0;
    for (int i = 0; i < debugRunConfig.capture_step_count; ++i) {
        if (write == 0 || debugRunConfig.capture_steps[i] != debugRunConfig.capture_steps[write - 1]) {
            debugRunConfig.capture_steps[write++] = debugRunConfig.capture_steps[i];
        }
    }
    debugRunConfig.capture_step_count = write;
    debugRunConfig.capture_steps_explicit = true;
    return true;
}

static const char *debug_option_value(int argc, char **argv, int *index,
                                      const char *argument, const char *prefix,
                                      const char *separate)
{
    size_t prefix_length = strlen(prefix);
    if (strncmp(argument, prefix, prefix_length) == 0) {
        return argument + prefix_length;
    }
    if (strcmp(argument, separate) == 0) {
        if (*index + 1 >= argc) {
            return NULL;
        }
        ++(*index);
        return argv[*index];
    }
    return (const char *)(uintptr_t)1;
}

static int debug_parse_argument(int argc, char **argv, int *index)
{
    const char *argument = argv[*index];
    if (strcmp(argument, "--debug-show-window") == 0) {
        debugRunConfig.show_window = true;
        return 1;
    }
    if (strcmp(argument, "--debug-matrix") == 0) {
        fprintf(stderr, "--debug-matrix must be run through the fps_ray launcher\n");
        return -1;
    }

    const char *value = debug_option_value(argc, argv, index, argument,
                                           "--debug-scenario=", "--debug-scenario");
    if (value != (const char *)(uintptr_t)1) {
        if (!debug_copy_option(debugRunConfig.scenario, sizeof(debugRunConfig.scenario),
                               value, "--debug-scenario")) return -1;
        debugRunConfig.requested = true;
        return 1;
    }

    value = debug_option_value(argc, argv, index, argument,
                               "--debug-steps=", "--debug-steps");
    if (value != (const char *)(uintptr_t)1) {
        if (!debug_parse_integer(value, DEBUG_MIN_STEPS, 10000,
                                 &debugRunConfig.steps, "--debug-steps")) return -1;
        return 1;
    }

    value = debug_option_value(argc, argv, index, argument,
                               "--debug-capture-steps=", "--debug-capture-steps");
    if (value != (const char *)(uintptr_t)1) {
        return debug_parse_capture_steps(value) ? 1 : -1;
    }

    value = debug_option_value(argc, argv, index, argument,
                               "--debug-output=", "--debug-output");
    if (value != (const char *)(uintptr_t)1) {
        return debug_copy_option(debugRunConfig.output_dir, sizeof(debugRunConfig.output_dir),
                                 value, "--debug-output") ? 1 : -1;
    }

    if (strncmp(argument, "--debug-", 8) == 0) {
        fprintf(stderr, "Unknown debug option '%s'\n", argument);
        return -1;
    }
    return 0;
}

static bool debug_run_requested(void)
{
    return debugRunConfig.requested;
}

static bool debug_show_window_requested(void)
{
    return debugRunConfig.show_window;
}

static const DebugScenarioDefinition *debug_find_scenario(const char *name)
{
    for (int i = 0; i < debugScenarioCount; ++i) {
        if (strcmp(debugScenarios[i].name, name) == 0) {
            return &debugScenarios[i];
        }
    }
    return NULL;
}

static bool debug_join_path(char *out, size_t capacity, const char *left, const char *right)
{
    if (!out || capacity == 0 || !left || !right) return false;
    size_t left_length = strlen(left);
    bool needs_separator = left_length > 0 && left[left_length - 1] != '/' && left[left_length - 1] != '\\';
    int written = snprintf(out, capacity, "%s%s%s", left, needs_separator ? "/" : "", right);
    return written >= 0 && (size_t)written < capacity;
}

static void debug_add_failure(DebugRunState *state, const char *format, ...)
{
    if (!state || !format || state->failure_count >= DEBUG_MAX_FAILURES) return;
    char message[DEBUG_MESSAGE_CAPACITY];
    va_list args;
    va_start(args, format);
    vsnprintf(message, sizeof(message), format, args);
    va_end(args);
    for (int i = 0; i < state->failure_count; ++i) {
        if (strcmp(state->failures[i], message) == 0) return;
    }
    snprintf(state->failures[state->failure_count++], DEBUG_MESSAGE_CAPACITY, "%s", message);
}

static int debug_add_cell(int gx, int gy, int gz, bool fixed, bool simulate,
                          Color color, int tag)
{
    float px = ((float)gx + 0.5f) * VOXEL_SIZE;
    float py = ((float)gy + 0.5f) * VOXEL_SIZE;
    float pz = ((float)gz + 0.5f) * VOXEL_SIZE;
    int index = addVoxel(px, py, pz, fixed, simulate, color, 0);
    if (index >= 0) voxels[index].debugClusterTag = tag;
    return index;
}

static int debug_add_test_cluster(bool simulate)
{
    int created = 0;
    Color color = simulate ? (Color){ 245, 145, 45, 255 } : (Color){ 80, 175, 245, 255 };
    for (int y = 12; y <= 13; ++y) {
        for (int z = -1; z <= 1; ++z) {
            for (int x = -1; x <= 1; ++x) {
                if (debug_add_cell(x, y, z, !simulate, simulate, color, DEBUG_TEST_TAG) >= 0) {
                    ++created;
                }
            }
        }
    }
    return created;
}

static int debug_add_support_column(void)
{
    int created = 0;
    Color support_color = (Color){ 120, 125, 135, 255 };
    for (int y = 0; y <= 11; ++y) {
        if (debug_add_cell(0, y, 0, true, false, support_color, DEBUG_SUPPORT_TAG) >= 0) {
            ++created;
        }
    }
    return created;
}

static int debug_count_tagged(int tag, int simulate_filter)
{
    int count = 0;
    for (int i = 0; i < voxel_count; ++i) {
        if (voxels[i].debugClusterTag != tag) continue;
        if (simulate_filter >= 0 && (voxels[i].simulate ? 1 : 0) != simulate_filter) continue;
        ++count;
    }
    return count;
}

static void debug_track_tagged_cells(DebugRunState *state, int tag)
{
    if (!state) return;
    state->tracked_cell_count = 0;
    for (int i = 0; i < voxel_count && state->tracked_cell_count < DEBUG_MAX_TRACKED_CELLS; ++i) {
        const Voxel *voxel = &voxels[i];
        if (voxel->debugClusterTag != tag) continue;
        state->tracked_cells[state->tracked_cell_count++] = (DebugTrackedCell){
            voxel->gx, voxel->gy, voxel->gz
        };
    }
}

static int debug_find_tagged_voxel(int tag, bool require_static)
{
    for (int i = 0; i < voxel_count; ++i) {
        if (voxels[i].debugClusterTag != tag) continue;
        if (require_static && voxels[i].simulate) continue;
        return i;
    }
    return -1;
}

static void debug_rebuild_world_state(void)
{
    rebuild_voxel_hash();
    rebuild_static_hash_if_dirty();
    rebuild_all_voxel_surfaces();
    rebuild_glue_constraints();
    refresh_static_voxel_beliefs();
    meshDirty = true;
}

static void debug_reset_world(void)
{
    clear_world_voxels();
    clear_pickups();
    glueConstraintCount = 0;
    voxelHashFramesSinceRebuild = 0;
    pbdTimeAccumulator = 0.0f;
    staticHashDirty = false;
    init_static_hash();
    table_cache_invalidate();
}

static bool debug_remove_tagged_voxels(int tag)
{
    bool removed = false;
    int i = 0;
    while (i < voxel_count) {
        if (voxels[i].debugClusterTag == tag) {
            remove_voxel_index(i);
            removed = true;
        } else {
            ++i;
        }
    }
    if (removed) debug_rebuild_world_state();
    return removed;
}

static bool debug_original_cluster_cells_are_clear(const DebugRunState *state)
{
    if (!state || state->tracked_cell_count <= 0) return false;
    for (int i = 0; i < state->tracked_cell_count; ++i) {
        const DebugTrackedCell *cell = &state->tracked_cells[i];
        if (table_get_static_only(cell->x, cell->y, cell->z) >= 0) return false;
    }
    return true;
}

#if defined(GRAPHICS_API_OPENGL_43)
static bool debug_gpu_static_cell_absent_in(const GpuStaticCell *cells, int x, int y, int z)
{
    if (!cells) return false;
    uint32_t hash = gpu_static_hash_coord(x, y, z);
    for (int probe = 0; probe < HASH_SIZE; ++probe) {
        const GpuStaticCell *entry = &cells[hash];
        if (entry->value[3] < 0) return true;
        if (entry->value[0] == x && entry->value[1] == y && entry->value[2] == z) return false;
        hash = (hash + 1u) & (HASH_SIZE - 1u);
    }
    return true;
}

static bool debug_gpu_cluster_cells_are_clear(const DebugRunState *state)
{
    if (!state || state->tracked_cell_count <= 0) return false;
    for (int i = 0; i < state->tracked_cell_count; ++i) {
        const DebugTrackedCell *cell = &state->tracked_cells[i];
        if (!debug_gpu_static_cell_absent_in(gpuPhysics.static_stage,
                                             cell->x, cell->y, cell->z)) return false;
    }
    return true;
}
#endif

static DebugSample debug_collect_sample(const char *phase, int step, float initial_centroid_y)
{
    DebugSample sample = { 0 };
    snprintf(sample.phase, sizeof(sample.phase), "%s", phase ? phase : "simulation");
    sample.step = step;
    sample.sim_particles = sim_particle_count;
    sample.min_y = FLT_MAX;
    sample.max_y = -FLT_MAX;
    sample.static_hash_generation = staticHashGeneration;
#if defined(GRAPHICS_API_OPENGL_43)
    sample.gpu_static_uploaded_generation = gpuPhysics.static_uploaded_generation;
#else
    sample.gpu_static_uploaded_generation = 0;
#endif

    float centroid_sum = 0.0f;
    float velocity_sum = 0.0f;
    for (int i = 0; i < voxel_count; ++i) {
        const Voxel *voxel = &voxels[i];
        if (voxel->debugClusterTag != DEBUG_TEST_TAG) continue;
        ++sample.tagged_voxels;
        if (voxel->simulate) ++sample.tagged_dynamic;
        else ++sample.tagged_static;

        Vector3 center = { 0.0f, 0.0f, 0.0f };
        int valid_corners = 0;
        for (int corner = 0; corner < VOXEL_CORNER_COUNT; ++corner) {
            Particle *particle = voxel->particles[corner];
            if (!particle || !particle->active || !v_isfinite(particle->pos) ||
                !v_isfinite(particle->vel)) {
                ++sample.invalid_particle_refs;
                continue;
            }
            center = v_add(center, particle->pos);
            ++valid_corners;
            if (voxel->simulate && particle->base_inv_mass <= 0.0f) {
                ++sample.zero_mass_corners;
            }
            if (voxel->simulate &&
                (particle->sim_index < 0 || particle->sim_index >= sim_particle_count ||
                 sim_particles[particle->sim_index] != particle)) {
                ++sample.invalid_particle_refs;
            }
        }
        if (valid_corners > 0) center = v_mul(center, 1.0f / (float)valid_corners);
        centroid_sum += center.y;
        velocity_sum += voxel->vel.y;
        if (center.y < sample.min_y) sample.min_y = center.y;
        if (center.y > sample.max_y) sample.max_y = center.y;

        if (voxel->simulate) {
            for (int face = 0; face < 6; ++face) {
                if (!voxel->glued_faces[face]) continue;
                int nx = voxel->gx + face_offsets[face][0];
                int ny = voxel->gy + face_offsets[face][1];
                int nz = voxel->gz + face_offsets[face][2];
                if (table_get_static_only(nx, ny, nz) >= 0) ++sample.static_glued_faces;
            }
        }
    }
    if (sample.tagged_voxels > 0) {
        sample.centroid_y = centroid_sum / (float)sample.tagged_voxels;
        sample.mean_velocity_y = velocity_sum / (float)sample.tagged_voxels;
    } else {
        sample.min_y = 0.0f;
        sample.max_y = 0.0f;
    }
    sample.centroid_drop = initial_centroid_y - sample.centroid_y;
    sample.physics_step_ms = physicsBackend.last_step_ms;
    return sample;
}

static DebugSample *debug_record_sample(DebugRunState *state, const char *phase, int step)
{
    if (!state || state->sample_count >= DEBUG_MAX_SAMPLES) return NULL;
    DebugSample *sample = &state->samples[state->sample_count++];
    *sample = debug_collect_sample(phase, step, state->initial_centroid_y);
    return sample;
}

static void debug_draw_overlay(const DebugRunState *state, const DebugSample *sample)
{
    if (!state || !sample) return;
    int box_height = 180;
    if (state->gpu_readback_performed) box_height += 30;
    if (state->overhang_impact) box_height += 30;
    DrawRectangle(18, 18, 650, box_height, Fade(BLACK, 0.78f));
    DrawText(TextFormat("SCENARIO %s", debugRunConfig.scenario), 34, 32, 22, RAYWHITE);
    DrawText(TextFormat("BACKEND requested=%s active=%s",
                        physics_backend_name(physicsBackend.requested),
                        physics_backend_name(physicsBackend.active)),
             34, 62, 19, LIGHTGRAY);
    DrawText(TextFormat("PHASE %s | STEP %d/%d", sample->phase, sample->step,
                        debugRunConfig.steps),
             34, 88, 19, LIGHTGRAY);
    DrawText(TextFormat("TAGGED %d dyn=%d static=%d | SIM PARTICLES %d",
                        sample->tagged_voxels, sample->tagged_dynamic,
                        sample->tagged_static, sample->sim_particles),
             34, 114, 19, LIGHTGRAY);
    DrawText(TextFormat("CENTROID Y %.3f | DROP %.3f | VEL Y %.3f",
                        sample->centroid_y, sample->centroid_drop,
                        sample->mean_velocity_y),
             34, 140, 19, LIGHTGRAY);
    Color result_color = state->failure_count == 0 ? LIME : RED;
    DrawText(TextFormat("%s | invalid=%d zeroMass=%d staticGlue=%d",
                        state->failure_count == 0 ? "PASS SO FAR" : "FAIL",
                        sample->invalid_particle_refs, sample->zero_mass_corners,
                        sample->static_glued_faces),
             34, 166, 19, result_color);
    if (state->gpu_readback_performed) {
        DrawText(TextFormat("GPU SSBO staticStale=%d staticDiff=%d simMismatch=%d missing=%d gpuZeroMass=%d",
                            state->gpu_static_ssbo_stale_cells,
                            state->gpu_static_ssbo_diff_entries,
                            state->gpu_sim_id_mismatches,
                            state->gpu_tagged_corners_missing_from_sim_ids,
                            state->gpu_zero_mass_corners),
                 34, 192, 17, result_color);
    }
    if (state->overhang_impact) {
        int y = state->gpu_readback_performed ? 222 : 192;
        DrawText(TextFormat("OVERHANG rootDrop30=%.3f tipDrop30=%.3f relative=%.3f anchor=%s",
                            state->overhang_root_drop_at_step_30,
                            state->overhang_tip_drop_at_step_30,
                            state->overhang_relative_deflection_at_step_30,
                            state->overhang_anchor_observed ? "yes" : "no"),
                 34, y, 17, result_color);
    }
}

static void debug_render_scene(DebugRunState *state, const DebugSample *sample)
{
    BeginTextureMode(state->target);
        ClearBackground((Color){ 32, 42, 58, 255 });
        BeginMode3D(state->camera);
            DrawPlane((Vector3){ 0.0f, 0.0f, 0.0f }, (Vector2){ 24.0f, 24.0f },
                      (Color){ 68, 75, 82, 255 });
            DrawGrid(48, VOXEL_SIZE);
            DrawVoxels(state->camera);
            DrawLine3D((Vector3){ -4.0f, state->initial_centroid_y, 0.0f },
                       (Vector3){ 4.0f, state->initial_centroid_y, 0.0f },
                       (Color){ 255, 225, 80, 210 });
        EndMode3D();
        debug_draw_overlay(state, sample);
    EndTextureMode();
}

static bool debug_capture(DebugRunState *state, const char *file_name,
                          const char *phase, int step)
{
    if (!state || !file_name || state->capture_count >= DEBUG_MAX_CAPTURES) return false;
    DebugSample *sample = debug_record_sample(state, phase, step);
    if (!sample) {
        debug_add_failure(state, "sample capacity exceeded while capturing %s", file_name);
        return false;
    }
    debug_render_scene(state, sample);

    char path[DEBUG_PATH_CAPACITY];
    if (!debug_join_path(path, sizeof(path), state->output_dir, file_name)) {
        debug_add_failure(state, "capture path is too long for %s", file_name);
        return false;
    }
    Image image = LoadImageFromTexture(state->target.texture);
    if (!IsImageValid(image)) {
        debug_add_failure(state, "unable to read render texture for %s", file_name);
        return false;
    }
    ImageFlipVertical(&image);
    bool exported = ExportImage(image, path);
    UnloadImage(image);
    if (!exported) {
        debug_add_failure(state, "unable to export screenshot %s", path);
        return false;
    }
    snprintf(state->captures[state->capture_count++], DEBUG_PATH_CAPACITY, "%s", path);
    return true;
}

static bool debug_should_capture_step(int step)
{
    for (int i = 0; i < debugRunConfig.capture_step_count; ++i) {
        if (debugRunConfig.capture_steps[i] == step) return true;
    }
    return false;
}

static bool debug_prepare_capture_steps(void)
{
    if (!debugRunConfig.capture_steps_explicit) {
        const int defaults[] = { 0, 1, 5, 15, 30, 60 };
        debugRunConfig.capture_step_count = 0;
        for (int i = 0; i < (int)(sizeof(defaults) / sizeof(defaults[0])); ++i) {
            if (defaults[i] <= debugRunConfig.steps) {
                debugRunConfig.capture_steps[debugRunConfig.capture_step_count++] = defaults[i];
            }
        }
        if (debugRunConfig.capture_steps[debugRunConfig.capture_step_count - 1] != debugRunConfig.steps) {
            debugRunConfig.capture_steps[debugRunConfig.capture_step_count++] = debugRunConfig.steps;
        }
    }
    for (int i = 0; i < debugRunConfig.capture_step_count; ++i) {
        if (debugRunConfig.capture_steps[i] > debugRunConfig.steps) {
            fprintf(stderr, "Debug capture step %d exceeds --debug-steps=%d\n",
                    debugRunConfig.capture_steps[i], debugRunConfig.steps);
            return false;
        }
    }
    return true;
}

static void debug_validate_tagged_dynamic(DebugRunState *state, const DebugSample *sample,
                                          const char *context)
{
    if (!state || !sample) return;
    if (sample->tagged_voxels != state->expected_voxels) {
        debug_add_failure(state, "%s: expected %d tagged voxels, found %d",
                          context, state->expected_voxels, sample->tagged_voxels);
    }
    if (sample->tagged_dynamic != state->expected_voxels || sample->tagged_static != 0) {
        debug_add_failure(state, "%s: tagged cluster is not fully dynamic (dynamic=%d static=%d)",
                          context, sample->tagged_dynamic, sample->tagged_static);
    }
    if (sample->invalid_particle_refs > 0) {
        debug_add_failure(state, "%s: %d invalid particle or sim-index references",
                          context, sample->invalid_particle_refs);
    }
    if (sample->zero_mass_corners > 0) {
        debug_add_failure(state, "%s: %d activated corners have zero base inverse mass",
                          context, sample->zero_mass_corners);
    }
    if (sample->static_glued_faces > 0) {
        debug_add_failure(state, "%s: %d activated faces remain glued to static voxels",
                          context, sample->static_glued_faces);
    }
}

static void debug_validate_overhang_activation(DebugRunState *state,
                                               const DebugSample *sample)
{
    if (!state || !sample) return;
    if (sample->tagged_dynamic != state->expected_voxels || sample->tagged_static != 0) {
        debug_add_failure(state,
                          "overhang activation is incomplete (dynamic=%d static=%d expected=%d)",
                          sample->tagged_dynamic, sample->tagged_static, state->expected_voxels);
    }
    int unexpected_invalid = sample->invalid_particle_refs - sample->zero_mass_corners;
    if (unexpected_invalid > 0) {
        debug_add_failure(state, "overhang activation has %d non-anchor particle reference errors",
                          unexpected_invalid);
    }
    if (sample->zero_mass_corners <= 0 || sample->static_glued_faces <= 0) {
        debug_add_failure(state,
                          "overhang activation did not preserve its static root anchor (zeroMass=%d staticGlue=%d)",
                          sample->zero_mass_corners, sample->static_glued_faces);
    } else {
        state->overhang_anchor_observed = true;
    }
}

static void debug_validate_all_sim_particles(DebugRunState *state, int step)
{
    for (int i = 0; i < sim_particle_count; ++i) {
        Particle *particle = sim_particles[i];
        if (!particle || !particle->active || particle->sim_index != i ||
            !v_isfinite(particle->pos) || !v_isfinite(particle->vel) ||
            !v_isfinite(particle->predicted_pos)) {
            debug_add_failure(state, "step %d: invalid global sim particle at index %d", step, i);
            return;
        }
    }
}

static void debug_validate_gpu_static_state(DebugRunState *state, int step)
{
    if (physicsBackend.active != PHYSICS_BACKEND_GPU_GL43 || step < 1 ||
        !state || !state->activation_succeeded || state->gpu_readback_performed) return;
#if defined(GRAPHICS_API_OPENGL_43)
    if (gpuPhysics.static_uploaded_generation != staticHashGeneration) {
        debug_add_failure(state,
                          "step %d: GPU static generation %llu does not match CPU generation %llu",
                          step,
                          (unsigned long long)gpuPhysics.static_uploaded_generation,
                          (unsigned long long)staticHashGeneration);
    }
    if (!debug_gpu_cluster_cells_are_clear(state)) {
        debug_add_failure(state, "step %d: activated cells remain in the staged GPU static hash", step);
    }
#else
    debug_add_failure(state, "step %d: GPU backend active in a non-GPU build", step);
#endif
}

static void debug_readback_gpu_buffers(DebugRunState *state, int step)
{
    if (!state || physicsBackend.active != PHYSICS_BACKEND_GPU_GL43 || step < 1 ||
        !state->activation_succeeded || state->gpu_readback_performed) return;
#if defined(GRAPHICS_API_OPENGL_43)
    state->gpu_readback_performed = true;
    if (gpuPhysics.static_cell_buffer.id == 0 || gpuPhysics.sim_id_buffer.id == 0 ||
        gpuPhysics.particle_buffer.id == 0 ||
        !gpuPhysics.static_stage || !gpuPhysics.sim_ids || !gpuPhysics.particle_stage ||
        !gpuPhysics.pool_to_compact) {
        debug_add_failure(state, "step %d: GPU diagnostic buffers are not initialized", step);
        return;
    }

    size_t static_bytes = (size_t)HASH_SIZE * sizeof(GpuStaticCell);
    GpuStaticCell *static_readback = malloc(static_bytes);
    size_t sim_bytes = (size_t)sim_particle_count * sizeof(uint32_t);
    uint32_t *sim_id_readback = sim_bytes ? malloc(sim_bytes) : NULL;
    size_t particle_bytes = (size_t)gpuPhysics.compact_count * sizeof(GpuParticleState);
    GpuParticleState *particle_readback = particle_bytes ? malloc(particle_bytes) : NULL;
    if (!static_readback || (sim_bytes && !sim_id_readback) ||
        (particle_bytes && !particle_readback)) {
        free(static_readback);
        free(sim_id_readback);
        free(particle_readback);
        debug_add_failure(state, "step %d: unable to allocate GPU diagnostic readback storage", step);
        return;
    }

    rlReadShaderBuffer(gpuPhysics.static_cell_buffer.id, static_readback,
                       (unsigned int)static_bytes, 0);
    if (sim_bytes) {
        rlReadShaderBuffer(gpuPhysics.sim_id_buffer.id, sim_id_readback,
                           (unsigned int)sim_bytes, 0);
    }
    if (particle_bytes) {
        rlReadShaderBuffer(gpuPhysics.particle_buffer.id, particle_readback,
                           (unsigned int)particle_bytes, 0);
    }

    for (int i = 0; i < HASH_SIZE; ++i) {
        if (memcmp(&static_readback[i], &gpuPhysics.static_stage[i],
                   sizeof(static_readback[i])) != 0) {
            ++state->gpu_static_ssbo_diff_entries;
        }
    }
    for (int i = 0; i < state->tracked_cell_count; ++i) {
        const DebugTrackedCell *cell = &state->tracked_cells[i];
        if (!debug_gpu_static_cell_absent_in(static_readback,
                                             cell->x, cell->y, cell->z)) {
            ++state->gpu_static_ssbo_stale_cells;
        }
    }

    for (int i = 0; i < sim_particle_count; ++i) {
        uint32_t actual = sim_id_readback[i];
        if (actual != gpuPhysics.sim_ids[i] || actual >= (uint32_t)gpuPhysics.compact_count) {
            ++state->gpu_sim_id_mismatches;
        }
    }
    for (int i = 0; i < voxel_count; ++i) {
        const Voxel *voxel = &voxels[i];
        if (!voxel->simulate || voxel->debugClusterTag != DEBUG_TEST_TAG) continue;
        for (int corner = 0; corner < VOXEL_CORNER_COUNT; ++corner) {
            const Particle *particle = voxel->particles[corner];
            if (!particle) {
                ++state->gpu_tagged_corners_missing_from_sim_ids;
                continue;
            }
            ptrdiff_t pool_id = particle - particles_pool;
            if (pool_id < 0 || pool_id >= particle_pool_count) {
                ++state->gpu_tagged_corners_missing_from_sim_ids;
                continue;
            }
            int compact_id = gpuPhysics.pool_to_compact[pool_id];
            if (compact_id < 0 || compact_id >= gpuPhysics.compact_count) {
                ++state->gpu_tagged_corners_missing_from_sim_ids;
                continue;
            }
            bool zero_mass =
                particle_readback[compact_id].predicted_base_inv_mass[3] <= 0.0f;
            if (zero_mass) {
                ++state->gpu_zero_mass_corners;
                if (state->overhang_impact) continue;
            }
            bool found = false;
            for (int sim = 0; sim < sim_particle_count; ++sim) {
                if (sim_id_readback[sim] == (uint32_t)compact_id) {
                    found = true;
                    break;
                }
            }
            if (!found) ++state->gpu_tagged_corners_missing_from_sim_ids;
        }
    }

    free(static_readback);
    free(sim_id_readback);
    free(particle_readback);

    if (state->gpu_static_ssbo_diff_entries > 0) {
        debug_add_failure(state, "step %d: %d GPU static SSBO entries differ from the CPU staging table",
                          step, state->gpu_static_ssbo_diff_entries);
    }
    if (state->gpu_static_ssbo_stale_cells > 0) {
        debug_add_failure(state, "step %d: %d activated cluster cells remain in the GPU static SSBO",
                          step, state->gpu_static_ssbo_stale_cells);
    }
    if (state->gpu_sim_id_mismatches > 0) {
        debug_add_failure(state, "step %d: %d uploaded simId values are invalid or differ from staging",
                          step, state->gpu_sim_id_mismatches);
    }
    if (state->gpu_tagged_corners_missing_from_sim_ids > 0) {
        debug_add_failure(state, "step %d: %d tagged dynamic corners are missing from the uploaded simId list",
                          step, state->gpu_tagged_corners_missing_from_sim_ids);
    }
    if (!state->overhang_impact && state->gpu_zero_mass_corners > 0) {
        debug_add_failure(state, "step %d: %d tagged dynamic corners have zero mass in the GPU particle SSBO",
                          step, state->gpu_zero_mass_corners);
    } else if (state->overhang_impact && state->gpu_zero_mass_corners <= 0) {
        debug_add_failure(state, "step %d: GPU particle SSBO does not retain the overhang root anchors",
                          step);
    }
#else
    (void)state;
    (void)step;
#endif
}

static bool debug_activate_test_cluster(DebugRunState *state)
{
    int seed = debug_find_tagged_voxel(DEBUG_TEST_TAG, true);
    if (seed < 0) {
        debug_add_failure(state, "unable to find a static test-cluster seed for activation");
        return false;
    }
    bool activated = activate_static_voxel_for_tether(seed, 0, ACTIVATION_TETHER_BELIEF);
    state->activation_succeeded = activated;
    if (!activated) {
        debug_add_failure(state, "tether activation rejected the unsupported test cluster");
        return false;
    }
    rebuild_static_hash_if_dirty();
    if (!debug_original_cluster_cells_are_clear(state)) {
        debug_add_failure(state, "activated cells remain in the CPU static hash");
    }
    return true;
}

static void debug_set_initial_centroid(DebugRunState *state)
{
    DebugSample sample = debug_collect_sample("initial", 0, 0.0f);
    state->initial_centroid_y = sample.centroid_y;
}

static bool debug_tagged_original_column_y(int tag, int original_x, float *out_y)
{
    if (!out_y) return false;
    float sum = 0.0f;
    int count = 0;
    for (int i = 0; i < voxel_count; ++i) {
        const Voxel *voxel = &voxels[i];
        if (voxel->debugClusterTag != tag || voxel->orig_min_gx != original_x) continue;
        float center_y = 0.0f;
        int corners = 0;
        for (int c = 0; c < VOXEL_CORNER_COUNT; ++c) {
            const Particle *particle = voxel->particles[c];
            if (!particle || !particle->active || !v_isfinite(particle->pos)) continue;
            center_y += particle->pos.y;
            ++corners;
        }
        if (corners <= 0) continue;
        sum += center_y / (float)corners;
        ++count;
    }
    if (count <= 0) return false;
    *out_y = sum / (float)count;
    return true;
}

static bool debug_setup_dynamic_freefall(DebugRunState *state)
{
    debug_reset_world();
    int created = debug_add_test_cluster(true);
    debug_track_tagged_cells(state, DEBUG_TEST_TAG);
    debug_rebuild_world_state();
    state->activation_succeeded = true;
    debug_set_initial_centroid(state);
    if (created != state->expected_voxels) {
        debug_add_failure(state, "dynamic-freefall created %d/%d voxels",
                          created, state->expected_voxels);
        return false;
    }
    return true;
}

static bool debug_setup_activation_floating(DebugRunState *state)
{
    debug_reset_world();
    int created = debug_add_test_cluster(false);
    debug_track_tagged_cells(state, DEBUG_TEST_TAG);
    debug_rebuild_world_state();
    debug_set_initial_centroid(state);
    if (created != state->expected_voxels) {
        debug_add_failure(state, "activation-floating created %d/%d voxels",
                          created, state->expected_voxels);
        return false;
    }
    debug_capture(state, "phase-01-floating-static.png", "floating-static", 0);
    if (!debug_activate_test_cluster(state)) return false;
    debug_set_initial_centroid(state);
    return true;
}

static bool debug_setup_sleep_wake_floating(DebugRunState *state)
{
    debug_reset_world();
    int support_created = debug_add_support_column();
    debug_rebuild_world_state();
    int cluster_created = debug_add_test_cluster(true);
    debug_track_tagged_cells(state, DEBUG_TEST_TAG);
    debug_rebuild_world_state();
    debug_set_initial_centroid(state);

    if (support_created != 12 || cluster_created != state->expected_voxels) {
        debug_add_failure(state, "sleep-wake setup created support=%d/12 cluster=%d/%d",
                          support_created, cluster_created, state->expected_voxels);
        return false;
    }
    debug_capture(state, "phase-01-supported-dynamic.png", "supported-dynamic", 0);

    for (int i = 0; i < voxel_count; ++i) {
        if (voxels[i].debugClusterTag == DEBUG_TEST_TAG && voxels[i].simulate) {
            voxels[i].sleepFrames = VOXEL_DEACTIVATION_FRAMES - 1;
            set_voxel_velocity(&voxels[i], (Vector3){ 0.0f, 0.0f, 0.0f });
        }
    }
    state->sleeping_conversion_succeeded = deactivate_sleeping_voxels();
    rebuild_static_hash_if_dirty();
    if (!state->sleeping_conversion_succeeded ||
        debug_count_tagged(DEBUG_TEST_TAG, 0) != state->expected_voxels) {
        debug_add_failure(state, "sleep conversion did not produce %d tagged static voxels",
                          state->expected_voxels);
        return false;
    }
    debug_capture(state, "phase-02-sleeping-static.png", "sleeping-static", 0);

    state->support_removal_succeeded = debug_remove_tagged_voxels(DEBUG_SUPPORT_TAG);
    if (!state->support_removal_succeeded || debug_count_tagged(DEBUG_SUPPORT_TAG, -1) != 0) {
        debug_add_failure(state, "temporary support column was not completely removed");
        return false;
    }
    debug_capture(state, "phase-03-unsupported-static.png", "unsupported-static", 0);

    if (!debug_activate_test_cluster(state)) return false;
    debug_set_initial_centroid(state);
    return true;
}

static bool debug_setup_tether_throw_floating(DebugRunState *state)
{
    debug_reset_world();
    int support_created = debug_add_support_column();
    debug_rebuild_world_state();
    int cluster_created = debug_add_test_cluster(true);
    debug_track_tagged_cells(state, DEBUG_TEST_TAG);
    debug_rebuild_world_state();
    debug_set_initial_centroid(state);

    if (support_created != 12 || cluster_created != state->expected_voxels) {
        debug_add_failure(state, "tether-throw setup created support=%d/12 cluster=%d/%d",
                          support_created, cluster_created, state->expected_voxels);
        return false;
    }
    debug_capture(state, "phase-01-supported-dynamic.png", "supported-dynamic", 0);

    for (int i = 0; i < voxel_count; ++i) {
        if (voxels[i].debugClusterTag == DEBUG_TEST_TAG && voxels[i].simulate) {
            voxels[i].sleepFrames = VOXEL_DEACTIVATION_FRAMES - 1;
            set_voxel_velocity(&voxels[i], (Vector3){ 0.0f, 0.0f, 0.0f });
        }
    }
    state->sleeping_conversion_succeeded = deactivate_sleeping_voxels();
    rebuild_static_hash_if_dirty();
    if (!state->sleeping_conversion_succeeded ||
        debug_count_tagged(DEBUG_TEST_TAG, 0) != state->expected_voxels) {
        debug_add_failure(state, "sleep conversion did not produce %d tagged static voxels",
                          state->expected_voxels);
        return false;
    }
    debug_capture(state, "phase-02-sleeping-static.png", "sleeping-static", 0);

    state->support_removal_succeeded = debug_remove_tagged_voxels(DEBUG_SUPPORT_TAG);
    if (!state->support_removal_succeeded || debug_count_tagged(DEBUG_SUPPORT_TAG, -1) != 0) {
        debug_add_failure(state, "temporary support column was not completely removed");
        return false;
    }
    debug_capture(state, "phase-03-unsupported-static.png", "unsupported-static", 0);

    int thrown_idx = debug_add_cell(-5, 12, 0, false, true,
                                    (Color){ 255, 210, 70, 255 }, DEBUG_THROW_TAG);
    if (thrown_idx < 0) {
        debug_add_failure(state, "unable to create the tether-thrown test voxel");
        return false;
    }
    debug_rebuild_world_state();
    Voxel *thrown = &voxels[thrown_idx];
    detach_all_glue_faces(thrown_idx);
    deactivate_glue_constraints_for_voxel(thrown_idx);
    thrown->owner = 0;
    thrown->activator = 0;
    thrown->wasTethered = true;
    thrown->activationBelief = 1.0f;

    memset(tetherTag, 0, sizeof(tetherTag));
    for (int i = 0; i < MAX_PLAYERS; ++i) {
        players[i].tetherHolding = false;
        players[i].tetherVoxel = -1;
    }
    players[0].yaw = -90.0f;
    players[0].pitch = 0.0f;
    players[0].tetherHolding = true;
    players[0].tetherVoxel = thrown_idx;
    prepare_tether_forces();
    release_tether(0);

    state->delayed_activation = true;
    state->activation_step = -1;
    debug_set_initial_centroid(state);
    debug_capture(state, "phase-04-tether-released.png", "tether-released", 0);
    return true;
}

static bool debug_setup_overhang_impact(DebugRunState *state)
{
    debug_reset_world();
    int support_created = 0;
    int overhang_created = 0;
    Color support_color = (Color){ 115, 120, 135, 255 };
    Color overhang_color = (Color){ 80, 175, 245, 255 };

    for (int y = 0; y <= 8; ++y) {
        for (int z = -1; z <= 1; ++z) {
            for (int x = -2; x <= 0; ++x) {
                if (debug_add_cell(x, y, z, true, false,
                                   support_color, DEBUG_SUPPORT_TAG) >= 0) {
                    ++support_created;
                }
            }
        }
    }
    for (int z = -1; z <= 1; ++z) {
        for (int x = 1; x <= 6; ++x) {
            if (debug_add_cell(x, 8, z, true, false,
                               overhang_color, DEBUG_TEST_TAG) >= 0) {
                ++overhang_created;
            }
        }
    }
    debug_track_tagged_cells(state, DEBUG_TEST_TAG);
    debug_rebuild_world_state();
    debug_set_initial_centroid(state);
    if (support_created != 81 || overhang_created != state->expected_voxels) {
        debug_add_failure(state, "overhang setup created support=%d/81 overhang=%d/%d",
                          support_created, overhang_created, state->expected_voxels);
        return false;
    }
    if (!debug_tagged_original_column_y(DEBUG_TEST_TAG, 1, &state->overhang_initial_root_y) ||
        !debug_tagged_original_column_y(DEBUG_TEST_TAG, 6, &state->overhang_initial_tip_y)) {
        debug_add_failure(state, "unable to measure the initial overhang root and tip");
        return false;
    }
    debug_capture(state, "phase-01-static-overhang.png", "static-overhang", 0);

    int impact_idx = debug_add_cell(5, 15, 0, false, true,
                                    (Color){ 255, 205, 65, 255 }, DEBUG_THROW_TAG);
    if (impact_idx < 0) {
        debug_add_failure(state, "unable to create the falling impact voxel");
        return false;
    }
    Voxel *impact = &voxels[impact_idx];
    impact->owner = 0;
    impact->activator = 0;
    impact->wasTethered = false;
    impact->isBullet = false;
    impact->activationBelief = 1.0f;
    set_voxel_velocity(impact, (Vector3){ 0.0f, -8.0f, 0.0f });
    debug_rebuild_world_state();

    state->delayed_activation = true;
    state->overhang_impact = true;
    state->activation_step = -1;
    debug_set_initial_centroid(state);
    debug_capture(state, "phase-02-falling-voxel.png", "falling-voxel", 0);
    return true;
}

static void debug_physics_step(void)
{
    activate_static_voxels_near_dynamic();
    prepare_tether_forces();
    rebuild_static_hash_if_dirty();
    simulate_voxel_pbd(PBD_MAX_STEP_DT);
    recycle_dead_voxels();
    rebuild_voxel_hash();
    deactivate_sleeping_voxels();
    rebuild_static_hash_if_dirty();
}

static void debug_validate_step(DebugRunState *state, int step)
{
    debug_validate_all_sim_particles(state, step);
    DebugSample sample = debug_collect_sample("validation", step, state->initial_centroid_y);
    if (state->delayed_activation) {
        if (!state->activation_succeeded &&
            sample.tagged_dynamic == state->expected_voxels && sample.tagged_static == 0) {
            state->activation_succeeded = true;
            state->activation_step = step;
        }
        if (state->activation_succeeded && step <= state->activation_step + 1) {
            if (state->overhang_impact) {
                debug_validate_overhang_activation(state, &sample);
            } else {
                debug_validate_tagged_dynamic(state, &sample, "gameplay-path activation");
            }
        }
    } else if (step <= 1) {
        debug_validate_tagged_dynamic(state, &sample, step == 0 ? "after activation" : "first physics step");
    }
    debug_validate_gpu_static_state(state, step);
    debug_readback_gpu_buffers(state, step);
    if (state->overhang_impact && step == 30) {
        float root_y = 0.0f;
        float tip_y = 0.0f;
        if (debug_tagged_original_column_y(DEBUG_TEST_TAG, 1, &root_y) &&
            debug_tagged_original_column_y(DEBUG_TEST_TAG, 6, &tip_y)) {
            state->overhang_root_drop_at_step_30 = state->overhang_initial_root_y - root_y;
            state->overhang_tip_drop_at_step_30 = state->overhang_initial_tip_y - tip_y;
            state->overhang_relative_deflection_at_step_30 =
                state->overhang_tip_drop_at_step_30 - state->overhang_root_drop_at_step_30;
            state->overhang_deflection_sampled = true;
        } else {
            debug_add_failure(state, "step 30: unable to measure overhang deflection");
        }
    }
    if (step == 60) {
        state->drop_sampled = true;
        state->drop_at_step_60 = sample.centroid_drop;
    }
}

static void debug_finalize_validation(DebugRunState *state)
{
    if (physicsBackend.requested == PHYSICS_BACKEND_GPU_GL43 &&
        physicsBackend.active != PHYSICS_BACKEND_GPU_GL43) {
        debug_add_failure(state, "forced GPU run fell back to %s: %s",
                          physics_backend_name(physicsBackend.active),
                          physicsBackend.fallback_reason);
    }
    if (!state->activation_succeeded) {
        debug_add_failure(state, "scenario did not complete its activation stage");
    }
    if (!state->drop_sampled) {
        debug_add_failure(state, "scenario did not reach the required 60-step drop sample");
    } else if (!state->overhang_impact &&
               (!isfinite(state->drop_at_step_60) || state->drop_at_step_60 < 0.5f)) {
        debug_add_failure(state, "cluster dropped only %.3f metres after 60 steps (minimum 0.500)",
                          state->drop_at_step_60);
    }
    if (state->overhang_impact) {
        if (!state->overhang_anchor_observed) {
            debug_add_failure(state, "overhang never retained a static root anchor after activation");
        }
        if (!state->overhang_deflection_sampled ||
            !isfinite(state->overhang_relative_deflection_at_step_30) ||
            state->overhang_relative_deflection_at_step_30 < 0.05f) {
            debug_add_failure(state,
                              "overhang relative tip deflection at step 30 was only %.3f metres (minimum 0.050)",
                              state->overhang_relative_deflection_at_step_30);
        }
    }
    if (strcmp(debugRunConfig.scenario, "sleep-wake-floating") == 0 ||
        (state->delayed_activation && !state->overhang_impact)) {
        if (!state->sleeping_conversion_succeeded) {
            debug_add_failure(state, "sleep-wake scenario never converted the cluster to static");
        }
        if (!state->support_removal_succeeded) {
            debug_add_failure(state, "sleep-wake scenario never removed its support");
        }
    }
}

static void debug_write_json_string(FILE *file, const char *value)
{
    fputc('"', file);
    const unsigned char *cursor = (const unsigned char *)(value ? value : "");
    while (*cursor) {
        unsigned char ch = *cursor++;
        switch (ch) {
            case '"': fputs("\\\"", file); break;
            case '\\': fputs("\\\\", file); break;
            case '\b': fputs("\\b", file); break;
            case '\f': fputs("\\f", file); break;
            case '\n': fputs("\\n", file); break;
            case '\r': fputs("\\r", file); break;
            case '\t': fputs("\\t", file); break;
            default:
                if (ch < 0x20) fprintf(file, "\\u%04x", (unsigned int)ch);
                else fputc((int)ch, file);
                break;
        }
    }
    fputc('"', file);
}

static bool debug_write_report(DebugRunState *state)
{
    FILE *file = fopen(state->report_path, "wb");
    if (!file) {
        debug_add_failure(state, "unable to create report %s: %s",
                          state->report_path, strerror(errno));
        return false;
    }
    bool passed = state->failure_count == 0;
    fputs("{\n  \"schemaVersion\": 1,\n  \"scenario\": ", file);
    debug_write_json_string(file, debugRunConfig.scenario);
    fputs(",\n  \"requestedBackend\": ", file);
    debug_write_json_string(file, physics_backend_name(physicsBackend.requested));
    fputs(",\n  \"activeBackend\": ", file);
    debug_write_json_string(file, physics_backend_name(physicsBackend.active));
    fputs(",\n  \"collisionProfile\": ", file);
    debug_write_json_string(file, collision_profile_name(collisionProfileActive));
    fprintf(file, ",\n  \"passed\": %s,\n  \"steps\": %d,\n  \"fixedDt\": %.9g,\n",
            passed ? "true" : "false", debugRunConfig.steps, (double)PBD_MAX_STEP_DT);
    fprintf(file, "  \"fallback\": %s,\n  \"fallbackReason\": ",
            physicsBackend.sticky_fallback ? "true" : "false");
    debug_write_json_string(file, physicsBackend.fallback_reason);
    fprintf(file,
            ",\n  \"activationSucceeded\": %s,\n  \"sleepingConversionSucceeded\": %s,\n"
            "  \"supportRemovalSucceeded\": %s,\n  \"activationStep\": %d,"
            "\n  \"dropAtStep60\": %.9g,\n",
            state->activation_succeeded ? "true" : "false",
            state->sleeping_conversion_succeeded ? "true" : "false",
            state->support_removal_succeeded ? "true" : "false",
            state->activation_step,
            (double)state->drop_at_step_60);

    fprintf(file,
            "  \"gpuDiagnostics\": {\"readbackPerformed\": %s,"
            " \"staticSsboDiffEntries\": %d, \"staticSsboStaleClusterCells\": %d,"
            " \"simIdMismatches\": %d, \"taggedCornersMissingFromSimIds\": %d,"
            " \"zeroMassTaggedCorners\": %d},\n",
            state->gpu_readback_performed ? "true" : "false",
            state->gpu_static_ssbo_diff_entries,
            state->gpu_static_ssbo_stale_cells,
            state->gpu_sim_id_mismatches,
            state->gpu_tagged_corners_missing_from_sim_ids,
            state->gpu_zero_mass_corners);

    fprintf(file,
            "  \"overhangDiagnostics\": {\"enabled\": %s, \"anchorObserved\": %s,"
            " \"deflectionSampled\": %s, \"rootDropAtStep30\": %.9g,"
            " \"tipDropAtStep30\": %.9g, \"relativeTipDeflectionAtStep30\": %.9g},\n",
            state->overhang_impact ? "true" : "false",
            state->overhang_anchor_observed ? "true" : "false",
            state->overhang_deflection_sampled ? "true" : "false",
            (double)state->overhang_root_drop_at_step_30,
            (double)state->overhang_tip_drop_at_step_30,
            (double)state->overhang_relative_deflection_at_step_30);

    fputs("  \"failures\": [", file);
    for (int i = 0; i < state->failure_count; ++i) {
        if (i > 0) fputs(", ", file);
        debug_write_json_string(file, state->failures[i]);
    }
    fputs("],\n  \"captures\": [", file);
    for (int i = 0; i < state->capture_count; ++i) {
        if (i > 0) fputs(", ", file);
        debug_write_json_string(file, state->captures[i]);
    }
    fputs("],\n  \"samples\": [\n", file);
    for (int i = 0; i < state->sample_count; ++i) {
        const DebugSample *sample = &state->samples[i];
        fputs("    {\"phase\": ", file);
        debug_write_json_string(file, sample->phase);
        fprintf(file,
                ", \"step\": %d, \"taggedVoxels\": %d, \"taggedDynamic\": %d,"
                " \"taggedStatic\": %d, \"simParticles\": %d,"
                " \"zeroMassCorners\": %d, \"invalidParticleRefs\": %d,"
                " \"staticGluedFaces\": %d, \"centroidY\": %.9g,"
                " \"centroidDrop\": %.9g, \"minY\": %.9g, \"maxY\": %.9g,"
                " \"meanVelocityY\": %.9g, \"staticHashGeneration\": %llu,"
                " \"gpuStaticUploadedGeneration\": %llu, \"physicsStepMs\": %.9g}%s\n",
                sample->step, sample->tagged_voxels, sample->tagged_dynamic,
                sample->tagged_static, sample->sim_particles,
                sample->zero_mass_corners, sample->invalid_particle_refs,
                sample->static_glued_faces, (double)sample->centroid_y,
                (double)sample->centroid_drop, (double)sample->min_y,
                (double)sample->max_y, (double)sample->mean_velocity_y,
                (unsigned long long)sample->static_hash_generation,
                (unsigned long long)sample->gpu_static_uploaded_generation,
                sample->physics_step_ms,
                (i + 1 < state->sample_count) ? "," : "");
    }
    fputs("  ]\n}\n", file);
    int write_error = ferror(file);
    int close_error = fclose(file);
    bool ok = write_error == 0 && close_error == 0;
    if (!ok) {
        debug_add_failure(state, "failed while writing report %s", state->report_path);
    }
    return ok;
}

static bool debug_prepare_output_directory(DebugRunState *state)
{
    if (debugRunConfig.output_dir[0]) {
        snprintf(state->output_dir, sizeof(state->output_dir), "%s", debugRunConfig.output_dir);
    } else {
        char artifacts[DEBUG_PATH_CAPACITY];
        char scenario_dir[DEBUG_PATH_CAPACITY];
        if (!debug_join_path(artifacts, sizeof(artifacts), GetApplicationDirectory(), "debug-artifacts") ||
            !debug_join_path(scenario_dir, sizeof(scenario_dir), artifacts, debugRunConfig.scenario) ||
            !debug_join_path(state->output_dir, sizeof(state->output_dir), scenario_dir,
                             physics_backend_name(physicsBackend.active))) {
            return false;
        }
    }
    if (!DirectoryExists(state->output_dir) && MakeDirectory(state->output_dir) != 0) {
        fprintf(stderr, "Unable to create debug output directory '%s'\n", state->output_dir);
        return false;
    }
    return debug_join_path(state->report_path, sizeof(state->report_path),
                           state->output_dir, "report.json");
}

static int run_debug_harness(void)
{
    const DebugScenarioDefinition *scenario = debug_find_scenario(debugRunConfig.scenario);
    if (!scenario) {
        fprintf(stderr, "Unknown debug scenario '%s'. Available scenarios:\n", debugRunConfig.scenario);
        for (int i = 0; i < debugScenarioCount; ++i) {
            fprintf(stderr, "  %s - %s\n", debugScenarios[i].name, debugScenarios[i].description);
        }
        return 2;
    }
    if (!debug_prepare_capture_steps()) return 2;

    DebugRunState state = { 0 };
    state.expected_voxels = DEBUG_EXPECTED_VOXELS;
    state.camera.position = (Vector3){ 9.0f, 8.5f, 12.0f };
    state.camera.target = (Vector3){ 0.0f, 4.0f, 0.0f };
    state.camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    state.camera.fovy = 45.0f;
    state.camera.projection = CAMERA_PERSPECTIVE;

    if (!debug_prepare_output_directory(&state)) return 4;
    state.target = LoadRenderTexture(DEBUG_RENDER_WIDTH, DEBUG_RENDER_HEIGHT);
    if (!IsRenderTextureValid(state.target)) {
        fprintf(stderr, "Unable to create debug render target\n");
        return 4;
    }

    bool setup_ok = scenario->setup(&state);
    if (!setup_ok && state.failure_count == 0) {
        debug_add_failure(&state, "scenario setup failed without a diagnostic");
    }

    if (setup_ok) {
        debug_validate_step(&state, 0);
        if (debug_should_capture_step(0)) {
            if (state.delayed_activation) {
                if (state.overhang_impact) {
                    debug_capture(&state, "step-0000-falling-voxel.png", "falling-voxel", 0);
                } else {
                    debug_capture(&state, "step-0000-tether-released.png", "tether-released", 0);
                }
            } else {
                debug_capture(&state, "step-0000-activated.png", "activated", 0);
            }
        }
        for (int step = 1; step <= debugRunConfig.steps; ++step) {
            debug_physics_step();
            debug_validate_step(&state, step);
            if (debug_should_capture_step(step)) {
                char file_name[64];
                snprintf(file_name, sizeof(file_name), "step-%04d.png", step);
                debug_capture(&state, file_name, "simulation", step);
            }
        }
        debug_finalize_validation(&state);
    }

    bool report_ok = debug_write_report(&state);
    bool passed = setup_ok && report_ok && state.failure_count == 0;
    fprintf(stderr,
            "debug-run scenario=%s requested=%s active=%s result=%s steps=%d drop60=%.6f captures=%d report=%s\n",
            debugRunConfig.scenario,
            physics_backend_name(physicsBackend.requested),
            physics_backend_name(physicsBackend.active),
            passed ? "PASS" : "FAIL",
            debugRunConfig.steps, state.drop_at_step_60,
            state.capture_count, state.report_path);
    for (int i = 0; i < state.failure_count; ++i) {
        fprintf(stderr, "debug-failure: %s\n", state.failures[i]);
    }
    UnloadRenderTexture(state.target);
    return passed ? 0 : 4;
}
