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
#define DEBUG_GIF_WIDTH 640
#define DEBUG_GIF_HEIGHT 360
#define DEBUG_GIF_HASH_SIZE 8191

typedef struct DebugRunConfig {
    bool requested;
    bool show_window;
    bool record_gif;
    bool capture_steps_explicit;
    int steps;
    int gif_fps;
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
    int tagged_sleeping;
    int sim_particles;
    int collision_particles;
    int collision_particles_culled;
    int static_occupied_cells;
    int static_surface_patches;
    int static_surface_cells;
    int static_surface_patch_refs;
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
    bool sleep_snap_original_test;
    bool sleep_snap_original_succeeded;
    bool sleep_snap_active_observed;
    bool support_removal_succeeded;
    bool drop_sampled;
    bool delayed_activation;
    bool overhang_impact;
    bool pillar_impact_drift_test;
    bool single_corner_hinge_test;
    bool tether_active_then_static_test;
    bool tether_thin_wall_ccd_test;
    bool tether_thin_wall_checked;
    bool tether_thin_wall_tunneled;
    bool tether_thin_wall_hit;
    bool tether_thin_wall_bounced;
    bool tether_thin_wall_wall_woke;
    bool tether_thin_wall_separation_checked;
    bool tether_thin_wall_received_impulse;
    bool tether_second_started;
    bool tether_second_survived;
    bool recycle_test;
    bool recycle_succeeded;
    bool falling_pillar_test;
    bool platform_lifecycle_test;
    bool restoration_deadline_test;
    bool restoration_deadline_hold_applied;
    bool overhang_anchor_observed;
    bool overhang_deflection_sampled;
    bool gpu_readback_performed;
    int activation_step;
    int sleep_snap_static_step;
    int recycle_queue_before;
    int recycle_queue_after;
    float initial_centroid_y;
    float drop_at_step_60;
    float overhang_initial_root_y;
    float overhang_initial_tip_y;
    float overhang_root_drop_at_step_30;
    float overhang_tip_drop_at_step_30;
    float overhang_relative_deflection_at_step_30;
    Vector3 drift_window_centroid;
    float drift_window_distance;
    float drift_window_speed;
    float max_hinge_angle_deg;
    bool hinge_broken;
    int tether_voxel_count_before_static;
    int tether_static_remaining;
    int tether_index_mismatch_step;
    int tether_expected_index;
    int tether_selected_index;
    Vector3 tether_static_start;
    float tether_static_pull_distance;
    uint64_t tether_thin_wall_identity;
    uint64_t tether_thin_wall_hits_before;
    float tether_thin_wall_initial_x;
    float tether_thin_wall_final_x;
    float tether_thin_wall_final_min_x;
    float tether_thin_wall_final_max_x;
    float tether_thin_wall_wall_min_x;
    float tether_thin_wall_wall_max_x;
    float tether_thin_wall_rebound_velocity_x;
    float tether_thin_wall_mean_velocity_x;
    float tether_thin_wall_momentum_x;
    int tether_thin_wall_shared_particles;
    int tether_thin_wall_matching_groups;
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
    char gif_path[DEBUG_PATH_CAPACITY];
    FILE *gif_file;
    int gif_frame_count;
    int gif_step_stride;
    int gif_delay_cs;
    bool gif_failed;
    int lifecycle_sleep_step;
    int lifecycle_recycle_step;
    int restoration_deadline_frame;
} DebugRunState;

typedef bool (*DebugScenarioSetupFn)(DebugRunState *state);

typedef struct DebugScenarioDefinition {
    const char *name;
    const char *description;
    DebugScenarioSetupFn setup;
} DebugScenarioDefinition;

static DebugRunConfig debugRunConfig = {
    .steps = DEBUG_DEFAULT_STEPS,
    .gif_fps = 20
};

static bool debug_setup_dynamic_freefall(DebugRunState *state);
static bool debug_setup_creative_reset_cache(DebugRunState *state);
static bool debug_setup_activation_floating(DebugRunState *state);
static bool debug_setup_sleep_wake_floating(DebugRunState *state);
static bool debug_setup_sleep_snap_original(DebugRunState *state);
static bool debug_setup_sleep_recycle(DebugRunState *state);
static bool debug_setup_falling_pillar_sleep(DebugRunState *state);
static bool debug_setup_three_leg_cut_platform(DebugRunState *state);
static bool debug_setup_tether_throw_floating(DebugRunState *state);
static bool debug_setup_overhang_impact(DebugRunState *state);
static bool debug_setup_pillar_impact_drift(DebugRunState *state);
static bool debug_setup_single_corner_hinge(DebugRunState *state);
static bool debug_setup_tether_active_then_static(DebugRunState *state);
static bool debug_setup_tether_thin_wall_ccd(DebugRunState *state);

static const DebugScenarioDefinition debugScenarios[] = {
    { "dynamic-freefall", "Dynamic cluster freefall control", debug_setup_dynamic_freefall },
    { "creative-reset-cache", "Creative reset invalidates cached static render and collision surfaces", debug_setup_creative_reset_cache },
    { "activation-floating", "Static cluster activated while unsupported", debug_setup_activation_floating },
    { "sleep-wake-floating", "Supported cluster sleeps, loses support, then activates", debug_setup_sleep_wake_floating },
    { "sleep-snap-original", "A calm 4x4x6 pillar near its original pose snaps back to static", debug_setup_sleep_snap_original },
    { "sleep-recycle", "Sleeping dynamic island enters delayed restoration atomically", debug_setup_sleep_recycle },
    { "falling-pillar-sleep", "Rotated pillar falls, settles, and sleeps without changing pose", debug_setup_falling_pillar_sleep },
    { "three-leg-cut-platform", "Four-legged platform with three legs removed completes sleep and recycling", debug_setup_three_leg_cut_platform },
    { "tether-throw-floating", "Released tether voxel wakes an unsupported sleeping cluster", debug_setup_tether_throw_floating },
    { "overhang-impact", "Falling dynamic voxel activates and bends a supported static overhang", debug_setup_overhang_impact },
    { "pillar-impact-drift", "Large undercut pillar impact exposes post-contact sliding and wiggle", debug_setup_pillar_impact_drift },
    { "single-corner-hinge", "A one-corner glue remnant must break beyond the hinge-angle limit", debug_setup_single_corner_hinge },
    { "tether-active-then-static", "Tethering an active voxel then a static voxel preserves and pulls both", debug_setup_tether_active_then_static },
    { "tether-thin-wall-ccd", "High-speed tether throw uses swept particle rays against a one-voxel wall", debug_setup_tether_thin_wall_ccd }
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
    if (strcmp(argument, "--debug-gif") == 0) {
        debugRunConfig.record_gif = true;
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
                               "--debug-gif-fps=", "--debug-gif-fps");
    if (value != (const char *)(uintptr_t)1) {
        if (!debug_parse_integer(value, 1, 60,
                                 &debugRunConfig.gif_fps, "--debug-gif-fps")) return -1;
        debugRunConfig.record_gif = true;
        return 1;
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

static int debug_count_tagged_sleeping(int tag)
{
    int count = 0;
    for (int i = 0; i < voxel_count; ++i) {
        if (voxels[i].debugClusterTag == tag && voxels[i].simulate && voxels[i].sleeping) {
            ++count;
        }
    }
    return count;
}

static void debug_transform_tagged_cluster_off_grid(int tag)
{
    Vector3 pivot = { 0.0f, 0.0f, 0.0f };
    int pivot_count = 0;
    int stamp = ++particle_sync_stamp;
    if (stamp == 0) {
        stamp = 1;
        particle_sync_stamp = 1;
    }
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        if (voxel->debugClusterTag != tag || !voxel->simulate) continue;
        for (int corner = 0; corner < VOXEL_CORNER_COUNT; ++corner) {
            Particle *particle = voxel->particles[corner];
            if (!particle || particle->sync_stamp == stamp) continue;
            particle->sync_stamp = stamp;
            pivot = v_add(pivot, particle->pos);
            ++pivot_count;
        }
    }
    if (pivot_count <= 0) return;
    pivot = v_mul(pivot, 1.0f / (float)pivot_count);

    // Keep the current grid-derived glue graph connected while still proving
    // that the sleep transition preserves a genuinely non-grid pose.
    const float angle = 3.0f * DEG2RAD;
    const float cs = cosf(angle);
    const float sn = sinf(angle);
    const Vector3 offset = { VOXEL_SIZE * 0.18f, 0.0f, VOXEL_SIZE * 0.13f };
    stamp = ++particle_sync_stamp;
    if (stamp == 0) {
        stamp = 1;
        particle_sync_stamp = 1;
    }
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        if (voxel->debugClusterTag != tag || !voxel->simulate) continue;
        for (int corner = 0; corner < VOXEL_CORNER_COUNT; ++corner) {
            Particle *particle = voxel->particles[corner];
            if (!particle || particle->sync_stamp == stamp) continue;
            particle->sync_stamp = stamp;
            Vector3 local = v_sub(particle->pos, pivot);
            Vector3 rotated = {
                cs * local.x - sn * local.z,
                local.y,
                sn * local.x + cs * local.z
            };
            particle->pos = v_add(v_add(pivot, rotated), offset);
            particle->prev_pos = particle->pos;
            particle->predicted_pos = particle->pos;
            particle->vel = (Vector3){ 0.0f, 0.0f, 0.0f };
        }
        Vector3 center = { 0.0f, 0.0f, 0.0f };
        for (int corner = 0; corner < VOXEL_CORNER_COUNT; ++corner) {
            center = v_add(center, voxel->particles[corner]->pos);
        }
        voxel->pos = v_mul(center, 1.0f / (float)VOXEL_CORNER_COUNT);
        voxel->vel = (Vector3){ 0.0f, 0.0f, 0.0f };
    }
    rebuild_voxel_hash();
}

static void debug_place_tagged_cluster_near_original(int tag, Vector3 offset)
{
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        if (voxel->debugClusterTag != tag || !voxel->simulate) continue;
        float min_x = (float)voxel->orig_min_gx * VOXEL_SIZE;
        float min_y = (float)voxel->orig_min_gy * VOXEL_SIZE;
        float min_z = (float)voxel->orig_min_gz * VOXEL_SIZE;
        float max_x = (float)(voxel->orig_max_gx + 1) * VOXEL_SIZE;
        float max_y = (float)(voxel->orig_max_gy + 1) * VOXEL_SIZE;
        float max_z = (float)(voxel->orig_max_gz + 1) * VOXEL_SIZE;
        Vector3 center = { 0.0f, 0.0f, 0.0f };
        for (int corner = 0; corner < VOXEL_CORNER_COUNT; ++corner) {
            Particle *particle = voxel->particles[corner];
            if (!particle) continue;
            Vector3 original = {
                corner_signs[corner][0] < 0 ? min_x : max_x,
                corner_signs[corner][1] < 0 ? min_y : max_y,
                corner_signs[corner][2] < 0 ? min_z : max_z
            };
            particle->pos = v_add(original, offset);
            particle->prev_pos = particle->pos;
            particle->predicted_pos = particle->pos;
            particle->vel = (Vector3){ 0.0f, 0.0f, 0.0f };
            center = v_add(center, particle->pos);
        }
        voxel->pos = v_mul(center, 1.0f / (float)VOXEL_CORNER_COUNT);
        voxel->vel = (Vector3){ 0.0f, 0.0f, 0.0f };
    }
    rebuild_voxel_hash();
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

#if defined(GRAPHICS_API_OPENGL_43) || (defined(__APPLE__) && defined(FPS_GPU_METAL))
static bool debug_gpu_static_cell_absent_in(const GpuStaticCell *cells, int x, int y, int z)
{
    if (!cells) return false;
    uint32_t hash = gpu_static_hash_coord(x, y, z);
    for (int probe = 0; probe < HASH_SIZE; ++probe) {
        const GpuStaticCell *entry = &cells[hash];
        if (entry->value[3] == -1) return true;
        if (entry->value[0] == x && entry->value[1] == y && entry->value[2] == z) {
            // Non-negative values are exterior surface-patch chains indexed in
            // the empty cell beside a solid.  Only <= -2 encodes authoritative
            // static occupancy for embedded-particle recovery.
            return entry->value[3] >= 0;
        }
        hash = (hash + 1u) & (HASH_SIZE - 1u);
    }
    return true;
}

static bool debug_gpu_cluster_cells_are_clear(const DebugRunState *state)
{
    if (!state || state->tracked_cell_count <= 0) return false;
    for (int i = 0; i < state->tracked_cell_count; ++i) {
        const DebugTrackedCell *cell = &state->tracked_cells[i];
        // Some activation scenarios deliberately retain a static anchor.  Only
        // require cells that are absent from authoritative CPU occupancy to be
        // absent from the staged GPU table.
        if (table_get_static_only(cell->x, cell->y, cell->z) >= 0) continue;
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
    if (collisionTopologyDirty) rebuild_particle_collision_metadata();
    (void)ensure_static_surface_cache();
    sample.collision_particles = collision_particle_count;
    sample.collision_particles_culled = sim_particle_count - collision_particle_count;
    sample.static_occupied_cells = staticSurfaceOccupiedCellCount;
    sample.static_surface_patches = patchCount;
    sample.static_surface_cells = staticSurfaceCellCount;
    sample.static_surface_patch_refs = staticSurfacePatchRefCount;
    sample.min_y = FLT_MAX;
    sample.max_y = -FLT_MAX;
    sample.static_hash_generation = staticHashGeneration;
#if defined(GRAPHICS_API_OPENGL_43) || (defined(__APPLE__) && defined(FPS_GPU_METAL))
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
        if (voxel->simulate && voxel->sleeping) ++sample.tagged_sleeping;

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
    DrawText(TextFormat("TAGGED %d dyn=%d static=%d | COLLISION %d/%d",
                        sample->tagged_voxels, sample->tagged_dynamic,
                        sample->tagged_static, sample->collision_particles,
                        sample->sim_particles),
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
            if (state->sleep_snap_original_test) {
                for (int i = 0; i < voxel_count; ++i) {
                    Voxel *voxel = &voxels[i];
                    if (!voxel->simulate || voxel->debugClusterTag != DEBUG_TEST_TAG) continue;
                    for (int corner = 0; corner < VOXEL_CORNER_COUNT; ++corner) {
                        Particle *particle = voxel->particles[corner];
                        if (particle) {
                            DrawSphere(particle->pos, VOXEL_SIZE * 0.08f,
                                       (Color){ 255, 75, 55, 255 });
                        }
                    }
                }
            }
            DrawLine3D((Vector3){ -4.0f, state->initial_centroid_y, 0.0f },
                       (Vector3){ 4.0f, state->initial_centroid_y, 0.0f },
                       (Color){ 255, 225, 80, 210 });
        EndMode3D();
        debug_draw_overlay(state, sample);
    EndTextureMode();
}

typedef struct DebugGifBitWriter {
    FILE *file;
    unsigned char block[255];
    int block_size;
    uint32_t bits;
    int bit_count;
    bool failed;
} DebugGifBitWriter;

static int debugGifHashPrefix[DEBUG_GIF_HASH_SIZE];
static int debugGifHashCode[DEBUG_GIF_HASH_SIZE];
static unsigned char debugGifHashSuffix[DEBUG_GIF_HASH_SIZE];

static void debug_gif_write_u16(FILE *file, int value)
{
    fputc(value & 0xff, file);
    fputc((value >> 8) & 0xff, file);
}

static void debug_gif_flush_block(DebugGifBitWriter *writer)
{
    if (!writer || writer->failed || writer->block_size <= 0) return;
    if (fputc(writer->block_size, writer->file) == EOF ||
        fwrite(writer->block, 1, (size_t)writer->block_size, writer->file) !=
            (size_t)writer->block_size) {
        writer->failed = true;
    }
    writer->block_size = 0;
}

static void debug_gif_write_data_byte(DebugGifBitWriter *writer, unsigned char value)
{
    if (!writer || writer->failed) return;
    writer->block[writer->block_size++] = value;
    if (writer->block_size == 255) debug_gif_flush_block(writer);
}

static void debug_gif_write_code(DebugGifBitWriter *writer, int code, int code_bits)
{
    writer->bits |= ((uint32_t)code) << writer->bit_count;
    writer->bit_count += code_bits;
    while (writer->bit_count >= 8) {
        debug_gif_write_data_byte(writer, (unsigned char)(writer->bits & 0xffu));
        writer->bits >>= 8;
        writer->bit_count -= 8;
    }
}

static void debug_gif_dictionary_reset(void)
{
    memset(debugGifHashCode, 0xff, sizeof(debugGifHashCode));
}

static int debug_gif_dictionary_find(int prefix, unsigned char suffix, int *slot_out)
{
    unsigned int slot = ((unsigned int)prefix * 257u + (unsigned int)suffix * 17u) %
                        DEBUG_GIF_HASH_SIZE;
    for (;;) {
        if (debugGifHashCode[slot] < 0) {
            if (slot_out) *slot_out = (int)slot;
            return -1;
        }
        if (debugGifHashPrefix[slot] == prefix && debugGifHashSuffix[slot] == suffix) {
            if (slot_out) *slot_out = (int)slot;
            return debugGifHashCode[slot];
        }
        slot = (slot + 1u) % DEBUG_GIF_HASH_SIZE;
    }
}

static bool debug_gif_write_lzw(FILE *file, const unsigned char *pixels, size_t count)
{
    if (!file || !pixels || count == 0) return false;
    const int clear_code = 256;
    const int end_code = 257;
    int next_code = 258;
    int code_bits = 9;
    int code_limit = 1 << code_bits;
    DebugGifBitWriter writer = { .file = file };
    fputc(8, file); // minimum LZW code size for the 256-colour palette
    debug_gif_dictionary_reset();
    debug_gif_write_code(&writer, clear_code, code_bits);

    int prefix = pixels[0];
    for (size_t i = 1; i < count; ++i) {
        unsigned char suffix = pixels[i];
        int slot = -1;
        int found = debug_gif_dictionary_find(prefix, suffix, &slot);
        if (found >= 0) {
            prefix = found;
            continue;
        }

        debug_gif_write_code(&writer, prefix, code_bits);
        if (next_code < 4096) {
            debugGifHashPrefix[slot] = prefix;
            debugGifHashSuffix[slot] = suffix;
            debugGifHashCode[slot] = next_code++;
            // The decoder learns a dictionary entry one emitted code after the
            // encoder, so widen only after crossing the current code boundary.
            if (next_code > code_limit && code_bits < 12) {
                ++code_bits;
                code_limit <<= 1;
            }
        } else {
            debug_gif_write_code(&writer, clear_code, code_bits);
            debug_gif_dictionary_reset();
            next_code = 258;
            code_bits = 9;
            code_limit = 1 << code_bits;
        }
        prefix = suffix;
    }
    debug_gif_write_code(&writer, prefix, code_bits);
    debug_gif_write_code(&writer, end_code, code_bits);
    if (writer.bit_count > 0) {
        debug_gif_write_data_byte(&writer, (unsigned char)(writer.bits & 0xffu));
    }
    debug_gif_flush_block(&writer);
    fputc(0, file); // end of image data sub-blocks
    return !writer.failed && !ferror(file);
}

static bool debug_gif_begin(DebugRunState *state)
{
    if (!debugRunConfig.record_gif) return true;
    if (!debug_join_path(state->gif_path, sizeof(state->gif_path),
                         state->output_dir, "simulation.gif")) {
        debug_add_failure(state, "GIF output path is too long");
        return false;
    }
    state->gif_file = fopen(state->gif_path, "wb");
    if (!state->gif_file) {
        debug_add_failure(state, "unable to create GIF %s: %s",
                          state->gif_path, strerror(errno));
        return false;
    }

    float source_fps = 1.0f / (PBD_MAX_STEP_DT);
    state->gif_step_stride = (int)roundf(source_fps / (float)debugRunConfig.gif_fps);
    if (state->gif_step_stride < 1) state->gif_step_stride = 1;
    state->gif_delay_cs = (int)roundf(100.0f * PBD_MAX_STEP_DT *
                                     (float)state->gif_step_stride);
    if (state->gif_delay_cs < 1) state->gif_delay_cs = 1;

    fwrite("GIF89a", 1, 6, state->gif_file);
    debug_gif_write_u16(state->gif_file, DEBUG_GIF_WIDTH);
    debug_gif_write_u16(state->gif_file, DEBUG_GIF_HEIGHT);
    fputc(0xf7, state->gif_file); // global 256-entry colour table
    fputc(0, state->gif_file);
    fputc(0, state->gif_file);
    for (int i = 0; i < 256; ++i) {
        fputc(((i >> 5) & 7) * 255 / 7, state->gif_file);
        fputc(((i >> 2) & 7) * 255 / 7, state->gif_file);
        fputc((i & 3) * 255 / 3, state->gif_file);
    }
    const unsigned char loop_extension[] = {
        0x21, 0xff, 0x0b, 'N','E','T','S','C','A','P','E','2','.','0',
        0x03, 0x01, 0x00, 0x00, 0x00
    };
    fwrite(loop_extension, 1, sizeof(loop_extension), state->gif_file);
    if (ferror(state->gif_file)) {
        debug_add_failure(state, "unable to write GIF header %s", state->gif_path);
        fclose(state->gif_file);
        state->gif_file = NULL;
        return false;
    }
    return true;
}

static bool debug_gif_capture(DebugRunState *state, int step)
{
    if (!state || !state->gif_file || state->gif_failed) {
        return false;
    }
    if (step != 0 && step != debugRunConfig.steps &&
        (step % state->gif_step_stride) != 0) return true;

    DebugSample sample = debug_collect_sample("animation", step, state->initial_centroid_y);
    debug_render_scene(state, &sample);
    Image frame = LoadImageFromTexture(state->target.texture);
    if (!IsImageValid(frame)) {
        debug_add_failure(state, "unable to read GIF frame at step %d", step);
        state->gif_failed = true;
        return false;
    }
    ImageFlipVertical(&frame);
    ImageResize(&frame, DEBUG_GIF_WIDTH, DEBUG_GIF_HEIGHT);
    ImageFormat(&frame, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
    unsigned char *rgba = (unsigned char *)frame.data;
    size_t pixel_count = (size_t)DEBUG_GIF_WIDTH * (size_t)DEBUG_GIF_HEIGHT;
    unsigned char *indexed = (unsigned char *)malloc(pixel_count);
    if (!indexed) {
        UnloadImage(frame);
        debug_add_failure(state, "unable to allocate GIF frame buffer");
        state->gif_failed = true;
        return false;
    }
    for (size_t i = 0; i < pixel_count; ++i) {
        indexed[i] = (unsigned char)((rgba[i * 4 + 0] & 0xe0u) |
                                     ((rgba[i * 4 + 1] & 0xe0u) >> 3) |
                                     (rgba[i * 4 + 2] >> 6));
    }
    UnloadImage(frame);

    fputc(0x21, state->gif_file); fputc(0xf9, state->gif_file); fputc(0x04, state->gif_file);
    fputc(0x00, state->gif_file);
    debug_gif_write_u16(state->gif_file, state->gif_delay_cs);
    fputc(0x00, state->gif_file); fputc(0x00, state->gif_file);
    fputc(0x2c, state->gif_file);
    debug_gif_write_u16(state->gif_file, 0); debug_gif_write_u16(state->gif_file, 0);
    debug_gif_write_u16(state->gif_file, DEBUG_GIF_WIDTH);
    debug_gif_write_u16(state->gif_file, DEBUG_GIF_HEIGHT);
    fputc(0x00, state->gif_file);
    bool written = debug_gif_write_lzw(state->gif_file, indexed, pixel_count);
    free(indexed);
    if (!written) {
        debug_add_failure(state, "unable to encode GIF frame at step %d", step);
        state->gif_failed = true;
        return false;
    }
    ++state->gif_frame_count;
    return true;
}

static bool debug_gif_end(DebugRunState *state)
{
    if (!state || !state->gif_file) return !debugRunConfig.record_gif || state->gif_failed;
    fputc(0x3b, state->gif_file);
    bool ok = !state->gif_failed && !ferror(state->gif_file) && fclose(state->gif_file) == 0;
    state->gif_file = NULL;
    if (!ok && !state->gif_failed) {
        debug_add_failure(state, "unable to finalize GIF %s", state->gif_path);
    }
    return ok;
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
    if (sample->tagged_sleeping != 0) {
        debug_add_failure(state, "%s: %d tagged voxels remain sleeping", context,
                          sample->tagged_sleeping);
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
    if (!physics_backend_is_gpu(physicsBackend.active) || step < 1 ||
        !state || !state->activation_succeeded || state->gpu_readback_performed) return;
#if defined(GRAPHICS_API_OPENGL_43) || (defined(__APPLE__) && defined(FPS_GPU_METAL))
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
    if (!state || !physics_backend_is_gpu(physicsBackend.active) || step < 1 ||
        !state->activation_succeeded || state->gpu_readback_performed) return;
#if defined(GRAPHICS_API_OPENGL_43) || (defined(__APPLE__) && defined(FPS_GPU_METAL))
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

    gpu_buffer_read(&gpuPhysics.static_cell_buffer, static_readback, static_bytes, 0);
    if (sim_bytes) {
        gpu_buffer_read(&gpuPhysics.sim_id_buffer, sim_id_readback, sim_bytes, 0);
    }
    if (particle_bytes) {
        gpu_buffer_read(&gpuPhysics.particle_buffer, particle_readback, particle_bytes, 0);
    }

    for (int i = 0; i < HASH_SIZE; ++i) {
        if (memcmp(&static_readback[i], &gpuPhysics.static_stage[i],
                   sizeof(static_readback[i])) != 0) {
            ++state->gpu_static_ssbo_diff_entries;
        }
    }
    for (int i = 0; i < state->tracked_cell_count; ++i) {
        const DebugTrackedCell *cell = &state->tracked_cells[i];
        if (table_get_static_only(cell->x, cell->y, cell->z) >= 0) continue;
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
    int sleeping_seed = -1;
    for (int i = 0; i < voxel_count; ++i) {
        if (voxels[i].debugClusterTag == DEBUG_TEST_TAG &&
            voxels[i].simulate && voxels[i].sleeping) {
            sleeping_seed = i;
            break;
        }
    }
    if (sleeping_seed >= 0) {
        bool activated = wake_sleeping_cluster_in_place(sleeping_seed);
        state->activation_succeeded = activated;
        if (!activated) {
            debug_add_failure(state, "unable to wake sleeping test cluster in place");
        }
        return activated;
    }

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

static bool debug_setup_sleep_snap_original(DebugRunState *state)
{
    debug_reset_world();
    state->sleep_snap_original_test = true;
    state->sleep_snap_static_step = -1;
    state->expected_voxels = 4 * 4 * 6;
    int cluster_created = 0;
    Color pillar_color = (Color){ 80, 175, 245, 255 };
    for (int y = 0; y < 6; ++y) {
        for (int z = -2; z < 2; ++z) {
            for (int x = -2; x < 2; ++x) {
                if (debug_add_cell(x, y, z, true, false,
                                   pillar_color, DEBUG_TEST_TAG) >= 0) {
                    ++cluster_created;
                }
            }
        }
    }
    debug_track_tagged_cells(state, DEBUG_TEST_TAG);
    debug_rebuild_world_state();

    if (cluster_created != state->expected_voxels) {
        debug_add_failure(state, "sleep-snap setup created pillar=%d/%d",
                          cluster_created, state->expected_voxels);
        return false;
    }
    state->camera.position = (Vector3){ 6.5f, 4.5f, 8.0f };
    state->camera.target = (Vector3){ 0.0f, 1.5f, 0.0f };
    debug_set_initial_centroid(state);
    debug_capture(state, "phase-01-original-static.png", "original-static", 0);
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
    debug_transform_tagged_cluster_off_grid(DEBUG_TEST_TAG);
    debug_capture(state, "phase-01-supported-dynamic.png", "supported-off-grid-dynamic", 0);

    Vector3 *positions_before_sleep = malloc((size_t)particle_pool_count * sizeof(*positions_before_sleep));
    if (!positions_before_sleep) {
        debug_add_failure(state, "unable to allocate no-snap pose snapshot");
        return false;
    }
    for (int i = 0; i < particle_pool_count; ++i) {
        positions_before_sleep[i] = particles_pool[i].pos;
    }

    for (int i = 0; i < voxel_count; ++i) {
        if (voxels[i].debugClusterTag == DEBUG_TEST_TAG && voxels[i].simulate) {
            voxels[i].sleepFrames = VOXEL_DEACTIVATION_FRAMES - 1;
            set_voxel_velocity(&voxels[i], (Vector3){ 0.0f, 0.0f, 0.0f });
        }
    }
    state->sleeping_conversion_succeeded = deactivate_sleeping_voxels();
    rebuild_static_hash_if_dirty();
    if (!state->sleeping_conversion_succeeded ||
        debug_count_tagged_sleeping(DEBUG_TEST_TAG) != state->expected_voxels) {
        debug_add_failure(state, "sleep conversion did not produce %d frozen dynamic voxels",
                          state->expected_voxels);
        free(positions_before_sleep);
        return false;
    }
    float max_sleep_pose_delta = 0.0f;
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        if (voxel->debugClusterTag != DEBUG_TEST_TAG) continue;
        for (int corner = 0; corner < VOXEL_CORNER_COUNT; ++corner) {
            Particle *particle = voxel->particles[corner];
            ptrdiff_t pool_id = particle ? particle - particles_pool : -1;
            if (pool_id < 0 || pool_id >= particle_pool_count) continue;
            float delta = v_length(v_sub(particle->pos, positions_before_sleep[pool_id]));
            if (delta > max_sleep_pose_delta) max_sleep_pose_delta = delta;
        }
    }
    free(positions_before_sleep);
    if (max_sleep_pose_delta > 1e-6f) {
        debug_add_failure(state, "sleep transition changed the preserved pose by %.9g metres",
                          (double)max_sleep_pose_delta);
        return false;
    }
    debug_capture(state, "phase-02-sleeping-dynamic.png", "sleeping-dynamic", 0);

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

static bool debug_setup_creative_reset_cache(DebugRunState *state)
{
    debug_reset_world();
    int created = debug_add_test_cluster(false);
    debug_rebuild_world_state();
    if (created != state->expected_voxels || !ensure_static_surface_cache() ||
        patchCount <= 0 || staticSurfaceOccupiedCellCount <= 0) {
        debug_add_failure(state, "unable to populate the pre-reset static surface cache");
        return false;
    }

    uint64_t populated_generation = staticSurfaceGeneration;
    ResetCreative();
    if (!ensure_static_surface_cache()) {
        debug_add_failure(state, "creative reset could not rebuild the empty surface cache");
        return false;
    }
    if (staticSurfaceGeneration == populated_generation) {
        debug_add_failure(state, "creative reset did not advance the static surface generation");
        return false;
    }
    if (voxel_count != 0 || patchCount != 0 || staticSurfaceOccupiedCellCount != 0 ||
        staticSurfaceCellCount != 0 || staticSurfacePatchRefCount != 0) {
        debug_add_failure(state,
                          "creative reset retained stale geometry: voxels=%d patches=%d occupied=%d cells=%d refs=%d",
                          voxel_count, patchCount, staticSurfaceOccupiedCellCount,
                          staticSurfaceCellCount, staticSurfacePatchRefCount);
        return false;
    }

    // Finish as a normal deterministic simulation so the scenario also checks
    // that rebuilding the world after the empty creative reset remains valid.
    return debug_setup_dynamic_freefall(state);
}

static bool debug_setup_sleep_recycle(DebugRunState *state)
{
    debug_reset_world();
    state->recycle_test = true;
    int support_created = debug_add_support_column();
    debug_rebuild_world_state();
    int cluster_created = debug_add_test_cluster(true);
    debug_track_tagged_cells(state, DEBUG_TEST_TAG);
    debug_rebuild_world_state();

    if (support_created != 12 || cluster_created != state->expected_voxels) {
        debug_add_failure(state, "sleep-recycle setup created support=%d/12 cluster=%d/%d",
                          support_created, cluster_created, state->expected_voxels);
        return false;
    }

    debug_transform_tagged_cluster_off_grid(DEBUG_TEST_TAG);
    for (int i = 0; i < voxel_count; ++i) {
        if (voxels[i].debugClusterTag == DEBUG_TEST_TAG && voxels[i].simulate) {
            voxels[i].sleepFrames = VOXEL_DEACTIVATION_FRAMES - 1;
            set_voxel_velocity(&voxels[i], (Vector3){ 0.0f, 0.0f, 0.0f });
        }
    }
    state->sleeping_conversion_succeeded = deactivate_sleeping_voxels();
    if (!state->sleeping_conversion_succeeded ||
        debug_count_tagged_sleeping(DEBUG_TEST_TAG) != state->expected_voxels) {
        debug_add_failure(state, "sleep-recycle did not freeze %d dynamic voxels",
                          state->expected_voxels);
        return false;
    }
    debug_capture(state, "phase-01-sleeping-island.png", "sleeping-before-recycle", 0);

    for (int i = 0; i < voxel_count; ++i) {
        if (voxels[i].debugClusterTag == DEBUG_TEST_TAG && voxels[i].sleeping) {
            voxels[i].lifeFrames = RECYCLE_DYNAMIC_MAX_FRAMES;
        }
    }
    state->recycle_queue_before = recycleQueueCount;
    recycle_dead_voxels();
    state->recycle_queue_after = recycleQueueCount;
    state->recycle_succeeded =
        debug_count_tagged(DEBUG_TEST_TAG, -1) == 0 &&
        (state->recycle_queue_after - state->recycle_queue_before) == state->expected_voxels;
    if (!state->recycle_succeeded) {
        debug_add_failure(state,
                          "sleep-recycle left tagged=%d and queued=%d/%d snapshots",
                          debug_count_tagged(DEBUG_TEST_TAG, -1),
                          state->recycle_queue_after - state->recycle_queue_before,
                          state->expected_voxels);
        return false;
    }
    debug_capture(state, "phase-02-recycled-island.png", "after-atomic-recycle", 0);
    return true;
}

static bool debug_setup_falling_pillar_sleep(DebugRunState *state)
{
    debug_reset_world();
    state->falling_pillar_test = true;
    state->recycle_queue_before = recycleQueueCount;
    state->lifecycle_sleep_step = -1;
    state->lifecycle_recycle_step = -1;
    const int min_y = 1;
    const int height = 10;
    int created = 0;
    Color pillar_color = (Color){ 225, 155, 65, 255 };
    for (int y = min_y; y < min_y + height; ++y) {
        for (int z = -1; z <= 1; ++z) {
            for (int x = -1; x <= 1; ++x) {
                if (debug_add_cell(x, y, z, false, true,
                                   pillar_color, DEBUG_TEST_TAG) >= 0) {
                    ++created;
                }
            }
        }
    }
    state->expected_voxels = created;
    debug_track_tagged_cells(state, DEBUG_TEST_TAG);
    debug_rebuild_world_state();

    if (created != 3 * height * 3) {
        debug_add_failure(state, "falling-pillar setup created %d/%d voxels",
                          created, 3 * height * 3);
        return false;
    }

    Vector3 pivot = { 0.0f, 0.0f, 0.0f };
    int pivot_count = 0;
    int stamp = ++particle_sync_stamp;
    if (stamp == 0) {
        stamp = 1;
        particle_sync_stamp = 1;
    }
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        if (voxel->debugClusterTag != DEBUG_TEST_TAG) continue;
        for (int corner = 0; corner < VOXEL_CORNER_COUNT; ++corner) {
            Particle *particle = voxel->particles[corner];
            if (!particle || particle->sync_stamp == stamp) continue;
            particle->sync_stamp = stamp;
            pivot = v_add(pivot, particle->pos);
            ++pivot_count;
        }
    }
    if (pivot_count <= 0) {
        debug_add_failure(state, "falling-pillar setup found no particles");
        return false;
    }
    pivot = v_mul(pivot, 1.0f / (float)pivot_count);

    const float angle = 0.0f;
    const float cs = cosf(angle);
    const float sn = sinf(angle);
    const float angular_speed = 0.75f;
    stamp = ++particle_sync_stamp;
    if (stamp == 0) {
        stamp = 1;
        particle_sync_stamp = 1;
    }
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        if (voxel->debugClusterTag != DEBUG_TEST_TAG) continue;
        for (int corner = 0; corner < VOXEL_CORNER_COUNT; ++corner) {
            Particle *particle = voxel->particles[corner];
            if (!particle || particle->sync_stamp == stamp) continue;
            particle->sync_stamp = stamp;
            Vector3 local = v_sub(particle->pos, pivot);
            Vector3 rotated = {
                cs * local.x - sn * local.y,
                sn * local.x + cs * local.y,
                local.z
            };
            particle->pos = v_add(pivot, rotated);
            particle->prev_pos = particle->pos;
            particle->predicted_pos = particle->pos;
            particle->vel = (Vector3){
                -angular_speed * rotated.y,
                 angular_speed * rotated.x,
                 0.0f
            };
        }
        Vector3 center = { 0.0f, 0.0f, 0.0f };
        Vector3 velocity = { 0.0f, 0.0f, 0.0f };
        for (int corner = 0; corner < VOXEL_CORNER_COUNT; ++corner) {
            center = v_add(center, voxel->particles[corner]->pos);
            velocity = v_add(velocity, voxel->particles[corner]->vel);
        }
        voxel->pos = v_mul(center, 1.0f / (float)VOXEL_CORNER_COUNT);
        voxel->vel = v_mul(velocity, 1.0f / (float)VOXEL_CORNER_COUNT);
        voxel->wake_timer = 0;
    }
    // Model a late activation joining an older unsupported island.  The young
    // member must not extend the maximum awake lifetime of the whole structure.
    for (int i = 0; i < voxel_count; ++i) {
        if (voxels[i].debugClusterTag == DEBUG_TEST_TAG) {
            voxels[i].lifeFrames = -RECYCLE_DYNAMIC_MAX_FRAMES;
            break;
        }
    }
    rebuild_voxel_hash();
    state->activation_succeeded = true;
    debug_set_initial_centroid(state);
    debug_capture(state, "phase-01-elevated-rotated-pillar.png",
                  "elevated-rotated-pillar", 0);
    return true;
}

static bool debug_setup_three_leg_cut_platform(DebugRunState *state)
{
    debug_reset_world();
    state->platform_lifecycle_test = true;
    state->restoration_deadline_test = true;
    state->recycle_queue_before = recycleQueueCount;
    state->lifecycle_sleep_step = -1;
    state->lifecycle_recycle_step = -1;
    Color deck_color = (Color){ 90, 220, 150, 255 };
    Color leg_color = (Color){ 200, 210, 70, 255 };
    int created = 0;

    // Exact gameplay construction: one-cell 12x12 deck, three underside ribs,
    // and only the southwest notched 2x2 leg after the other three are cut.
    for (int z = -6; z <= 5; ++z) {
        for (int x = -6; x <= 5; ++x) {
            if (debug_add_cell(x, 5, z, true, false,
                               deck_color, DEBUG_TEST_TAG) >= 0) ++created;
        }
    }
    for (int x = -6; x <= 5; ++x) {
        if (x == 0) continue;
        for (int z = -2; z <= 2; z += 2) {
            if (debug_add_cell(x, 4, z, true, false,
                               deck_color, DEBUG_TEST_TAG) >= 0) ++created;
        }
    }
    for (int y = 0; y <= 4; ++y) {
        for (int z = -6; z <= -5; ++z) {
            for (int x = -6; x <= -5; ++x) {
                if (y == 4 && x == -5 && z == -5) continue;
                if (debug_add_cell(x, y, z, true, false,
                                   leg_color, DEBUG_TEST_TAG) >= 0) ++created;
            }
        }
    }
    state->expected_voxels = created;
    debug_track_tagged_cells(state, DEBUG_TEST_TAG);
    debug_rebuild_world_state();
    debug_set_initial_centroid(state);
    debug_capture(state, "phase-01-three-legs-cut.png", "three-legs-cut-static", 0);
    if (created != 196) {
        debug_add_failure(state, "three-leg platform created %d/196 voxels", created);
        return false;
    }
    if (!debug_activate_test_cluster(state)) return false;
    state->restoration_deadline_frame =
        recycleFrameCounter + RECYCLE_STATIC_RESTORE_DELAY;
    int dynamic_count = debug_count_tagged(DEBUG_TEST_TAG, 1);
    if (dynamic_count != state->expected_voxels) {
        debug_add_failure(state, "three-leg platform activated %d/%d voxels",
                          dynamic_count, state->expected_voxels);
        return false;
    }
    debug_set_initial_centroid(state);
    state->camera.position = (Vector3){ 7.0f, 6.0f, 9.0f };
    state->camera.target = (Vector3){ 0.0f, 1.0f, 0.0f };
    debug_capture(state, "phase-02-platform-activated.png", "platform-activated", 0);
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
        debug_count_tagged_sleeping(DEBUG_TEST_TAG) != state->expected_voxels) {
        debug_add_failure(state, "sleep conversion did not produce %d frozen dynamic voxels",
                          state->expected_voxels);
        return false;
    }
    debug_capture(state, "phase-02-sleeping-dynamic.png", "sleeping-dynamic", 0);

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
        players[i].tetherVoxelIdentity = 0;
    }
    players[0].yaw = -90.0f;
    players[0].pitch = 0.0f;
    players[0].tetherHolding = true;
    players[0].tetherVoxel = thrown_idx;
    players[0].tetherVoxelIdentity = voxels[thrown_idx].identity;
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

static bool debug_tagged_centroid_and_horizontal_speed(int tag, Vector3 *centroid_out,
                                                        float *speed_out)
{
    Vector3 centroid = { 0.0f, 0.0f, 0.0f };
    float speed_sum = 0.0f;
    int count = 0;
    for (int i = 0; i < voxel_count; ++i) {
        const Voxel *voxel = &voxels[i];
        if (voxel->debugClusterTag != tag || !voxel->simulate) continue;
        centroid = v_add(centroid, voxel->pos);
        speed_sum += sqrtf(voxel->vel.x * voxel->vel.x + voxel->vel.z * voxel->vel.z);
        ++count;
    }
    if (count <= 0) return false;
    if (centroid_out) *centroid_out = v_mul(centroid, 1.0f / (float)count);
    if (speed_out) *speed_out = speed_sum / (float)count;
    return true;
}

static Vector3 debug_voxel_axis(const Voxel *voxel, int axis)
{
    Vector3 result = { 0.0f, 0.0f, 0.0f };
    if (!voxel || axis < 0 || axis > 2) return result;
    const int pairs[3][4][2] = {
        { {0,1}, {2,3}, {4,5}, {6,7} },
        { {0,2}, {1,3}, {4,6}, {5,7} },
        { {0,4}, {1,5}, {2,6}, {3,7} }
    };
    for (int pair = 0; pair < 4; ++pair) {
        Particle *a = voxel->particles[pairs[axis][pair][0]];
        Particle *b = voxel->particles[pairs[axis][pair][1]];
        if (a && b) result = v_add(result, v_sub(b->pos, a->pos));
    }
    return v_norm(result);
}

static float debug_voxel_relative_angle_deg(const Voxel *a, const Voxel *b)
{
    float minimum_dot = 1.0f;
    for (int axis = 0; axis < 3; ++axis) {
        float dot = clampf(v_dot(debug_voxel_axis(a, axis), debug_voxel_axis(b, axis)),
                           -1.0f, 1.0f);
        if (dot < minimum_dot) minimum_dot = dot;
    }
    return acosf(minimum_dot) / DEG2RAD;
}

static bool debug_setup_pillar_impact_drift(DebugRunState *state)
{
    debug_reset_world();
    state->pillar_impact_drift_test = true;
    Color color = (Color){ 225, 145, 55, 255 };
    int created = 0;
    for (int y = 12; y < 32; ++y) {
        for (int z = -2; z <= 2; ++z) {
            for (int x = -2; x <= 2; ++x) {
                if (debug_add_cell(x, y, z, true, false, color, DEBUG_TEST_TAG) >= 0) {
                    ++created;
                }
            }
        }
    }
    state->expected_voxels = created;
    debug_track_tagged_cells(state, DEBUG_TEST_TAG);
    debug_rebuild_world_state();
    if (created != 500 || !debug_activate_test_cluster(state)) {
        debug_add_failure(state, "pillar-impact setup activated %d/500 voxels", created);
        return false;
    }

    Vector3 pivot = { 0.0f, 3.0f, 0.0f };
    const float angle = 4.0f * DEG2RAD;
    const float cs = cosf(angle), sn = sinf(angle);
    const float angular_speed = 0.45f;
    int stamp = ++particle_sync_stamp;
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        if (voxel->debugClusterTag != DEBUG_TEST_TAG || !voxel->simulate) continue;
        for (int corner = 0; corner < VOXEL_CORNER_COUNT; ++corner) {
            Particle *particle = voxel->particles[corner];
            if (!particle || particle->sync_stamp == stamp) continue;
            particle->sync_stamp = stamp;
            Vector3 local = v_sub(particle->pos, pivot);
            Vector3 rotated = { cs * local.x - sn * local.y,
                                sn * local.x + cs * local.y, local.z };
            particle->pos = v_add(pivot, rotated);
            particle->prev_pos = particle->pos;
            particle->predicted_pos = particle->pos;
            particle->vel = (Vector3){ -angular_speed * rotated.y,
                                       angular_speed * rotated.x, 0.0f };
        }
        Vector3 center = { 0.0f, 0.0f, 0.0f };
        Vector3 velocity = { 0.0f, 0.0f, 0.0f };
        for (int corner = 0; corner < VOXEL_CORNER_COUNT; ++corner) {
            center = v_add(center, voxel->particles[corner]->pos);
            velocity = v_add(velocity, voxel->particles[corner]->vel);
        }
        voxel->pos = v_mul(center, 0.125f);
        voxel->vel = v_mul(velocity, 0.125f);
    }
    collisionTopologyDirty = true;
    rebuild_voxel_hash();
    debug_set_initial_centroid(state);
    state->camera.position = (Vector3){ 10.0f, 6.0f, 13.0f };
    state->camera.target = (Vector3){ 0.0f, 3.5f, 0.0f };
    debug_capture(state, "phase-01-undercut-pillar.png", "undercut-pillar", 0);
    return true;
}

static bool debug_setup_single_corner_hinge(DebugRunState *state)
{
    debug_reset_world();
    state->single_corner_hinge_test = true;
    int lower_idx = debug_add_cell(0, 8, 0, true, true,
                                   (Color){ 105, 120, 145, 255 }, DEBUG_SUPPORT_TAG);
    int upper_idx = debug_add_cell(0, 9, 0, false, true,
                                   (Color){ 245, 150, 55, 255 }, DEBUG_TEST_TAG);
    if (lower_idx < 0 || upper_idx < 0) {
        debug_add_failure(state, "single-corner hinge could not create its two voxels");
        return false;
    }
    debug_rebuild_world_state();
    Voxel *lower = &voxels[lower_idx];
    Voxel *upper = &voxels[upper_idx];
    const int upper_face = 3;
    const int lower_face = 2;
    int pivot_corner = face_corner_indices[upper_face][0];
    Particle *pivot_particle = upper->particles[pivot_corner];
    for (int c = 1; c < 4; ++c) {
        int corner = face_corner_indices[upper_face][c];
        Particle *old = upper->particles[corner];
        Particle *clone = particle_clone(old);
        if (!clone) {
            debug_add_failure(state, "single-corner hinge could not split face corner %d", c);
            return false;
        }
        clone->inv_mass = 1.0f;
        clone->base_inv_mass = 1.0f;
        if (clone->sim_index < 0) sim_particles_add(clone);
        upper->particles[corner] = clone;
        particle_release(old);
    }
    upper->glued_faces[upper_face] = true;
    lower->glued_faces[lower_face] = true;

    const float angle = 75.0f * DEG2RAD;
    const float cs = cosf(angle), sn = sinf(angle);
    Vector3 pivot = pivot_particle->pos;
    for (int corner = 0; corner < VOXEL_CORNER_COUNT; ++corner) {
        Particle *particle = upper->particles[corner];
        if (!particle || particle == pivot_particle) continue;
        Vector3 local = v_sub(particle->pos, pivot);
        Vector3 rotated = { cs * local.x - sn * local.y,
                            sn * local.x + cs * local.y, local.z };
        particle->pos = v_add(pivot, rotated);
        particle->prev_pos = particle->pos;
        particle->predicted_pos = particle->pos;
        particle->vel = (Vector3){ 0.0f, 0.0f, 0.0f };
    }
    Vector3 center = { 0.0f, 0.0f, 0.0f };
    for (int corner = 0; corner < VOXEL_CORNER_COUNT; ++corner) {
        center = v_add(center, upper->particles[corner]->pos);
    }
    upper->pos = v_mul(center, 0.125f);
    upper->vel = (Vector3){ 0.0f, 0.0f, 0.0f };
    state->expected_voxels = 1;
    state->activation_succeeded = true;
    state->max_hinge_angle_deg = debug_voxel_relative_angle_deg(lower, upper);
    collisionTopologyDirty = true;
    rebuild_voxel_hash();
    debug_set_initial_centroid(state);
    state->camera.position = (Vector3){ 2.2f, 2.6f, 3.0f };
    state->camera.target = (Vector3){ 0.15f, 2.25f, 0.15f };
    debug_capture(state, "phase-01-single-corner-hinge.png", "single-corner-hinge", 0);
    return true;
}

static bool debug_setup_tether_active_then_static(DebugRunState *state)
{
    debug_reset_world();
    state->tether_active_then_static_test = true;
    state->tether_index_mismatch_step = -1;
    Color static_color = (Color){ 75, 175, 245, 255 };
    for (int z = -1; z <= 1; ++z) {
        if (debug_add_cell(8, 8, z, true, false, static_color, DEBUG_SUPPORT_TAG) < 0) {
            debug_add_failure(state, "tether sequence could not create static target %d", z);
            return false;
        }
    }
    int active_idx = debug_add_cell(0, 8, 0, false, true,
                                    (Color){ 245, 145, 45, 255 }, DEBUG_THROW_TAG);
    if (active_idx < 0) {
        debug_add_failure(state, "tether sequence could not create active target");
        return false;
    }
    debug_track_tagged_cells(state, DEBUG_SUPPORT_TAG);
    debug_rebuild_world_state();

    memset(tetherTag, 0, sizeof(tetherTag));
    for (int i = 0; i < MAX_PLAYERS; ++i) {
        players[i].tetherHolding = false;
        players[i].tetherVoxel = -1;
        players[i].tetherVoxelIdentity = 0;
        players[i].respawn_timer = (i == 0) ? 0.0f : 1000.0f;
    }
    players[0].pos = (Vector3){ 0.0f, 2.0f, 2.0f };
    Vector3 first_aim = v_norm(v_sub(voxels[active_idx].pos, players[0].pos));
    dir_to_yaw_pitch(first_aim, &players[0].yaw, &players[0].pitch);
    int first_ray_hit = first_voxel_hit((Ray){ players[0].pos, player_forward(&players[0]) },
                                        TETHER_RANGE, -1);
    start_tether(0);
    if (!players[0].tetherHolding || players[0].tetherVoxel != active_idx) {
        debug_add_failure(state,
                          "first tether did not select the active voxel (ray=%d selected=%d expected=%d table=%d)",
                          first_ray_hit, players[0].tetherVoxel, active_idx,
                          table_get(0, 8, 0));
        return false;
    }
    prepare_tether_forces();
    state->expected_voxels = 0;
    state->activation_succeeded = true;
    debug_set_initial_centroid(state);
    state->camera.position = (Vector3){ 4.0f, 3.4f, 5.0f };
    state->camera.target = (Vector3){ 1.0f, 2.1f, 0.0f };
    debug_capture(state, "phase-01-active-tethered.png", "active-tethered", 0);
    return true;
}

static bool debug_setup_tether_thin_wall_ccd(DebugRunState *state)
{
    debug_reset_world();
    tetherThrowCcdHitCount = 0;
    state->tether_thin_wall_ccd_test = true;
    state->tether_thin_wall_wall_min_x = 0.0f;
    state->tether_thin_wall_wall_max_x = VOXEL_SIZE;

    Color wall_color = (Color){ 105, 125, 155, 255 };
    int wall_created = 0;
    // One grid cell thick along the throw axis, but broad enough that no corner
    // can miss around an edge.  Its front face is x=0.
    for (int y = 1; y <= 7; ++y) {
        for (int z = -3; z <= 3; ++z) {
            if (debug_add_cell(0, y, z, true, false,
                               wall_color, DEBUG_SUPPORT_TAG) >= 0) {
                ++wall_created;
            }
        }
    }
    int thrown_index = debug_add_cell(-3, 4, 0, false, true,
                                      (Color){ 255, 195, 55, 255 }, DEBUG_TEST_TAG);
    if (wall_created != 49 || thrown_index < 0) {
        debug_add_failure(state,
                          "thin-wall setup created wall=%d/49 thrown=%d",
                          wall_created, thrown_index);
        return false;
    }
    debug_rebuild_world_state();

    Voxel *thrown = &voxels[thrown_index];
    thrown->owner = 0;
    thrown->activator = 0;
    thrown->wasTethered = true;
    thrown->activationBelief = 1.0f;
    // A strongly accelerated tether can already carry velocity when release
    // adds its normal impulse.  300 m/s guarantees a one-cell wall lies fully
    // between two endpoint-only collision samples.
    set_voxel_velocity(thrown, (Vector3){ 300.0f, 0.0f, 0.0f });
    state->tether_thin_wall_identity = thrown->identity;
    state->tether_thin_wall_initial_x = thrown->pos.x;
    state->tether_thin_wall_hits_before = tetherThrowCcdHitCount;

    players[0].yaw = -90.0f;
    players[0].pitch = 0.0f;
    players[0].tetherHolding = true;
    players[0].tetherVoxel = thrown_index;
    players[0].tetherVoxelIdentity = thrown->identity;
    release_tether(0);

    if (tetherThrowCcdEnabled && thrown->tetherThrowCcdFrames == 0) {
        debug_add_failure(state, "released high-speed tether voxel did not enter CCD tracking");
        return false;
    }
    state->expected_voxels = 1;
    state->activation_succeeded = true;
    debug_set_initial_centroid(state);
    state->camera.position = (Vector3){ -4.5f, 4.2f, 7.0f };
    state->camera.target = (Vector3){ 0.0f, 2.2f, 0.0f };
    debug_capture(state, "phase-01-tether-released.png", "tether-released", 0);
    return true;
}

static void debug_physics_step(void)
{
    activate_static_voxels_near_dynamic();
    prepare_tether_forces();
    wake_sleeping_voxels_near_awake();
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
    if (state->sleep_snap_original_test) {
        if (step == 1 && !state->activation_succeeded) {
            if (!debug_activate_test_cluster(state)) return;
            debug_rebuild_world_state();
            debug_place_tagged_cluster_near_original(
                DEBUG_TEST_TAG, (Vector3){ VOXEL_SIZE * 0.05f, 0.0f, 0.0f });
            for (int i = 0; i < voxel_count; ++i) {
                Voxel *voxel = &voxels[i];
                if (voxel->debugClusterTag != DEBUG_TEST_TAG || !voxel->simulate) continue;
                voxel->wake_timer = 120;
                voxel->sleepFrames = 0;
            }
            state->activation_step = step;
            state->sleep_snap_active_observed =
                debug_count_tagged(DEBUG_TEST_TAG, 1) == state->expected_voxels;
            debug_capture(state, "phase-02-activated-dynamic.png", "activated-dynamic", step);
        } else if (step > 1 && !state->sleep_snap_original_succeeded) {
            if (debug_count_tagged(DEBUG_TEST_TAG, 1) == state->expected_voxels) {
                state->sleep_snap_active_observed = true;
            }
            if (step == 20) {
                // Return the deliberately motionless test island to a small
                // residual offset, then let the production calm-island path
                // make the static-versus-in-place-sleep decision.
                debug_place_tagged_cluster_near_original(
                    DEBUG_TEST_TAG, (Vector3){ VOXEL_SIZE * 0.05f, 0.0f, 0.0f });
                for (int i = 0; i < voxel_count; ++i) {
                    Voxel *voxel = &voxels[i];
                    if (voxel->debugClusterTag != DEBUG_TEST_TAG || !voxel->simulate) continue;
                    voxel->wake_timer = 0;
                    voxel->sleepFrames = VOXEL_DEACTIVATION_FRAMES - 1;
                }
                state->sleeping_conversion_succeeded = deactivate_sleeping_voxels();
                int tagged_static = debug_count_tagged(DEBUG_TEST_TAG, 0);
                int tagged_dynamic = debug_count_tagged(DEBUG_TEST_TAG, 1);
                bool original_cells_static = true;
                for (int i = 0; i < state->tracked_cell_count; ++i) {
                    DebugTrackedCell *cell = &state->tracked_cells[i];
                    int idx = table_get_static_only(cell->x, cell->y, cell->z);
                    if (idx < 0 || idx >= voxel_count ||
                        voxels[idx].debugClusterTag != DEBUG_TEST_TAG) {
                        original_cells_static = false;
                        break;
                    }
                }
                state->sleep_snap_original_succeeded =
                    state->sleeping_conversion_succeeded &&
                    tagged_static == state->expected_voxels && tagged_dynamic == 0 &&
                    original_cells_static;
                state->sleep_snap_static_step = state->sleep_snap_original_succeeded ? step : -1;
                if (!state->sleep_snap_original_succeeded) {
                    debug_add_failure(state,
                                      "near-origin island did not snap atomically "
                                      "(static=%d/%d dynamic=%d originalCells=%s)",
                                      tagged_static, state->expected_voxels, tagged_dynamic,
                                      original_cells_static ? "true" : "false");
                }
                debug_capture(state, "phase-03-snapped-static.png", "snapped-original-static", step);
            }
        }
        return;
    }
    if (state->tether_thin_wall_ccd_test) {
        if (step == 1) {
            state->tether_thin_wall_checked = true;
            int thrown_index = find_voxel_by_identity(state->tether_thin_wall_identity);
            if (thrown_index < 0) {
                debug_add_failure(state, "tether-thrown voxel disappeared before thin-wall validation");
                return;
            }
            Voxel *thrown = &voxels[thrown_index];
            VoxelWorldBounds bounds;
            voxel_particle_world_bounds(thrown, &bounds);
            state->tether_thin_wall_final_x = thrown->pos.x;
            state->tether_thin_wall_final_min_x = bounds.minx;
            state->tether_thin_wall_final_max_x = bounds.maxx;
            state->tether_thin_wall_tunneled =
                bounds.minx > state->tether_thin_wall_wall_max_x + VOXEL_SIZE * 0.05f;
            state->tether_thin_wall_hit =
                tetherThrowCcdHitCount > state->tether_thin_wall_hits_before;
            state->tether_thin_wall_rebound_velocity_x = thrown->vel.x;
            state->tether_thin_wall_bounced = thrown->vel.x < -1.0f;

            if (tetherThrowCcdEnabled) {
                if (!state->tether_thin_wall_hit) {
                    debug_add_failure(state,
                                      "swept tether collision did not report the crossed wall");
                }
                if (bounds.maxx > state->tether_thin_wall_wall_min_x + VOXEL_SIZE * 0.05f) {
                    debug_add_failure(state,
                                      "CCD left tether voxel through wall (bounds %.3f..%.3f wall %.3f..%.3f)",
                                      bounds.minx, bounds.maxx,
                                      state->tether_thin_wall_wall_min_x,
                                      state->tether_thin_wall_wall_max_x);
                }
                if (!state->tether_thin_wall_bounced) {
                    debug_add_failure(state,
                                      "CCD stopped but did not rebound tether voxel (velocityX %.3f)",
                                      thrown->vel.x);
                }
            } else if (state->tether_thin_wall_tunneled) {
                debug_add_failure(state,
                                  "reproduced endpoint tunneling: tether voxel bounds %.3f..%.3f crossed wall %.3f..%.3f",
                                  bounds.minx, bounds.maxx,
                                  state->tether_thin_wall_wall_min_x,
                                  state->tether_thin_wall_wall_max_x);
            } else {
                debug_add_failure(state,
                                  "CCD-disabled run did not reproduce thin-wall tunneling");
            }
        }
        if (tetherThrowCcdEnabled && step >= 2 &&
            !state->tether_thin_wall_separation_checked) {
            int wall_dynamic = debug_count_tagged(DEBUG_SUPPORT_TAG, 1);
            if (wall_dynamic > 0) {
                state->tether_thin_wall_wall_woke = true;
                state->tether_thin_wall_separation_checked = true;
                int thrown_index = find_voxel_by_identity(state->tether_thin_wall_identity);
                if (thrown_index < 0) {
                    debug_add_failure(state,
                                      "tether voxel disappeared while the impacted wall woke");
                    return;
                }
                if (collisionTopologyDirty) rebuild_particle_collision_metadata();
                Voxel *thrown = &voxels[thrown_index];
                float wall_velocity_x = 0.0f;
                float wall_momentum_x = 0.0f;
                int wall_velocity_count = 0;
                for (int wall_index = 0; wall_index < voxel_count; ++wall_index) {
                    Voxel *wall = &voxels[wall_index];
                    if (!wall->simulate || wall->debugClusterTag != DEBUG_SUPPORT_TAG) continue;
                    wall_velocity_x += wall->vel.x;
                    wall_momentum_x += wall->vel.x * voxel_mass_units(wall);
                    ++wall_velocity_count;
                }
                if (wall_velocity_count > 0) {
                    state->tether_thin_wall_mean_velocity_x =
                        wall_velocity_x / (float)wall_velocity_count;
                    state->tether_thin_wall_momentum_x = wall_momentum_x;
                    state->tether_thin_wall_received_impulse = wall_momentum_x > 1.0f;
                }
                for (int thrown_corner = 0; thrown_corner < VOXEL_CORNER_COUNT;
                     ++thrown_corner) {
                    Particle *thrown_particle = thrown->particles[thrown_corner];
                    bool shared = false;
                    bool matching_group = false;
                    for (int wall_index = 0; wall_index < voxel_count; ++wall_index) {
                        Voxel *wall = &voxels[wall_index];
                        if (!wall->simulate || wall->debugClusterTag != DEBUG_SUPPORT_TAG) continue;
                        for (int wall_corner = 0; wall_corner < VOXEL_CORNER_COUNT;
                             ++wall_corner) {
                            Particle *wall_particle = wall->particles[wall_corner];
                            if (wall_particle == thrown_particle) shared = true;
                            if (thrown_particle && wall_particle &&
                                thrown_particle->collision_group >= 0 &&
                                thrown_particle->collision_group == wall_particle->collision_group) {
                                matching_group = true;
                            }
                        }
                    }
                    if (shared) ++state->tether_thin_wall_shared_particles;
                    if (matching_group) ++state->tether_thin_wall_matching_groups;
                }
                if (state->tether_thin_wall_shared_particles != 0) {
                    debug_add_failure(state,
                                      "woken wall adopted %d tether-voxel particles",
                                      state->tether_thin_wall_shared_particles);
                }
                if (state->tether_thin_wall_matching_groups != 0) {
                    debug_add_failure(state,
                                      "woken wall assigned %d tether corners to its collision group",
                                      state->tether_thin_wall_matching_groups);
                }
                if (thrown->glueEligible) {
                    debug_add_failure(state,
                                      "tether impactor became glue-eligible during wall activation");
                }
                if (!state->tether_thin_wall_received_impulse) {
                    debug_add_failure(state,
                                      "woken wall received no forward collision impulse "
                                      "(mean velocityX %.3f momentumX %.3f)",
                                      state->tether_thin_wall_mean_velocity_x,
                                      state->tether_thin_wall_momentum_x);
                }
            }
        }
        return;
    }
    if (state->tether_active_then_static_test) {
        if (step == 4) {
            release_tether(0);
            prepare_tether_forces();
            players[0].pos = (Vector3){ 0.0f, 2.0f, 0.0f };
            Vector3 second_target = { 8.5f * VOXEL_SIZE, 8.5f * VOXEL_SIZE,
                                      0.5f * VOXEL_SIZE };
            Vector3 second_aim = v_norm(v_sub(second_target, players[0].pos));
            dir_to_yaw_pitch(second_aim, &players[0].yaw, &players[0].pitch);
            state->tether_voxel_count_before_static = voxel_count;
            start_tether(0);
            state->tether_second_started = players[0].tetherHolding &&
                players[0].tetherVoxel >= 0 && players[0].tetherVoxel < voxel_count;
            if (state->tether_second_started) {
                prepare_tether_forces();
                if (!resolve_tether_voxel(&players[0])) {
                    debug_add_failure(state, "second tether lost its target during static conversion");
                    return;
                }
                Voxel *selected = &voxels[players[0].tetherVoxel];
                state->tether_static_start = selected->pos;
                selected->debugClusterTag = DEBUG_TEST_TAG;
                state->tether_static_remaining = debug_count_tagged(DEBUG_SUPPORT_TAG, 0);
                // Move the hand well beyond its initial hold point so the test
                // distinguishes a live tether from a merely surviving voxel.
                players[0].pos = (Vector3){ -2.0f, 2.0f, 0.0f };
                prepare_tether_forces();
            }
        } else if (step > 4 && step <= 24 && state->tether_second_started) {
            resolve_tether_voxel(&players[0]);
            int selected = players[0].tetherVoxel;
            int actual = debug_find_tagged_voxel(DEBUG_TEST_TAG, false);
            state->tether_expected_index = actual;
            state->tether_selected_index = selected;
            if (selected != actual && state->tether_index_mismatch_step < 0)
                state->tether_index_mismatch_step = step;
            if (players[0].tetherHolding && selected >= 0 && selected < voxel_count &&
                selected == actual && voxels[selected].simulate) {
                state->tether_second_survived = true;
                Vector3 travel = v_sub(voxels[selected].pos, state->tether_static_start);
                float distance = v_length(travel);
                if (distance > state->tether_static_pull_distance)
                    state->tether_static_pull_distance = distance;
            }
        } else if (step == 25 && players[0].tetherHolding) {
            release_tether(0);
            prepare_tether_forces();
        }
    }
    bool lifecycle_test = state->falling_pillar_test || state->platform_lifecycle_test;
    if (lifecycle_test) {
        int tagged_total = sample.tagged_dynamic + sample.tagged_static;
        if (state->lifecycle_sleep_step < 0 &&
            sample.tagged_sleeping == state->expected_voxels) {
            state->lifecycle_sleep_step = step;
            if (state->restoration_deadline_test &&
                !state->restoration_deadline_hold_applied) {
                // Keep this island sleeping beyond its ordinary recycle age.
                // Its queued-original deadline must still remove it before the
                // first static source voxel is restored.
                for (int i = 0; i < voxel_count; ++i) {
                    if (voxels[i].debugClusterTag == DEBUG_TEST_TAG &&
                        voxels[i].simulate && voxels[i].sleeping) {
                        voxels[i].lifeFrames = -RECYCLE_STATIC_RESTORE_DELAY;
                    }
                }
                state->restoration_deadline_hold_applied = true;
            }
        }
        if (state->lifecycle_recycle_step < 0 && step > 0 && tagged_total == 0) {
            state->lifecycle_recycle_step = step;
        }
        if (tagged_total > 0 && sample.tagged_static != 0) {
            debug_add_failure(state,
                              "step %d: lifecycle conversion produced %d grid-static voxels",
                              step, sample.tagged_static);
        }
        if (step == RECYCLE_DYNAMIC_MAX_FRAMES + 1 &&
            state->lifecycle_sleep_step < 0 && state->lifecycle_recycle_step < 0) {
            debug_add_failure(state,
                              "step %d: active lifecycle group never entered sleep",
                              step);
        }
        if (!state->restoration_deadline_test &&
            state->lifecycle_sleep_step >= 0 &&
            step == state->lifecycle_sleep_step + RECYCLE_DYNAMIC_MAX_FRAMES + 1) {
            int queued = recycleQueueCount - state->recycle_queue_before;
            if (tagged_total != 0 || queued != state->expected_voxels) {
                debug_add_failure(state,
                                  "step %d: sleeping lifecycle recycle left %d voxels and queued %d/%d",
                                  step, tagged_total, queued, state->expected_voxels);
            }
        }
        if (state->restoration_deadline_test &&
            step == state->restoration_deadline_frame && tagged_total != 0) {
            debug_add_failure(state,
                              "step %d: restoration deadline left %d dynamic source voxels",
                              step, tagged_total);
        }
    }
    if (state->recycle_test) {
        if (sample.tagged_dynamic != 0 || sample.tagged_static != 0 ||
            sample.tagged_sleeping != 0) {
            debug_add_failure(state, "step %d: recycled sleeping island reappeared", step);
        }
        debug_validate_gpu_static_state(state, step);
        debug_readback_gpu_buffers(state, step);
        return;
    }
    if (state->delayed_activation) {
        if (!state->activation_succeeded &&
            sample.tagged_dynamic == state->expected_voxels && sample.tagged_static == 0 &&
            sample.tagged_sleeping == 0) {
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
    } else if (step <= 1 && !state->tether_active_then_static_test &&
               !state->sleep_snap_original_test) {
        debug_validate_tagged_dynamic(state, &sample, step == 0 ? "after activation" : "first physics step");
    }
    // The snap scenario intentionally has no dynamic particles after setup, so
    // a requested GPU backend has no dispatch that could upload the new static
    // generation. Its assertions inspect the authoritative CPU static tables.
    if (!state->single_corner_hinge_test && !state->sleep_snap_original_test) {
        debug_validate_gpu_static_state(state, step);
        debug_readback_gpu_buffers(state, step);
    }
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
    if (state->pillar_impact_drift_test && step >= 240) {
        Vector3 centroid = { 0.0f, 0.0f, 0.0f };
        float horizontal_speed = 0.0f;
        if (debug_tagged_centroid_and_horizontal_speed(DEBUG_TEST_TAG, &centroid,
                                                       &horizontal_speed)) {
            if (step == 240) state->drift_window_centroid = centroid;
            float dx = centroid.x - state->drift_window_centroid.x;
            float dz = centroid.z - state->drift_window_centroid.z;
            state->drift_window_distance = sqrtf(dx * dx + dz * dz);
            state->drift_window_speed = horizontal_speed;
        }
    }
    if (state->single_corner_hinge_test) {
        int lower_idx = debug_find_tagged_voxel(DEBUG_SUPPORT_TAG, false);
        int upper_idx = debug_find_tagged_voxel(DEBUG_TEST_TAG, false);
        if (lower_idx >= 0 && upper_idx >= 0) {
            float angle = debug_voxel_relative_angle_deg(&voxels[lower_idx], &voxels[upper_idx]);
            if (angle > state->max_hinge_angle_deg) state->max_hinge_angle_deg = angle;
            state->hinge_broken = !voxels[lower_idx].glued_faces[2] ||
                                  !voxels[upper_idx].glued_faces[3];
        }
    }
}

static void debug_finalize_validation(DebugRunState *state)
{
    if (physics_backend_is_gpu(physicsBackend.requested) &&
        physicsBackend.active != physicsBackend.requested) {
        debug_add_failure(state, "forced GPU run fell back to %s: %s",
                          physics_backend_name(physicsBackend.active),
                          physicsBackend.fallback_reason);
    }
    if (state->sleep_snap_original_test) {
        if (!state->sleep_snap_active_observed) {
            debug_add_failure(state, "near-origin test never visibly entered the dynamic state");
        }
        if (!state->sleep_snap_original_succeeded) {
            debug_add_failure(state, "near-origin sleep transition did not restore static voxels");
        }
        return;
    }
    if (state->pillar_impact_drift_test) {
        if (!isfinite(state->drift_window_distance) || state->drift_window_distance > 1.00f) {
            debug_add_failure(state,
                              "settled pillar drifted %.3f metres after step 240 (maximum 1.000)",
                              state->drift_window_distance);
        }
        if (!isfinite(state->drift_window_speed) || state->drift_window_speed > 0.15f) {
            debug_add_failure(state,
                              "settled pillar retained %.3f m/s mean horizontal speed (maximum 0.150)",
                              state->drift_window_speed);
        }
        return;
    }
    if (state->tether_thin_wall_ccd_test) {
        if (!state->tether_thin_wall_checked) {
            debug_add_failure(state, "thin-wall tether collision was never validated");
        }
        if (tetherThrowCcdEnabled &&
            (!state->tether_thin_wall_hit || state->tether_thin_wall_tunneled)) {
            debug_add_failure(state, "tether throw CCD did not stop the voxel at the wall");
        }
        if (tetherThrowCcdEnabled && !state->tether_thin_wall_bounced) {
            debug_add_failure(state, "tether throw CCD did not produce a rebound");
        }
        if (tetherThrowCcdEnabled && !state->tether_thin_wall_received_impulse) {
            debug_add_failure(state, "tether throw CCD did not transfer momentum to the wall");
        }
        if (tetherThrowCcdEnabled &&
            (!state->tether_thin_wall_wall_woke ||
             !state->tether_thin_wall_separation_checked)) {
            debug_add_failure(state,
                              "thin wall did not wake for post-impact group validation");
        }
        return;
    }
    if (state->tether_active_then_static_test) {
        if (!state->tether_second_started) {
            debug_add_failure(state, "second tether did not acquire the static voxel");
        }
        if (!state->tether_second_survived || voxel_count != state->tether_voxel_count_before_static) {
            debug_add_failure(state,
                              "static tether target disappeared (survived=%d voxels=%d/%d)",
                              state->tether_second_survived ? 1 : 0, voxel_count,
                              state->tether_voxel_count_before_static);
        }
        if (state->tether_static_remaining != 2) {
            debug_add_failure(state, "static tether removed %d neighbouring voxels (remaining=%d/2)",
                              2 - state->tether_static_remaining, state->tether_static_remaining);
        }
        if (!isfinite(state->tether_static_pull_distance) ||
            state->tether_static_pull_distance < 0.25f) {
            debug_add_failure(state, "second tether pulled its target only %.3f metres",
                              state->tether_static_pull_distance);
        }
        return;
    }
    if (state->single_corner_hinge_test) {
        if (state->max_hinge_angle_deg <= GLUE_BREAK_HINGE_ANGLE_DEG) {
            debug_add_failure(state, "hinge reached only %.1f degrees; reproduction was not exercised",
                              state->max_hinge_angle_deg);
        }
        if (!state->hinge_broken) {
            debug_add_failure(state,
                              "single-corner glue survived %.1f degrees (limit %.1f)",
                              state->max_hinge_angle_deg, (double)GLUE_BREAK_HINGE_ANGLE_DEG);
        }
        return;
    }
    if (state->recycle_test) {
        if (!state->sleeping_conversion_succeeded || !state->recycle_succeeded) {
            debug_add_failure(state, "sleeping island did not complete atomic recycling");
        }
        if ((state->recycle_queue_after - state->recycle_queue_before) !=
            state->expected_voxels) {
            debug_add_failure(state, "sleeping island restoration snapshots were not preserved");
        }
        return;
    }
    if (!state->activation_succeeded) {
        debug_add_failure(state, "scenario did not complete its activation stage");
    }
    if (!state->drop_sampled) {
        debug_add_failure(state, "scenario did not reach the required 60-step drop sample");
    } else if (!state->overhang_impact && !state->platform_lifecycle_test &&
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
            debug_add_failure(state, "sleep-wake scenario never froze the cluster in place");
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
    fprintf(file, ",\n  \"passed\": %s,\n  \"steps\": %d,\n  \"fixedDt\": %.9g,\n",
            passed ? "true" : "false", debugRunConfig.steps, (double)PBD_MAX_STEP_DT);
    fprintf(file,
            "  \"animation\": {\"requested\": %s, \"format\": \"gif\","
            " \"frames\": %d, \"fpsRequested\": %d, \"stepStride\": %d, \"path\": ",
            debugRunConfig.record_gif ? "true" : "false",
            state->gif_frame_count, debugRunConfig.gif_fps, state->gif_step_stride);
    debug_write_json_string(file, state->gif_path);
    fputs("},\n", file);
    fprintf(file, "  \"fallback\": %s,\n  \"fallbackReason\": ",
            physicsBackend.sticky_fallback ? "true" : "false");
    debug_write_json_string(file, physicsBackend.fallback_reason);
    fprintf(file,
            ",\n  \"activationSucceeded\": %s,\n  \"sleepingConversionSucceeded\": %s,\n"
            "  \"sleepSnapOriginalTest\": %s,\n  \"sleepSnapActiveObserved\": %s,\n"
            "  \"sleepSnapOriginalSucceeded\": %s,\n  \"sleepSnapStaticStep\": %d,\n"
            "  \"supportRemovalSucceeded\": %s,\n  \"activationStep\": %d,"
            "\n  \"dropAtStep60\": %.9g,\n",
            state->activation_succeeded ? "true" : "false",
            state->sleeping_conversion_succeeded ? "true" : "false",
            state->sleep_snap_original_test ? "true" : "false",
            state->sleep_snap_active_observed ? "true" : "false",
            state->sleep_snap_original_succeeded ? "true" : "false",
            state->sleep_snap_static_step,
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

    fprintf(file,
            "  \"impactDiagnostics\": {\"pillarDriftEnabled\": %s,"
            " \"driftAfterStep240\": %.9g, \"meanHorizontalSpeed\": %.9g,"
            " \"singleCornerHingeEnabled\": %s, \"maxHingeAngleDeg\": %.9g,"
            " \"hingeBroken\": %s},\n",
            state->pillar_impact_drift_test ? "true" : "false",
            (double)state->drift_window_distance, (double)state->drift_window_speed,
            state->single_corner_hinge_test ? "true" : "false",
            (double)state->max_hinge_angle_deg,
            state->hinge_broken ? "true" : "false");

    fprintf(file,
            "  \"tetherSequenceDiagnostics\": {\"enabled\": %s,"
            " \"secondStarted\": %s, \"secondSurvived\": %s,"
            " \"staticRemaining\": %d, \"pullDistance\": %.9g,"
            " \"indexMismatchStep\": %d, \"expectedIndex\": %d,"
            " \"selectedIndex\": %d},\n",
            state->tether_active_then_static_test ? "true" : "false",
            state->tether_second_started ? "true" : "false",
            state->tether_second_survived ? "true" : "false",
            state->tether_static_remaining,
            (double)state->tether_static_pull_distance,
            state->tether_index_mismatch_step, state->tether_expected_index,
            state->tether_selected_index);

    fprintf(file,
            "  \"tetherThinWallDiagnostics\": {\"enabled\": %s,"
            " \"ccdOptionEnabled\": %s, \"checked\": %s, \"hit\": %s,"
            " \"tunneled\": %s, \"bounced\": %s, \"wallWoke\": %s,"
            " \"separationChecked\": %s, \"wallReceivedImpulse\": %s,"
            " \"sharedParticles\": %d,"
            " \"matchingCollisionGroups\": %d, \"reboundVelocityX\": %.9g,"
            " \"wallMeanVelocityX\": %.9g, \"wallMomentumX\": %.9g,"
            " \"initialX\": %.9g, \"finalX\": %.9g,"
            " \"finalMinX\": %.9g, \"finalMaxX\": %.9g,"
            " \"wallMinX\": %.9g, \"wallMaxX\": %.9g},\n",
            state->tether_thin_wall_ccd_test ? "true" : "false",
            tetherThrowCcdEnabled ? "true" : "false",
            state->tether_thin_wall_checked ? "true" : "false",
            state->tether_thin_wall_hit ? "true" : "false",
            state->tether_thin_wall_tunneled ? "true" : "false",
            state->tether_thin_wall_bounced ? "true" : "false",
            state->tether_thin_wall_wall_woke ? "true" : "false",
            state->tether_thin_wall_separation_checked ? "true" : "false",
            state->tether_thin_wall_received_impulse ? "true" : "false",
            state->tether_thin_wall_shared_particles,
            state->tether_thin_wall_matching_groups,
            (double)state->tether_thin_wall_rebound_velocity_x,
            (double)state->tether_thin_wall_mean_velocity_x,
            (double)state->tether_thin_wall_momentum_x,
            (double)state->tether_thin_wall_initial_x,
            (double)state->tether_thin_wall_final_x,
            (double)state->tether_thin_wall_final_min_x,
            (double)state->tether_thin_wall_final_max_x,
            (double)state->tether_thin_wall_wall_min_x,
            (double)state->tether_thin_wall_wall_max_x);

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
                " \"taggedStatic\": %d, \"taggedSleeping\": %d, \"simParticles\": %d,"
                " \"collisionParticles\": %d, \"collisionParticlesCulled\": %d,"
                " \"staticOccupiedCells\": %d, \"staticSurfacePatches\": %d,"
                " \"staticSurfaceCells\": %d, \"staticSurfacePatchRefs\": %d,"
                " \"zeroMassCorners\": %d, \"invalidParticleRefs\": %d,"
                " \"staticGluedFaces\": %d, \"centroidY\": %.9g,"
                " \"centroidDrop\": %.9g, \"minY\": %.9g, \"maxY\": %.9g,"
                " \"meanVelocityY\": %.9g, \"staticHashGeneration\": %llu,"
                " \"gpuStaticUploadedGeneration\": %llu, \"physicsStepMs\": %.9g}%s\n",
                sample->step, sample->tagged_voxels, sample->tagged_dynamic,
                sample->tagged_static, sample->tagged_sleeping, sample->sim_particles,
                sample->collision_particles, sample->collision_particles_culled,
                sample->static_occupied_cells, sample->static_surface_patches,
                sample->static_surface_cells, sample->static_surface_patch_refs,
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
        debug_gif_begin(&state);
        debug_validate_step(&state, 0);
        debug_gif_capture(&state, 0);
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
            debug_gif_capture(&state, step);
            if (debug_should_capture_step(step)) {
                char file_name[64];
                snprintf(file_name, sizeof(file_name), "step-%04d.png", step);
                debug_capture(&state, file_name, "simulation", step);
            }
        }
        debug_finalize_validation(&state);
    }

    debug_gif_end(&state);
    bool report_ok = debug_write_report(&state);
    bool passed = setup_ok && report_ok && state.failure_count == 0;
    fprintf(stderr,
            "debug-run scenario=%s requested=%s active=%s result=%s steps=%d drop60=%.6f captures=%d gif=%s frames=%d report=%s\n",
            debugRunConfig.scenario,
            physics_backend_name(physicsBackend.requested),
            physics_backend_name(physicsBackend.active),
            passed ? "PASS" : "FAIL",
            debugRunConfig.steps, state.drop_at_step_60,
            state.capture_count,
            debugRunConfig.record_gif ? state.gif_path : "disabled",
            state.gif_frame_count, state.report_path);
    for (int i = 0; i < state.failure_count; ++i) {
        fprintf(stderr, "debug-failure: %s\n", state.failures[i]);
    }
    UnloadRenderTexture(state.target);
    return passed ? 0 : 4;
}
