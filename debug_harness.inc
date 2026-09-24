#include <errno.h>
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

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
#define DEBUG_CONTAINER_TAG 60
#define DEBUG_FLUID_TAG 59
#define DEBUG_EXPECTED_VOXELS 18
#define DEBUG_MIN_STEPS 60
#define DEBUG_DEFAULT_STEPS 60
#define DEBUG_RENDER_WIDTH 1280
#define DEBUG_RENDER_HEIGHT 720
#define DEBUG_GIF_WIDTH 640
#define DEBUG_GIF_HEIGHT 360
#define DEBUG_GIF_HASH_SIZE 8191
#define DEBUG_REFINE_AT_STEP 120

typedef struct DebugRunConfig {
    bool requested;
    bool show_window;
    bool record_gif;
    bool capture_steps_explicit;
    int steps;
    int gif_fps;
    int capture_steps[DEBUG_MAX_CAPTURE_STEPS];
    int capture_step_count;
    int pbf_side;
    char scenario[64];
    char output_dir[DEBUG_PATH_CAPACITY];
} DebugRunConfig;

typedef struct DebugSample {
    char phase[48];
    int step;
    int tagged_voxels;
    int tagged_dynamic;
    int tagged_rendered;
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
    int fluid_particles;
    int invalid_fluid_particles;
    float fluid_centroid_y;
    float mean_density_ratio;
    float max_density_ratio;
    float mean_fluid_speed;
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
    int expected_static_leftovers;
    float effective_pbd_size;
    bool active_cube_drop;
    bool adaptive_cube_drop;
    bool rip_test;
    FILE *rip_history_file;
    char rip_history_path[DEBUG_PATH_CAPACITY];
    FILE *rip_strain_levels_file;
    char rip_strain_levels_path[DEBUG_PATH_CAPACITY];
    int rip_substeps;
    int rip_first_break_step;
    float rip_left_x;
    float rip_right_x;
    float rip_grip_depth;
    float rip_commanded_strain;
    int rip_broken_nodes;
    int rip_broken_leaves;
    FILE *curvature_levels_file;
    char curvature_levels_path[DEBUG_PATH_CAPACITY];
    bool refinement_demo;
    bool refinement_complete;
    bool coarsening_demo;
    bool coarsening_complete;
    float coarsening_projection_rms;
    int refinement_step;
    int refinement_elapsed;
    Particle *refine_grid[27];
    float refine_inv_mass[27];
    int refine_children[8];
    bool activation_succeeded;
    bool sleeping_conversion_succeeded;
    bool sleep_snap_original_test;
    bool sleep_snap_rebuild_test;
    bool sleep_snap_original_succeeded;
    bool sleep_snap_active_observed;
    bool sleep_snap_sleeping_observed;
    bool sleep_snap_impact_activation_requested;
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
    int sleep_snap_queue_at_static;
    uint64_t sleep_snap_static_identity_signature;
    bool sleep_snap_rebuilt_after_static;
    float sleep_snap_max_original_error;
    float sleep_snap_last_speed;
    float sleep_snap_last_strain;
    float sleep_snap_last_shear;
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
    bool voxel_fracture_test;
    int fracture_strain_voxel;
    int fracture_shear_voxel;
    int fracture_control_voxel;
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
    bool pbf_container;
    bool pbf_dynamic_container;
    int expected_fluid_particles;
    int expected_container_voxels;
    float initial_fluid_centroid_y;
    float fluid_drop_at_step_120;
    float container_min_x;
    float container_max_x;
    float container_min_z;
    float container_max_z;
    float container_floor_top;
    int escaped_fluid_particles;
    int below_container_particles;
    bool pbf_scaling_benchmark;
    int benchmark_warmup_steps;
    int benchmark_measured_steps;
    double benchmark_physics_total_ms;
    double benchmark_render_total_ms;
    double benchmark_frame_total_ms;
    double benchmark_gpu_compute_total_ms;
    double benchmark_gpu_hull_total_ms;
    double benchmark_gpu_pack_total_ms;
    double benchmark_gpu_upload_total_ms;
    double benchmark_gpu_readback_total_ms;
    uint64_t benchmark_gpu_upload_bytes;
    uint64_t benchmark_gpu_readback_bytes;
    int benchmark_gpu_dispatches;
    int benchmark_gpu_hash_capacity;
    uint32_t benchmark_gpu_neighbor_overflow;
    bool benchmark_gpu_hull;
    bool expect_gpu_neighbor_overflow;
    bool arena_platform_fluid_test;
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
    bool swarmer_wall_test;
    bool goliath_flight_test;
    bool goliath_behavior_test;
    int enemy_bot_index;
    int enemy_human_index;
    int swarmer_inside_wall_frames;
    float swarmer_min_z;
    float goliath_start_y;
    float goliath_max_y;
    float goliath_y_after_landing;
    int goliath_hover_frames;
    int goliath_forced_jumps;
    int goliath_visible_frames;
    int goliath_voxels_launched;
    int goliath_min_armor_count;
    int goliath_final_armor_count;
} DebugRunState;

typedef bool (*DebugScenarioSetupFn)(DebugRunState *state);

typedef struct DebugScenarioDefinition {
    const char *name;
    const char *description;
    DebugScenarioSetupFn setup;
    int default_steps;
    int minimum_steps;
} DebugScenarioDefinition;

static DebugRunConfig debugRunConfig = {
    .steps = DEBUG_DEFAULT_STEPS,
    .gif_fps = 20,
    .pbf_side = 4
};

static bool debug_setup_dynamic_freefall(DebugRunState *state);
static bool debug_setup_creative_reset_cache(DebugRunState *state);
static bool debug_setup_activation_floating(DebugRunState *state);
static bool debug_setup_sleep_wake_floating(DebugRunState *state);
static bool debug_setup_sleep_snap_original(DebugRunState *state);
static bool debug_setup_sleep_snap_no_rebuild(DebugRunState *state);
static bool debug_setup_sleep_recycle(DebugRunState *state);
static bool debug_setup_falling_pillar_sleep(DebugRunState *state);
static bool debug_setup_three_leg_cut_platform(DebugRunState *state);
static bool debug_setup_tether_throw_floating(DebugRunState *state);
static bool debug_setup_overhang_impact(DebugRunState *state);
static bool debug_setup_pillar_impact_drift(DebugRunState *state);
static bool debug_setup_single_corner_hinge(DebugRunState *state);
static bool debug_setup_tether_active_then_static(DebugRunState *state);
static bool debug_setup_tether_thin_wall_ccd(DebugRunState *state);
static bool debug_setup_pbf_container(DebugRunState *state);
static bool debug_setup_tall_pillar_snap(DebugRunState *state);
static bool debug_setup_swarmer_wall(DebugRunState *state);
static bool debug_setup_goliath_flight(DebugRunState *state);
static bool debug_setup_goliath_behavior(DebugRunState *state);
static bool debug_setup_pbf_dynamic_container(DebugRunState *state);
static bool debug_setup_pbf_scaling(DebugRunState *state);
static bool debug_setup_pbf_dense_neighbors(DebugRunState *state);
static bool debug_setup_cube_stress(DebugRunState *state);
static bool debug_setup_active_cube_drop(DebugRunState *state);
static bool debug_setup_adaptive_cube_drop(DebugRunState *state);
static bool debug_setup_rip_test(DebugRunState *state);
static bool debug_setup_vgs_refinement(DebugRunState *state);
static bool debug_setup_vgs_coarsening(DebugRunState *state);
static bool debug_setup_voxel_strain_shear_fracture(DebugRunState *state);
static bool debug_setup_arena_platform_fluid(DebugRunState *state);

static const DebugScenarioDefinition debugScenarios[] = {
    { "dynamic-freefall", "Dynamic cluster freefall control", debug_setup_dynamic_freefall, 60, 60 },
    { "creative-reset-cache", "Creative reset invalidates cached static render and collision surfaces", debug_setup_creative_reset_cache, 60, 60 },
    { "activation-floating", "Static cluster activated while unsupported", debug_setup_activation_floating, 60, 60 },
    { "sleep-wake-floating", "Supported cluster sleeps, loses support, then activates", debug_setup_sleep_wake_floating, 60, 60 },
    { "sleep-snap-original", "A calm 4x4x6 pillar near its original pose snaps back to static", debug_setup_sleep_snap_original, 60, 60 },
    { "tall-pillar-snap", "Top half of a tall pillar settles and restores to static", debug_setup_tall_pillar_snap, 120, 60 },
    { "sleep-snap-no-rebuild", "A restored pillar remains static after its queued rebuild deadline", debug_setup_sleep_snap_no_rebuild, 60, 60 },
    { "sleep-recycle", "Sleeping dynamic island enters delayed restoration atomically", debug_setup_sleep_recycle, 60, 60 },
    { "falling-pillar-sleep", "Rotated pillar falls, settles, and sleeps without changing pose", debug_setup_falling_pillar_sleep, 60, 60 },
    { "three-leg-cut-platform", "Four-legged platform with three legs removed completes sleep and recycling", debug_setup_three_leg_cut_platform, 60, 60 },
    { "tether-throw-floating", "Released tether voxel wakes an unsupported sleeping cluster", debug_setup_tether_throw_floating, 60, 60 },
    { "overhang-impact", "Falling dynamic voxel activates and bends a supported static overhang", debug_setup_overhang_impact, 60, 60 },
    { "pillar-impact-drift", "Large undercut pillar impact exposes post-contact sliding and wiggle", debug_setup_pillar_impact_drift, 60, 60 },
    { "single-corner-hinge", "Check one-corner glue behavior beyond the hinge-angle limit", debug_setup_single_corner_hinge, 60, 60 },
    { "tether-active-then-static", "Tethering an active voxel then a static voxel preserves and pulls both", debug_setup_tether_active_then_static, 60, 60 },
    { "tether-thin-wall-ccd", "High-speed tether throw uses swept particle rays against a one-voxel wall", debug_setup_tether_thin_wall_ccd, 60, 60 },
    { "pbf-container", "High-viscosity PBF falls into an immutable container", debug_setup_pbf_container, 240, 120 },
    { "pbf-dynamic-container", "PBF falls onto an elevated active PBD floor", debug_setup_pbf_dynamic_container, 240, 120 },
    { "pbf-scaling", "Timed PBF physics plus convex-hull rendering", debug_setup_pbf_scaling, 120, 60 },
    { "pbf-dense-neighbors", "Dense PBF neighborhood overflow fallback", debug_setup_pbf_dense_neighbors, 60, 60 },
    { "swarmer-wall", "Swarmer charges a visible target through a blocking wall", debug_setup_swarmer_wall, 90, 60 },
    { "goliath-flight", "Goliath jumps with attached armor and must land", debug_setup_goliath_flight, 120, 90 },
    { "goliath-behavior", "Goliath maintains armor, detects enemy, and launches voxels", debug_setup_goliath_behavior, 180, 60 },
    { "cube-stress", "Stress test with scalable dynamic cube (size set by FPS_CUBE_SIZE)", debug_setup_cube_stress, 60, 60 },
    { "active-cube-drop", "Activate a static 10x10x10 cube and drop it", debug_setup_active_cube_drop, 120, 60 },
    { "adaptive-cube-drop", "Curvature-driven octree cube drop", debug_setup_adaptive_cube_drop, 240, 60 },
    { "rip-test", "Pull both faces of an adaptive cube to test fracture", debug_setup_rip_test, 480, 60 },
    { "vgs-refinement", "Drop a root VGS cube, then refine into eight shared-particle children", debug_setup_vgs_refinement, 300, 180 },
    { "vgs-coarsening", "Drop eight child VGS cubes, then project them onto one root", debug_setup_vgs_coarsening, 300, 180 },
    { "voxel-strain-shear-fracture", "Check VGS deactivation behavior under high strain and shear", debug_setup_voxel_strain_shear_fracture, 60, 60 },
    { "arena-platform-fluid", "Test arena central platform activated with sleeping fluid on top", debug_setup_arena_platform_fluid, 60, 60 }
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
                               "--debug-pbf-side=", "--debug-pbf-side");
    if (value != (const char *)(uintptr_t)1) {
        return debug_parse_integer(value, 1, 25,
                                   &debugRunConfig.pbf_side, "--debug-pbf-side") ? 1 : -1;
    }

    value = debug_option_value(argc, argv, index, argument,
                               "--debug-output=", "--debug-output");
    if (value != (const char *)(uintptr_t)1) {
        return debug_copy_option(debugRunConfig.output_dir, sizeof(debugRunConfig.output_dir),
                                 value, "--debug-output") ? 1 : -1;
    }

    value = debug_option_value(argc, argv, index, argument,
                               "--debug-fluid-render=", "--debug-fluid-render");
    if (value != (const char *)(uintptr_t)1) {
        if (strcmp(value, "surface") == 0) {
            fluidRenderMode = FLUID_RENDER_SURFACE;
        } else if (strcmp(value, "particles") == 0) {
            fluidRenderMode = FLUID_RENDER_PARTICLES;
        } else {
            fprintf(stderr, "Invalid --debug-fluid-render value '%s' (expected 'surface' or 'particles')\n", value);
            return -1;
        }
        return 1;
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

static bool debug_make_directory_recursive(const char *path)
{
    if (!path || !path[0] || strlen(path) >= DEBUG_PATH_CAPACITY) return false;
    char copy[DEBUG_PATH_CAPACITY];
    strcpy(copy, path);
    for (char *cursor = copy + 1; *cursor; ++cursor) {
        if (*cursor != '/' && *cursor != '\\') continue;
#ifdef _WIN32
        if (cursor == copy + 2 && copy[1] == ':') continue;
#endif
        char separator = *cursor;
        *cursor = '\0';
#ifdef _WIN32
        bool created = _mkdir(copy) == 0 || errno == EEXIST;
#else
        bool created = mkdir(copy, 0755) == 0 || errno == EEXIST;
#endif
        *cursor = separator;
        if (!created) return false;
    }
#ifdef _WIN32
    return (_mkdir(copy) == 0 || errno == EEXIST) && DirectoryExists(path);
#else
    return (mkdir(copy, 0755) == 0 || errno == EEXIST) && DirectoryExists(path);
#endif
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

static bool debug_tilt_tagged_cube_for_corner_drop(int tag)
{
    Vector3 pivot = {0};
    int count = 0;
    int stamp = ++particle_sync_stamp;
    if (stamp == 0) stamp = particle_sync_stamp = 1;
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        if (!voxel->simulate || voxel->debugClusterTag != tag) continue;
        for (int c = 0; c < 8; ++c) {
            Particle *particle = voxel->particles[c];
            if (!particle || particle->sync_stamp == stamp) continue;
            particle->sync_stamp = stamp;
            pivot = v_add(pivot, particle->pos);
            ++count;
        }
    }
    if (!count) return false;
    pivot = v_mul(pivot, 1.0f / (float)count);

    const float angle = 45.0f * DEG2RAD;
    const float cs = cosf(angle), sn = sinf(angle);
    stamp = ++particle_sync_stamp;
    if (stamp == 0) stamp = particle_sync_stamp = 1;
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        if (!voxel->simulate || voxel->debugClusterTag != tag) continue;
        Vector3 center = {0};
        for (int c = 0; c < 8; ++c) {
            Particle *particle = voxel->particles[c];
            if (particle->sync_stamp != stamp) {
                particle->sync_stamp = stamp;
                Vector3 local = v_sub(particle->pos, pivot);
                // Apply 45 degrees about X, then 45 degrees about Z.
                float y = cs * local.y - sn * local.z;
                Vector3 rotated = {cs * local.x - sn * y,
                                   sn * local.x + cs * y,
                                   sn * local.y + cs * local.z};
                particle->pos = v_add(pivot, rotated);
                particle->prev_pos = particle->predicted_pos = particle->pos;
                particle->vel = (Vector3){0};
            }
            center = v_add(center, particle->pos);
        }
        voxel->pos = v_mul(center, 0.125f);
        voxel->vel = (Vector3){0};
    }
    rebuild_voxel_hash();
    collisionTopologyDirty = true;
    meshDirty = true;
    return true;
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
    debugColorVoxelsByVgsTemporalCurvature = false;
    debugRefinementCoarse = false;
    debugAdaptiveOctree = false;
    debugRipTest = false;
    debug_reset_vgs_temporal_curvature();
    pbdSolidVoxelSize = VOXEL_SIZE;
    clear_world_voxels();
    clear_pickups();
    voxelHashFramesSinceRebuild = 0;
    pbdTimeAccumulator = 0.0f;
    staticHashDirty = false;
    init_static_hash();
    table_cache_invalidate();
}

static void debug_glue_pbd_faces(int new_index, int neighbor_index, int face)
{
    Voxel *a = &voxels[new_index];
    Voxel *b = &voxels[neighbor_index];
    int opposite = opposite_face[face];
    for (int corner = 0; corner < 4; ++corner) {
        int ai = face_corner_indices[face][corner];
        int bi = face_corner_indices[opposite][corner];
        Particle *old = a->particles[ai];
        Particle *shared = b->particles[bi];
        if (old == shared) continue;
        if (old->glue_count > 0) --old->glue_count;
        particle_release(old);
        particle_retain(shared);
        ++shared->glue_count;
        a->particles[ai] = shared;
    }
    a->glued_faces[face] = b->glued_faces[opposite] = true;
    a->glued_neighbors[face] = neighbor_index;
    b->glued_neighbors[opposite] = new_index;
    collisionTopologyDirty = true;
}

static bool debug_remove_tagged_voxels(int tag);

static bool debug_build_cube_vgs_hierarchy(const int *lattice, int side, float leaf_edge)
{
    clear_vgs_hierarchy();
    if (side == 1) {
        vgsHierarchy.leaf_count = 1;
        vgsHierarchy.leaf_history = calloc((size_t)voxel_count, sizeof(VgsCurvatureHistory));
        return vgsHierarchy.leaf_history != NULL;
    }
    int levels = 0, node_count = 0;
    for (int width = 1; width < side; width *= 2) {
        node_count += width * width * width;
        ++levels;
    }
    VgsHierarchyNode *nodes = calloc((size_t)node_count, sizeof(*nodes));
    if (!nodes) return false;
    vgsHierarchy.nodes = nodes;
    vgsHierarchy.levels = levels;
    vgsHierarchy.node_count = node_count;
    vgsHierarchy.leaf_count = side * side * side;
    vgsHierarchy.node_history = calloc((size_t)node_count, sizeof(VgsCurvatureHistory));
    vgsHierarchy.leaf_history = calloc((size_t)voxel_count, sizeof(VgsCurvatureHistory));
    if (!vgsHierarchy.node_history || !vgsHierarchy.leaf_history) {
        clear_vgs_hierarchy();
        return false;
    }
    int offset = 0;
    for (int level = 0, grid = 1; level < levels; ++level, grid *= 2) {
        vgsHierarchy.level_start[level] = offset;
        vgsHierarchy.level_count[level] = grid * grid * grid;
        offset += grid * grid * grid;
    }
    for (int level = 0, grid = 1; level < levels; ++level, grid *= 2) {
        int span = side / grid;
        int child_grid = grid * 2;
        for (int y = 0; y < grid; ++y) for (int z = 0; z < grid; ++z) for (int x = 0; x < grid; ++x) {
            int id = vgsHierarchy.level_start[level] + (y * grid + z) * grid + x;
            VgsHierarchyNode *node = &nodes[id];
            node->parent = level ? vgsHierarchy.level_start[level - 1] +
                (((y / 2) * (grid / 2) + z / 2) * (grid / 2) + x / 2) : -1;
            node->rest_edge = span * leaf_edge;
            node->rest_volume = node->rest_edge * node->rest_edge * node->rest_edge;
            node->active = true;
            node->vgs_enabled = true;
            for (int c = 0; c < 8; ++c) {
                int cx = x * span + (corner_signs[c][0] > 0 ? span - 1 : 0);
                int cy = y * span + (corner_signs[c][1] > 0 ? span - 1 : 0);
                int cz = z * span + (corner_signs[c][2] > 0 ? span - 1 : 0);
                node->particles[c] = voxels[lattice[(cy * side + cz) * side + cx]].particles[c];
                int child_x = x * 2 + (c & 1);
                int child_y = y * 2 + ((c >> 1) & 1);
                int child_z = z * 2 + ((c >> 2) & 1);
                int child_flat = (child_y * child_grid + child_z) * child_grid + child_x;
                node->children[c] = level + 1 < levels
                    ? vgsHierarchy.level_start[level + 1] + child_flat
                    : -(lattice[child_flat] + 1);
            }
        }
    }
    return true;
}

static bool debug_validate_cube_vgs_hierarchy(int side)
{
    int expected = (side * side * side - 1) / 7;
    if (vgsHierarchy.node_count != expected || vgsHierarchy.leaf_count != side * side * side)
        return false;
    if (side > 1 && fabsf(vgsHierarchy.nodes[0].rest_edge - 5.0f) > 1e-5f)
        return false;
    for (int i = 0; i < vgsHierarchy.node_count; ++i) {
        VgsHierarchyNode *node = &vgsHierarchy.nodes[i];
        if (!node->active || node->rest_edge <= 0.0f ||
            fabsf(node->rest_volume - node->rest_edge * node->rest_edge * node->rest_edge) > 1e-4f)
            return false;
        for (int c = 0; c < 8; ++c) {
            Particle *p = node->particles[c];
            if (!p || (uintptr_t)p < (uintptr_t)particles_pool ||
                (uintptr_t)p >= (uintptr_t)(particles_pool + particle_pool_count) || !p->active)
                return false;
            int child = node->children[c];
            if (child >= 0) {
                if (child >= vgsHierarchy.node_count || vgsHierarchy.nodes[child].parent != i)
                    return false;
            } else if (-child - 1 >= voxel_count || !voxels[-child - 1].simulate) {
                return false;
            }
        }
    }
    // A broken child leaves the parent's VGS state independent.
    if (vgsHierarchy.node_count) {
        Voxel *leaf = &voxels[-vgsHierarchy.nodes[vgsHierarchy.node_count - 1].children[0] - 1];
        leaf->vgs_active = false;
        bool independent = vgsHierarchy.nodes[0].active &&
                           vgsHierarchy.nodes[0].vgs_enabled;
        leaf->vgs_active = true;
        if (!independent) return false;
    }
    return true;
}

static bool debug_convert_cube_to_pbd(DebugRunState *state, float requested_size)
{
    state->expected_voxels = 0;
    state->expected_static_leftovers = 1000;
    const double size = (double)requested_size;
    const double cube_edge = 10.0 * (double)VOXEL_SIZE;
    if (!isfinite(size) || size <= 0.0 || size > cube_edge) {
        debug_add_failure(state, "PBD size %.9g must be positive and no larger than the %.9g m cube",
                          size, cube_edge);
        return false;
    }
    // Keep the original cube bounds and use a uniform lattice. The requested
    // size is an upper bound when it does not divide the 5 m cube exactly.
    double minimum_side = ceil(cube_edge / size - 1e-6);
    int side = 1;
    while ((double)side < minimum_side && side < 64) side *= 2;
    if ((double)side < minimum_side || side > 32) {
        double lower_side = fmax((double)side, minimum_side);
        double required = lower_side * lower_side * lower_side;
        debug_add_failure(state, "PBD size %.9g requires at least %.0f voxels and %.0f particles (caps %d/%d)",
                          size, required, required * 8.0, MAX_VOXELS, MAX_PARTICLES);
        return false;
    }
    float effective_size = (float)(cube_edge / (double)side);
    int64_t dynamic_count = (int64_t)side * side * side;
    int64_t needed_particles = 8 * dynamic_count;
    if (dynamic_count > MAX_VOXELS || needed_particles > MAX_PARTICLES ||
        dynamic_count > INT_MAX) {
        debug_add_failure(state, "PBD size %.9g requires %lld voxels and %lld particles (caps %d/%d)",
                          size, (long long)dynamic_count,
                          (long long)needed_particles, MAX_VOXELS, MAX_PARTICLES);
        return false;
    }

    int *lattice = malloc((size_t)dynamic_count * sizeof(*lattice));
    float *particle_mass = calloc((size_t)MAX_PARTICLES, sizeof(*particle_mass));
    if (!lattice || !particle_mass) {
        free(lattice); free(particle_mass);
        debug_add_failure(state, "PBD size %.9g: unable to allocate conversion workspace", size);
        return false;
    }
    debug_remove_tagged_voxels(DEBUG_TEST_TAG);
    const Color color = { 245, 145, 45, 255 };
    pbdSolidVoxelSize = effective_size;
    for (int y = 0; y < side; ++y) for (int z = 0; z < side; ++z) for (int x = 0; x < side; ++x) {
        int flat = (y * side + z) * side + x;
        float px = -5.0f * VOXEL_SIZE + ((float)x + 0.5f) * pbdSolidVoxelSize;
        float py = 18.0f * VOXEL_SIZE + ((float)y + 0.5f) * pbdSolidVoxelSize +
                   (state->refinement_demo ? 15.0f : 0.0f);
        float pz = -5.0f * VOXEL_SIZE + ((float)z + 0.5f) * pbdSolidVoxelSize;
        int idx = addVoxel(px, py, pz, false, true, color, 0);
        if (idx < 0) goto conversion_failed;
        voxels[idx].debugClusterTag = DEBUG_TEST_TAG;
        lattice[flat] = idx;
        if (x > 0) debug_glue_pbd_faces(idx, lattice[flat - 1], 1);
        if (z > 0) debug_glue_pbd_faces(idx, lattice[flat - side], 5);
        if (y > 0) debug_glue_pbd_faces(idx, lattice[flat - side * side], 3);
    }
    // Each voxel contributes eight corner masses. Shared particles receive the
    // sum of contributions, preserving mass per physical volume after gluing.
    float contribution = (pbdSolidVoxelSize / VOXEL_SIZE);
    contribution = contribution * contribution * contribution;
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *v = &voxels[i];
        if (!v->simulate || v->debugClusterTag != DEBUG_TEST_TAG) continue;
        for (int c = 0; c < 8; ++c) {
            ptrdiff_t pool = v->particles[c] - particles_pool;
            particle_mass[pool] += contribution;
        }
    }
    for (int i = 0; i < particle_pool_count; ++i) {
        if (particle_mass[i] <= 0.0f) continue;
        particles_pool[i].inv_mass = particles_pool[i].base_inv_mass = 1.0f / particle_mass[i];
    }
    if (!debug_build_cube_vgs_hierarchy(lattice, side, effective_size)) {
        debug_add_failure(state, "PBD size %.9g: unable to allocate VGS hierarchy", size);
        goto conversion_failed;
    }
    if (!debug_validate_cube_vgs_hierarchy(side)) {
        debug_add_failure(state, "PBD size %.9g: invalid VGS octree", size);
        goto conversion_failed;
    }
    if (state->refinement_demo || state->coarsening_demo) {
        if (side != 2 || vgsHierarchy.node_count != 1) {
            debug_add_failure(state, "scale transition demo requires a 2x2x2 child lattice");
            goto conversion_failed;
        }
        for (int c = 0; c < 8; ++c) {
            int cx = c & 1, cy = (c >> 1) & 1, cz = (c >> 2) & 1;
            int idx = lattice[(cy * 2 + cz) * 2 + cx];
            state->refine_children[c] = idx;
            for (int k = 0; k < 8; ++k) {
                int gx = cx + (k & 1);
                int gy = cy + ((k >> 1) & 1);
                int gz = cz + ((k >> 2) & 1);
                int r = gx + 3 * gy + 9 * gz;
                Particle *particle = voxels[idx].particles[k];
                if (state->refine_grid[r] && state->refine_grid[r] != particle) {
                    debug_add_failure(state, "refinement grid point %d is not shared", r);
                    goto conversion_failed;
                }
                state->refine_grid[r] = particle;
                state->refine_inv_mass[r] = particle->base_inv_mass;
            }
        }
        for (int r = 0; r < 27; ++r) {
            if (!state->refine_grid[r]) {
                debug_add_failure(state, "refinement grid point %d is missing", r);
                goto conversion_failed;
            }
        }
    }
    double total_mass = 0.0;
    for (int i = 0; i < particle_pool_count; ++i) {
        if (particle_mass[i] > 0.0f) total_mass += 1.0 / particles_pool[i].base_inv_mass;
    }
    double expected_mass = 8.0 * 10.0 * 10.0 * 10.0;
    if (fabs(total_mass - expected_mass) > expected_mass * 1e-4) {
        debug_add_failure(state, "PBD size %.9g mass %.6g differs from expected %.6g",
                          size, total_mass, expected_mass);
    }
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *v = &voxels[i];
        if (!v->simulate || v->debugClusterTag != DEBUG_TEST_TAG) continue;
        if (fabsf(v->rest_edge - pbdSolidVoxelSize) > 1e-5f) {
            debug_add_failure(state, "PBD voxel %d has the wrong rest edge", i);
            break;
        }
        for (int c = 0; c < 8; ++c) {
            Particle *p = v->particles[c];
            Vector3 target = {
                v->pos.x + corner_signs[c][0] * 0.5f * pbdSolidVoxelSize,
                v->pos.y + corner_signs[c][1] * 0.5f * pbdSolidVoxelSize,
                v->pos.z + corner_signs[c][2] * 0.5f * pbdSolidVoxelSize
            };
            if (!p || v_length(v_sub(p->pos, target)) > 1e-4f) {
                debug_add_failure(state, "PBD voxel %d corner %d does not match its lattice position", i, c);
                break;
            }
        }
    }
    state->expected_voxels = (int)dynamic_count;
    state->expected_static_leftovers = 0;
    state->effective_pbd_size = effective_size;
    state->activation_succeeded = true;
    // The gameplay glue rebuild uses static-grid coordinates. Preserve the
    // explicit PBD-lattice links created above when refreshing debug caches.
    rebuild_voxel_hash();
    rebuild_static_hash_if_dirty();
    rebuild_all_voxel_surfaces();
    refresh_static_voxel_beliefs();
    meshDirty = true;
    free(lattice); free(particle_mass);
    return true;

conversion_failed:
    clear_vgs_hierarchy();
    debug_add_failure(state, "PBD size %.9g: cube conversion failed after preflight", size);
    free(lattice); free(particle_mass);
    return false;
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
    bool cube_dynamic_only = strcmp(debugRunConfig.scenario, "active-cube-drop") == 0 ||
                             strcmp(debugRunConfig.scenario, "adaptive-cube-drop") == 0;
    for (int i = 0; i < voxel_count; ++i) {
        const Voxel *voxel = &voxels[i];
        if (voxel->debugClusterTag != DEBUG_TEST_TAG) continue;
        ++sample.tagged_voxels;
        if (voxel->simulate) {
            ++sample.tagged_dynamic;
            if (voxel->vgs_active && (!vgsHierarchy.adaptive || voxel->simulate_dofs))
                ++sample.tagged_rendered;
        }
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
            if (voxel->simulate && particle->base_inv_mass <= 0.0f &&
                !vgsHierarchy.adaptive) {
                ++sample.zero_mass_corners;
            }
            if (voxel->simulate && (!vgsHierarchy.adaptive || particle->base_inv_mass > 0.0f) &&
                (particle->sim_index < 0 || particle->sim_index >= sim_particle_count ||
                 sim_particles[particle->sim_index] != particle)) {
                ++sample.invalid_particle_refs;
            }
        }
        if (valid_corners > 0) center = v_mul(center, 1.0f / (float)valid_corners);
        if (!cube_dynamic_only || voxel->simulate) {
            centroid_sum += center.y;
            velocity_sum += voxel->vel.y;
        }
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
    if (vgsHierarchy.adaptive)
        for (int i = 0; i < vgsHierarchy.node_count; ++i)
            if (vgsHierarchy.nodes[i].active && vgsHierarchy.nodes[i].vgs_enabled &&
                !vgs_node_has_active_children(&vgsHierarchy.nodes[i])) ++sample.tagged_rendered;
    int centroid_count = cube_dynamic_only ? sample.tagged_dynamic : sample.tagged_voxels;
    if (centroid_count > 0) {
        sample.centroid_y = centroid_sum / (float)centroid_count;
        sample.mean_velocity_y = velocity_sum / (float)centroid_count;
    } else {
        sample.min_y = 0.0f;
        sample.max_y = 0.0f;
    }
    sample.centroid_drop = initial_centroid_y - sample.centroid_y;

    float fluid_y_sum = 0.0f;
    float density_ratio_sum = 0.0f;
    float speed_sum = 0.0f;
    for (int i = 0; i < fluid_particle_count; ++i) {
        Particle *particle = fluid_particles[i];
        ++sample.fluid_particles;
        if (!particle || !particle->active || particle->material != PARTICLE_MATERIAL_FLUID ||
            particle->fluid_index != i || !v_isfinite(particle->pos) ||
            !v_isfinite(particle->vel) || !isfinite(particle->pbf_density)) {
            ++sample.invalid_fluid_particles;
            continue;
        }
        float density_ratio = particle->pbf_density / PBF_REST_DENSITY;
        float speed = v_length(particle->vel);
        fluid_y_sum += particle->pos.y;
        density_ratio_sum += density_ratio;
        speed_sum += speed;
        if (density_ratio > sample.max_density_ratio) {
            sample.max_density_ratio = density_ratio;
        }
    }
    int valid_fluid = sample.fluid_particles - sample.invalid_fluid_particles;
    if (valid_fluid > 0) {
        float inverse_count = 1.0f / (float)valid_fluid;
        sample.fluid_centroid_y = fluid_y_sum * inverse_count;
        sample.mean_density_ratio = density_ratio_sum * inverse_count;
        sample.mean_fluid_speed = speed_sum * inverse_count;
    }
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
    if (state->active_cube_drop && debugColorVoxelsByVgsTemporalCurvature) box_height += 30;
    if (state->refinement_demo || state->coarsening_demo) box_height += 30;
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
    if (state->active_cube_drop && debugColorVoxelsByVgsTemporalCurvature) {
        int y = state->gpu_readback_performed ? 222 : 192;
        DrawText("VGS curvature (s^-2)", 34, y, 17, RAYWHITE);
        DrawRectangle(270, y + 2, 24, 16, ColorFromHSV(120.0f, 0.85f, 0.9f));
        DrawText("0", 298, y, 17, LIGHTGRAY);
        DrawRectangle(328, y + 2, 24, 16, ColorFromHSV(60.0f, 0.85f, 0.9f));
        DrawText("50", 356, y, 17, LIGHTGRAY);
        DrawRectangle(421, y + 2, 24, 16, ColorFromHSV(0.0f, 0.85f, 0.9f));
        DrawText("100+", 449, y, 17, LIGHTGRAY);
    }
    if (state->refinement_demo || state->coarsening_demo) {
        int y = state->gpu_readback_performed ? 222 : 192;
        bool coarse = state->refinement_demo ? !state->refinement_complete :
                      state->coarsening_complete;
        DrawText(TextFormat("VGS %s | %s at step %d",
                            coarse ? "ROOT ONLY" :
                            (state->coarsening_demo ? "8 CHILDREN" : "ROOT + 8 CHILDREN"),
                            state->coarsening_demo ? "coarsen" : "refine",
                            state->refinement_step), 34, y, 17, RAYWHITE);
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
    if (state->pbf_container) {
        int y = 192;
        DrawText(TextFormat("PBF particles=%d invalid=%d centroidY=%.3f meanSpeed=%.3f",
                            sample->fluid_particles, sample->invalid_fluid_particles,
                            sample->fluid_centroid_y, sample->mean_fluid_speed),
                 34, y, 17, result_color);
        DrawText(TextFormat("PBF density/rest mean=%.3f max=%.3f render=%s",
                            sample->mean_density_ratio, sample->max_density_ratio,
                            fluidRenderMode == FLUID_RENDER_SURFACE ? "surface" : "particles"),
                 34, y + 26, 17, result_color);
        if (fluidRenderMode == FLUID_RENDER_SURFACE) {
            DrawText(TextFormat("HULL cells=%d triangles=%d grid=%.3fm",
                                fluidHullCellCount, fluidHullMesh.triangleCount,
                                PBF_HULL_CELL_SIZE),
                     34, y + 52, 17, result_color);
        }
        DrawText(TextFormat("DYNAMIC FLUID CONTACTS %d",
                            atomic_load_explicit(&fluidDynamicCollisionCount,
                                                 memory_order_relaxed)),
                 34, y + 78, 17, result_color);
    }
}

static void debug_render_scene(DebugRunState *state, const DebugSample *sample)
{
    voxel_render_frame++;
    BeginTextureMode(state->target);
        ClearBackground((Color){ 32, 42, 58, 255 });
        BeginMode3D(state->camera);
            DrawPlane((Vector3){ 0.0f, 0.0f, 0.0f }, (Vector2){ 24.0f, 24.0f },
                      (Color){ 68, 75, 82, 255 });
            DrawGrid(48, VOXEL_SIZE);
            DrawVoxels(state->camera);
            if (state->adaptive_cube_drop) {
                for (int level = 0; level < vgsHierarchy.levels; ++level) {
                    int first = vgsHierarchy.level_start[level];
                    int end = first + vgsHierarchy.level_count[level];
                    Color color = ColorFromHSV(235.0f *
                        (1.0f - (float)level / (float)vgsHierarchy.levels), 0.8f, 0.9f);
                    for (int i = first; i < end; ++i) {
                        VgsHierarchyNode *node = &vgsHierarchy.nodes[i];
                        if (!node->active || !node->vgs_enabled ||
                            vgs_node_has_active_children(node)) continue;
                        Vector3 center = {0};
                        for (int c = 0; c < 8; ++c) center = v_add(center, node->particles[c]->pos);
                        center = v_mul(center, 0.125f);
                        Vector3 x = v_sub(node->particles[1]->pos, node->particles[0]->pos);
                        Vector3 y = v_sub(node->particles[2]->pos, node->particles[0]->pos);
                        Vector3 z = v_sub(node->particles[4]->pos, node->particles[0]->pos);
                        Matrix transform = MatrixIdentity();
                        transform.m0 = x.x; transform.m1 = x.y; transform.m2 = x.z;
                        transform.m4 = y.x; transform.m5 = y.y; transform.m6 = y.z;
                        transform.m8 = z.x; transform.m9 = z.y; transform.m10 = z.z;
                        transform.m12 = center.x; transform.m13 = center.y; transform.m14 = center.z;
                        rlPushMatrix();
                        rlMultMatrixf(MatrixToFloat(transform));
                        DrawCubeV((Vector3){0}, (Vector3){1, 1, 1}, color);
                        if (level < 3)
                            DrawCubeWiresV((Vector3){0}, (Vector3){1, 1, 1}, DARKGRAY);
                        rlPopMatrix();
                    }
                }
            }
            if ((state->refinement_demo || state->coarsening_demo) &&
                debugRefinementCoarse &&
                vgsHierarchy.node_count == 1) {
                const VgsHierarchyNode *root = &vgsHierarchy.nodes[0];
                Vector3 center = { 0 };
                for (int c = 0; c < 8; ++c)
                    center = v_add(center, root->particles[c]->pos);
                center = v_mul(center, 0.125f);
                Vector3 size = {
                    v_length(v_sub(root->particles[1]->pos, root->particles[0]->pos)),
                    v_length(v_sub(root->particles[2]->pos, root->particles[0]->pos)),
                    v_length(v_sub(root->particles[4]->pos, root->particles[0]->pos))
                };
                DrawCubeV(center, size, (Color){ 40, 185, 220, 255 });
                DrawCubeWiresV(center, size, (Color){ 15, 75, 100, 255 });
            }
            if (state->swarmer_wall_test || state->goliath_flight_test) {
                draw_players();
            }
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
    int expected_total = state->expected_voxels + state->expected_static_leftovers;
    if (sample->tagged_voxels != expected_total) {
        debug_add_failure(state, "%s: expected %d tagged voxels, found %d",
                          context, expected_total, sample->tagged_voxels);
    }
    if (sample->tagged_dynamic != state->expected_voxels ||
        sample->tagged_static != state->expected_static_leftovers) {
        debug_add_failure(state, "%s: tagged counts dynamic=%d/%d static=%d/%d",
                          context, sample->tagged_dynamic, state->expected_voxels,
                          sample->tagged_static, state->expected_static_leftovers);
    }
    if (sample->tagged_sleeping != 0) {
        debug_add_failure(state, "%s: %d tagged voxels remain sleeping", context,
                          sample->tagged_sleeping);
    }
    if (sample->invalid_particle_refs > 0) {
        debug_add_failure(state, "%s: %d invalid particle or sim-index references",
                          context, sample->invalid_particle_refs);
    }
    if (sample->zero_mass_corners > 0 &&
        !(state->refinement_demo && !state->refinement_complete) &&
        !(state->coarsening_demo && state->coarsening_complete)) {
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
    // Adaptation follows the physics dispatch. Its new sim list is staged for
    // the next frame, so compare the lists only after that upload has happened.
    if (state->adaptive_cube_drop &&
        gpuPhysics.hierarchy_generation != vgsHierarchy.topology_generation) return;
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
            if (state->adaptive_cube_drop && particle->base_inv_mass <= 0.0f) continue;
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
    if (!state->overhang_impact &&
        !(state->refinement_demo && !state->refinement_complete) &&
        !(state->coarsening_demo && state->coarsening_complete) &&
        state->gpu_zero_mass_corners > 0) {
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

static bool debug_setup_cube_stress(DebugRunState *state)
{
    debug_reset_world();
    int size = 8;
    const char *size_env = getenv("FPS_CUBE_SIZE");
    if (size_env) {
        int parsed = atoi(size_env);
        if (parsed >= 2 && parsed <= 50) size = parsed;
    }
    int half = size / 2;
    int created = 0;
    Color color = (Color){ 245, 145, 45, 255 };
    for (int y = 20; y < 20 + size; ++y) {
        for (int z = -half; z < -half + size; ++z) {
            for (int x = -half; x < -half + size; ++x) {
                if (debug_add_cell(x, y, z, false, true, color, DEBUG_TEST_TAG) >= 0) {
                    ++created;
                }
            }
        }
    }
    state->expected_voxels = created;
    debug_track_tagged_cells(state, DEBUG_TEST_TAG);
    debug_rebuild_world_state();
    state->activation_succeeded = true;
    debug_set_initial_centroid(state);
    float dist = (float)size * 2.5f + 10.0f;
    state->camera.position = (Vector3){ dist, (float)(20 + size/2) + dist * 0.4f, dist };
    state->camera.target = (Vector3){ 0.0f, (float)(20 + size/2), 0.0f };
    return true;
}

static int debug_compare_curvature(const void *a, const void *b)
{
    float x = *(const float *)a, y = *(const float *)b;
    return (x > y) - (x < y);
}

static void debug_write_curvature_levels(DebugRunState *state, int step)
{
    if (!state->curvature_levels_file || !vgsHierarchy.leaf_history) return;
    int capacity = vgsHierarchy.leaf_count > vgsHierarchy.node_count
        ? vgsHierarchy.leaf_count : vgsHierarchy.node_count;
    float *values = malloc((size_t)capacity * sizeof(*values));
    if (!values) {
        debug_add_failure(state, "unable to allocate curvature distribution scratch");
        return;
    }
    for (int level = 0; level <= vgsHierarchy.levels; ++level) {
        int count = 0, above = 0, child_above = 0;
        int above_adaptive = 0, child_above_adaptive = 0;
        double child_sum = 0.0;
        float level_threshold = ldexpf(VGS_ADAPTIVE_CURVATURE_THRESHOLD,
                                      (VGS_ADAPTIVE_THRESHOLD_LEVEL_FACTOR == 4 ? 2 : 1) * level);
        if (level < vgsHierarchy.levels) {
            int first = vgsHierarchy.level_start[level];
            int end = first + vgsHierarchy.level_count[level];
            for (int i = first; i < end; ++i) {
                float value = vgsHierarchy.node_history[i].curvature;
                values[count++] = value;
                above += value > VGS_TEMPORAL_CURVATURE_RED_AT;
                above_adaptive += value > level_threshold;
                VgsHierarchyNode *node = &vgsHierarchy.nodes[i];
                float mean = 0.0f;
                for (int c = 0; c < 8; ++c) {
                    int child = node->children[c];
                    mean += child >= 0 ? vgsHierarchy.node_history[child].curvature :
                        vgsHierarchy.leaf_history[-child - 1].curvature;
                }
                mean *= 0.125f;
                child_sum += mean;
                child_above += mean > VGS_TEMPORAL_CURVATURE_RED_AT;
                child_above_adaptive += mean >
                    VGS_ADAPTIVE_THRESHOLD_LEVEL_FACTOR * level_threshold;
            }
        } else {
            for (int i = 0; i < voxel_count; ++i) {
                if (!voxels[i].simulate || voxels[i].debugClusterTag != DEBUG_TEST_TAG) continue;
                float value = vgsHierarchy.leaf_history[i].curvature;
                values[count++] = value;
                above += value > VGS_TEMPORAL_CURVATURE_RED_AT;
                above_adaptive += value > level_threshold;
            }
        }
        if (count > capacity) { debug_add_failure(state, "curvature level overflow"); break; }
        if (!count) continue;
        qsort(values, (size_t)count, sizeof(*values), debug_compare_curvature);
        fprintf(state->curvature_levels_file,
                "%d,%d,%d,%.9g,%.9g,%.9g,%.9g,%d,%.9g,%d,%d,%d\n",
                step, level, count,
                (double)values[(count - 1) / 2],
                (double)values[(int)((count - 1) * 0.90)],
                (double)values[(int)((count - 1) * 0.99)],
                (double)values[count - 1], above,
                count ? child_sum / count : 0.0, child_above,
                above_adaptive, child_above_adaptive);
    }
    free(values);
}

static bool debug_setup_active_cube_drop(DebugRunState *state)
{
    debug_reset_world();
    state->active_cube_drop = true;
    state->expected_voxels = 10 * 10 * 10;
    int created = 0;
    const Color color = { 245, 145, 45, 255 };
    for (int y = 18; y < 28; ++y) {
        for (int z = -5; z < 5; ++z) {
            for (int x = -5; x < 5; ++x) {
                if (debug_add_cell(x, y, z, true, false, color, DEBUG_TEST_TAG) >= 0) {
                    ++created;
                }
            }
        }
    }
    if (created != state->expected_voxels) {
        debug_add_failure(state, "active-cube-drop created %d/%d voxels",
                          created, state->expected_voxels);
        return false;
    }
    debug_track_tagged_cells(state, DEBUG_TEST_TAG);
    debug_rebuild_world_state();
    state->camera.position = (Vector3){ 18.0f, 16.0f, 20.0f };
    state->camera.target = (Vector3){ 0.0f, 7.0f, 0.0f };
    debug_capture(state, "phase-01-static-cube.png", "static-cube", 0);
    if (!debug_convert_cube_to_pbd(state, DEBUG_PBD_VOXEL_SIZE)) return false;
    if (!debug_init_vgs_temporal_curvature()) {
        debug_add_failure(state, "unable to allocate VGS temporal history");
        return false;
    }
    debugColorVoxelsByVgsTemporalCurvature = true;
    debug_sample_vgs_temporal_curvature(PBD_MAX_STEP_DT);
    if (!debug_join_path(state->curvature_levels_path, sizeof(state->curvature_levels_path),
                         state->output_dir, "curvature-by-level.csv")) return false;
    state->curvature_levels_file = fopen(state->curvature_levels_path, "w");
    if (!state->curvature_levels_file) return false;
    fputs("step,level,count,p50,p90,p99,max,above100,meanChildCurvature,childMeansAbove100,aboveAdaptiveThreshold,childMeansAboveAdaptiveThreshold\n",
          state->curvature_levels_file);
    debug_write_curvature_levels(state, 0);
    debug_set_initial_centroid(state);
    return true;
}

static bool debug_setup_adaptive_cube_drop(DebugRunState *state)
{
    state->adaptive_cube_drop = true;
    if (!debug_setup_active_cube_drop(state)) return false;
    if (getenv("FPS_AMR_CORNER_DROP") &&
        !debug_tilt_tagged_cube_for_corner_drop(DEBUG_TEST_TAG)) {
        debug_add_failure(state, "unable to tilt adaptive cube");
        return false;
    }
    if (!vgs_hierarchy_initialize_adaptive()) {
        debug_add_failure(state, "unable to initialize adaptive octree");
        return false;
    }
    // A diagnostic seed exercises a refined GPU constraint list even when the
    // normal curvature trigger has not fired yet.
    if (getenv("FPS_AMR_TEST_REFINE_ROOT") || getenv("FPS_AMR_TEST_BREAK_CHILD") ||
        getenv("FPS_AMR_TEST_BREAK_PARENT")) {
        vgs_hierarchy_transfer_refine(&vgsHierarchy.nodes[0], true);
        if (!vgs_hierarchy_rebuild_active_masses()) return false;
    }
    if (getenv("FPS_AMR_TEST_BREAK_CHILD")) {
        Particle *grid[27];
        if (!vgs_hierarchy_grid_points(&vgsHierarchy.nodes[0], grid)) return false;
        Particle *center = grid[13];
        center->pos.y += 3.5f;
        center->prev_pos = center->predicted_pos = center->pos;
    }
    if (getenv("FPS_AMR_TEST_BREAK_PARENT"))
        vgsHierarchy.nodes[0].vgs_enabled = false;
    debugAdaptiveOctree = true;
    debugColorVoxelsByVgsTemporalCurvature = false;
    return true;
}

static void debug_apply_rip_grips(DebugRunState *state)
{
    const float rest_edge = state->rip_right_x - state->rip_left_x;
    const int total_substeps = debugRunConfig.steps * PBD_SUBSTEPS;
    state->rip_commanded_strain = 2.0f * STRAIN_BREAK_THRESHOLD *
        (float)state->rip_substeps / (float)total_substeps;
    const float half_extension = 0.5f * rest_edge * state->rip_commanded_strain;
    const float rest_step = 0.5f * pbdSolidVoxelSize;
    for (int i = 0; i < particle_pool_count; ++i) {
        Particle *particle = &particles_pool[i];
        if (!particle->active) continue;
        float rest_x = (float)particle->rest_cell_x * rest_step;
        bool left = rest_x <= state->rip_left_x + state->rip_grip_depth + 1e-5f;
        bool right = rest_x >= state->rip_right_x - state->rip_grip_depth - 1e-5f;
        if (!left && !right) continue;
        float target_x = rest_x + (left ? -half_extension : half_extension);
        particle->pos.x = target_x;
        particle->prev_pos.x = target_x;
        particle->predicted_pos.x = target_x;
        particle->vel.x = 0.0f;
    }
    collisionTopologyDirty = true;
}

static bool debug_setup_rip_test(DebugRunState *state)
{
    state->adaptive_cube_drop = true;
    state->rip_test = true;
    state->rip_first_break_step = -1;
    if (!debug_setup_active_cube_drop(state)) return false;
    if (!vgs_hierarchy_initialize_adaptive()) return false;
    vgs_hierarchy_transfer_refine(&vgsHierarchy.nodes[0], true);
    if (!vgs_hierarchy_rebuild_active_masses()) return false;
    debugAdaptiveOctree = true;
    debugColorVoxelsByVgsTemporalCurvature = false;
    debugRipTest = true;
    state->rip_left_x = FLT_MAX;
    state->rip_right_x = -FLT_MAX;
    for (int c = 0; c < 8; ++c) {
        Particle *particle = vgsHierarchy.nodes[0].particles[c];
        if (particle->pos.x < state->rip_left_x) state->rip_left_x = particle->pos.x;
        if (particle->pos.x > state->rip_right_x) state->rip_right_x = particle->pos.x;
    }
    if (state->rip_right_x - state->rip_left_x < 4.9f) return false;
    state->rip_grip_depth = 0.25f * (state->rip_right_x - state->rip_left_x);
    state->camera.position = (Vector3){ 0.0f, 15.0f, 24.0f };
    state->camera.target = (Vector3){ 0.0f, 11.5f, 0.0f };
    if (!debug_join_path(state->rip_history_path, sizeof(state->rip_history_path),
                         state->output_dir, "rip-history.csv")) return false;
    state->rip_history_file = fopen(state->rip_history_path, "w");
    if (!state->rip_history_file) return false;
    fputs("step,commandedStrain,leftTargetX,rightTargetX,leftGripParticles,rightGripParticles,"
          "leftSimulatedGripParticles,rightSimulatedGripParticles,maxGripError,"
          "maxActiveVgsStrain,brokenNodes,brokenLeaves,enabledVgs,enabledWelds,"
          "refinements,coarsenings\n", state->rip_history_file);
    if (!debug_join_path(state->rip_strain_levels_path,
                         sizeof(state->rip_strain_levels_path),
                         state->output_dir, "rip-parent-child-strain.csv")) return false;
    state->rip_strain_levels_file = fopen(state->rip_strain_levels_path, "w");
    if (!state->rip_strain_levels_file) return false;
    fputs("step,parentLevel,relationships,parentMeanStrain,childMeanStrain,"
          "ratioOfMeans,meanNodeRatio,intactRelationships,intactParentMeanStrain,"
          "intactChildMeanStrain,intactRatioOfMeans\n",
          state->rip_strain_levels_file);
    debug_apply_rip_grips(state);
    return true;
}

static float debugRefineWeights[27][8];
static float debugCoarsenWeights[8][27];

static void debug_build_refinement_weights(void)
{
    static const float t[3] = { 0.0f, 0.5f, 1.0f };
    for (int gz = 0; gz < 3; ++gz) for (int gy = 0; gy < 3; ++gy) for (int gx = 0; gx < 3; ++gx) {
        int r = gx + 3 * gy + 9 * gz;
        for (int c = 0; c < 8; ++c) {
            float wx = (c & 1) ? t[gx] : 1.0f - t[gx];
            float wy = (c & 2) ? t[gy] : 1.0f - t[gy];
            float wz = (c & 4) ? t[gz] : 1.0f - t[gz];
            debugRefineWeights[r][c] = wx * wy * wz;
        }
    }
}

// For this fixed 3x3x3 lattice, A is the trilinear 27x8 refinement matrix.
// Build its left pseudoinverse once, then reuse it for positions and velocities.
static bool debug_build_coarsening_weights(void)
{
    double augmented[8][16] = { 0 };
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j)
            for (int r = 0; r < 27; ++r)
                augmented[i][j] += (double)debugRefineWeights[r][i] *
                                   (double)debugRefineWeights[r][j];
        augmented[i][8 + i] = 1.0;
    }
    for (int i = 0; i < 8; ++i) {
        int pivot = i;
        for (int row = i + 1; row < 8; ++row)
            if (fabs(augmented[row][i]) > fabs(augmented[pivot][i])) pivot = row;
        if (fabs(augmented[pivot][i]) < 1e-12) return false;
        if (pivot != i) for (int j = 0; j < 16; ++j) {
            double tmp = augmented[i][j];
            augmented[i][j] = augmented[pivot][j];
            augmented[pivot][j] = tmp;
        }
        double inverse_pivot = 1.0 / augmented[i][i];
        for (int j = 0; j < 16; ++j) augmented[i][j] *= inverse_pivot;
        for (int row = 0; row < 8; ++row) {
            if (row == i) continue;
            double factor = augmented[row][i];
            for (int j = 0; j < 16; ++j)
                augmented[row][j] -= factor * augmented[i][j];
        }
    }
    for (int i = 0; i < 8; ++i) for (int r = 0; r < 27; ++r) {
        double weight = 0.0;
        for (int j = 0; j < 8; ++j)
            weight += augmented[i][8 + j] * (double)debugRefineWeights[r][j];
        debugCoarsenWeights[i][r] = (float)weight;
    }
    // R*A must reproduce every possible root-corner state.
    for (int i = 0; i < 8; ++i) for (int j = 0; j < 8; ++j) {
        double value = 0.0;
        for (int r = 0; r < 27; ++r)
            value += (double)debugCoarsenWeights[i][r] * debugRefineWeights[r][j];
        if (fabs(value - (i == j ? 1.0 : 0.0)) > 1e-5) return false;
    }
    return true;
}

static bool debug_sync_refinement_grid(DebugRunState *state, bool refine)
{
    if (!state || (!state->refinement_demo && !state->coarsening_demo) ||
        vgsHierarchy.node_count != 1) return false;
    VgsHierarchyNode *root = &vgsHierarchy.nodes[0];
    Vector3 parent_pos[8], parent_vel[8];
    for (int c = 0; c < 8; ++c) {
        Particle *particle = root->particles[c];
        int r = 2 * (c & 1) + 6 * ((c >> 1) & 1) + 18 * ((c >> 2) & 1);
        if (!particle || state->refine_grid[r] != particle) return false;
        parent_pos[c] = particle->pos;
        parent_vel[c] = particle->vel;
    }
    for (int r = 0; r < 27; ++r) {
        Particle *particle = state->refine_grid[r];
        if (!particle) return false;
        int gx = r % 3, gy = (r / 3) % 3, gz = r / 9;
        bool parent_corner = gx != 1 && gy != 1 && gz != 1;
        if (!parent_corner) {
            Vector3 pos = { 0 }, vel = { 0 };
            for (int c = 0; c < 8; ++c) {
                float w = debugRefineWeights[r][c];
                pos = v_add(pos, v_mul(parent_pos[c], w));
                vel = v_add(vel, v_mul(parent_vel[c], w));
            }
            particle->pos = particle->prev_pos = particle->predicted_pos = pos;
            particle->vel = vel;
        }
        particle->inv_mass = particle->base_inv_mass = refine ?
            state->refine_inv_mass[r] : (parent_corner ? 0.001f : 0.0f);
    }
    for (int c = 0; c < 8; ++c) {
        Voxel *child = &voxels[state->refine_children[c]];
        child->simulate_dofs = refine;
        Vector3 center = { 0 }, velocity = { 0 };
        for (int k = 0; k < 8; ++k) {
            center = v_add(center, child->particles[k]->pos);
            velocity = v_add(velocity, child->particles[k]->vel);
        }
        child->pos = v_mul(center, 0.125f);
        child->vel = v_mul(velocity, 0.125f);
    }
    return true;
}

static bool debug_coarsen_grid(DebugRunState *state)
{
    if (!state || !state->coarsening_demo || vgsHierarchy.node_count != 1) return false;
    Vector3 fine_pos[27], fine_vel[27], coarse_pos[8], coarse_vel[8];
    for (int r = 0; r < 27; ++r) {
        Particle *particle = state->refine_grid[r];
        if (!particle) return false;
        fine_pos[r] = particle->pos;
        fine_vel[r] = particle->vel;
    }
    for (int c = 0; c < 8; ++c) {
        coarse_pos[c] = coarse_vel[c] = (Vector3){ 0 };
        for (int r = 0; r < 27; ++r) {
            float weight = debugCoarsenWeights[c][r];
            coarse_pos[c] = v_add(coarse_pos[c], v_mul(fine_pos[r], weight));
            coarse_vel[c] = v_add(coarse_vel[c], v_mul(fine_vel[r], weight));
        }
    }
    double error_squared = 0.0;
    for (int r = 0; r < 27; ++r) {
        Vector3 fitted = { 0 };
        for (int c = 0; c < 8; ++c)
            fitted = v_add(fitted, v_mul(coarse_pos[c], debugRefineWeights[r][c]));
        Vector3 error = v_sub(fitted, fine_pos[r]);
        error_squared += (double)v_dot(error, error);
    }
    state->coarsening_projection_rms = (float)sqrt(error_squared / 27.0);
    if (!isfinite(state->coarsening_projection_rms)) return false;
    for (int c = 0; c < 8; ++c) {
        Particle *particle = vgsHierarchy.nodes[0].particles[c];
        int r = 2 * (c & 1) + 6 * ((c >> 1) & 1) + 18 * ((c >> 2) & 1);
        if (!particle || particle != state->refine_grid[r]) return false;
        particle->pos = particle->prev_pos = particle->predicted_pos = coarse_pos[c];
        particle->vel = coarse_vel[c];
    }
    vgsHierarchy.nodes[0].active = true;
    return debug_sync_refinement_grid(state, false);
}

static bool debug_setup_vgs_refinement(DebugRunState *state)
{
    debug_reset_world();
    state->refinement_demo = true;
    state->refinement_step = DEBUG_REFINE_AT_STEP;
    const Color color = { 245, 145, 45, 255 };
    int created = 0;
    for (int y = 18; y < 28; ++y) for (int z = -5; z < 5; ++z) for (int x = -5; x < 5; ++x)
        if (debug_add_cell(x, y, z, true, false, color, DEBUG_TEST_TAG) >= 0) ++created;
    if (created != 1000) {
        debug_add_failure(state, "refinement demo created %d/1000 static cells", created);
        return false;
    }
    debug_track_tagged_cells(state, DEBUG_TEST_TAG);
    debug_rebuild_world_state();
    if (!debug_convert_cube_to_pbd(state, 2.5f)) return false;
    debug_build_refinement_weights();
    if (!debug_sync_refinement_grid(state, false)) {
        debug_add_failure(state, "refinement demo could not initialize the shared grid");
        return false;
    }
    debugRefinementCoarse = true;
    state->camera.position = (Vector3){ 34.0f, 31.0f, 35.0f };
    state->camera.target = (Vector3){ 0.0f, 15.0f, 0.0f };
    debug_set_initial_centroid(state);
    return true;
}

static bool debug_setup_vgs_coarsening(DebugRunState *state)
{
    debug_reset_world();
    state->coarsening_demo = true;
    state->refinement_step = DEBUG_REFINE_AT_STEP;
    const Color color = { 245, 145, 45, 255 };
    int created = 0;
    for (int y = 18; y < 28; ++y) for (int z = -5; z < 5; ++z) for (int x = -5; x < 5; ++x)
        if (debug_add_cell(x, y, z, true, false, color, DEBUG_TEST_TAG) >= 0) ++created;
    if (created != 1000) {
        debug_add_failure(state, "coarsening demo created %d/1000 static cells", created);
        return false;
    }
    debug_track_tagged_cells(state, DEBUG_TEST_TAG);
    debug_rebuild_world_state();
    if (!debug_convert_cube_to_pbd(state, 2.5f)) return false;
    debug_build_refinement_weights();
    if (!debug_build_coarsening_weights()) {
        debug_add_failure(state, "coarsening pseudoinverse is singular or inaccurate");
        return false;
    }
    vgsHierarchy.nodes[0].active = false;
    state->camera.position = (Vector3){ 20.0f, 18.0f, 22.0f };
    state->camera.target = (Vector3){ 0.0f, 6.0f, 0.0f };
    debug_set_initial_centroid(state);
    return true;
}

static bool debug_setup_voxel_strain_shear_fracture(DebugRunState *state)
{
    debug_reset_world();
    state->voxel_fracture_test = true;
    Color strain_color = (Color){ 245, 120, 60, 255 };
    Color shear_color = (Color){ 220, 70, 180, 255 };
    Color control_color = (Color){ 80, 200, 120, 255 };

    int idx_strain = debug_add_cell(-3, 8, 0, false, true, strain_color, DEBUG_TEST_TAG);
    int idx_shear = debug_add_cell(0, 8, 0, false, true, shear_color, DEBUG_TEST_TAG);
    int idx_ctrl = debug_add_cell(3, 8, 0, false, true, control_color, DEBUG_TEST_TAG);

    if (idx_strain < 0 || idx_shear < 0 || idx_ctrl < 0) {
        debug_add_failure(state, "fracture setup could not create test voxels");
        return false;
    }
    debug_track_tagged_cells(state, DEBUG_TEST_TAG);
    debug_rebuild_world_state();

    // 1. High Tensile Strain: stretch X coordinates by 35% (limit is 20%)
    Voxel *v_strain = &voxels[idx_strain];
    float stretch_x = VOXEL_SIZE * 0.35f * 0.5f;
    for (int corner = 0; corner < VOXEL_CORNER_COUNT; ++corner) {
        Particle *p = v_strain->particles[corner];
        if (corner_signs[corner][0] > 0) p->pos.x += stretch_x;
        else p->pos.x -= stretch_x;
        p->prev_pos = p->pos;
        p->predicted_pos = p->pos;
    }

    // 2. High Shear: skew X based on Y sign by 35% (limit is 20%)
    Voxel *v_shear = &voxels[idx_shear];
    float skew = VOXEL_SIZE * 0.35f * 0.5f;
    for (int corner = 0; corner < VOXEL_CORNER_COUNT; ++corner) {
        Particle *p = v_shear->particles[corner];
        if (corner_signs[corner][1] > 0) p->pos.x += skew;
        else p->pos.x -= skew;
        p->prev_pos = p->pos;
        p->predicted_pos = p->pos;
    }

    state->fracture_strain_voxel = idx_strain;
    state->fracture_shear_voxel = idx_shear;
    state->fracture_control_voxel = idx_ctrl;
    state->expected_voxels = 3;
    state->activation_succeeded = true;
    collisionTopologyDirty = true;
    rebuild_voxel_hash();
    debug_set_initial_centroid(state);
    state->camera.position = (Vector3){ 0.0f, 9.0f, 7.0f };
    state->camera.target = (Vector3){ 0.0f, 8.0f, 0.0f };
    debug_capture(state, "phase-01-fracture-setup.png", "fracture-setup", 0);
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
    debugLogVoxelDeactivation = true;
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
    int thrown_index = debug_add_cell(-5, 2, 0, false, true,
                                      (Color){ 255, 195, 55, 255 }, DEBUG_THROW_TAG);
    debug_rebuild_world_state();

    if (cluster_created != state->expected_voxels || thrown_index < 0) {
        debug_add_failure(state, "sleep-snap setup created pillar=%d/%d throw=%d",
                          cluster_created, state->expected_voxels, thrown_index);
        return false;
    }

    Voxel *thrown = &voxels[thrown_index];
    thrown->owner = 0;
    thrown->activator = 0;
    thrown->wasTethered = true;
    // Use the tether's direct-activation strength so this isolated grounded
    // pillar reliably wakes when struck, matching a gameplay structure whose
    // local support has already been weakened by damage.
    thrown->activationBelief = ACTIVATION_TETHER_BELIEF;
    players[0].yaw = -90.0f;
    players[0].pitch = 0.0f;
    players[0].tetherHolding = true;
    players[0].tetherVoxel = thrown_index;
    players[0].tetherVoxelIdentity = thrown->identity;
    release_tether(0);
    state->tether_thin_wall_hits_before = tetherThrowCcdHitCount;
    state->camera.position = (Vector3){ 6.5f, 4.5f, 8.0f };
    state->camera.target = (Vector3){ 0.0f, 1.5f, 0.0f };
    debug_set_initial_centroid(state);
    debug_capture(state, "phase-01-tether-thrown.png", "tether-thrown-static-pillar", 0);
    return true;
}

static bool debug_setup_tall_pillar_snap(DebugRunState *state)
{
    debug_reset_world();
    state->sleep_snap_original_test = true;
    debugLogVoxelDeactivation = true;
    state->sleep_snap_static_step = -1;
    const int height = 14;
    const int split_y = 7;
    int expected_test = 2 * 2 * (height - split_y);
    state->expected_voxels = expected_test;
    Color static_color = (Color){ 80, 175, 245, 255 };
    Color dynamic_color = (Color){ 245, 145, 45, 255 };

    // Bottom half static
    for (int y = 0; y < split_y; ++y) {
        for (int z = 0; z < 2; ++z) {
            for (int x = 0; x < 2; ++x) {
                debug_add_cell(x, y, z, true, false, static_color, DEBUG_SUPPORT_TAG);
            }
        }
    }
    // Top half dynamic (activated)
    for (int y = split_y; y < height; ++y) {
        for (int z = 0; z < 2; ++z) {
            for (int x = 0; x < 2; ++x) {
                int idx = debug_add_cell(x, y, z, false, true, dynamic_color, DEBUG_TEST_TAG);
                if (idx >= 0) {
                    voxels[idx].vel = (Vector3){ 0.8f, 0.0f, 0.0f };
                }
            }
        }
    }
    debug_track_tagged_cells(state, DEBUG_TEST_TAG);
    debug_rebuild_world_state();
    for (int i = 0; i < voxel_count; ++i) {
        if (voxels[i].debugClusterTag == DEBUG_TEST_TAG && voxels[i].gy == split_y) {
            break_face_link(&voxels[i], 3);
        }
    }
    for (int i = 0; i < voxel_count; ++i) {
        if (voxels[i].debugClusterTag == DEBUG_TEST_TAG && voxels[i].simulate) {
            set_voxel_velocity(&voxels[i], (Vector3){ 1.5f, 0.0f, 0.0f });
        }
    }
    state->camera.position = (Vector3){ 6.5f, 5.5f, 8.0f };
    state->camera.target = (Vector3){ 0.5f, 3.5f, 0.5f };
    debug_set_initial_centroid(state);
    state->activation_succeeded = true;
    state->sleep_snap_active_observed = true;
    return true;
}

static bool debug_setup_sleep_snap_no_rebuild(DebugRunState *state)
{
    if (!debug_setup_sleep_snap_original(state)) return false;
    state->sleep_snap_rebuild_test = true;
    int minimum_steps = RECYCLE_STATIC_RESTORE_DELAY + 200;
    if (debugRunConfig.steps < minimum_steps) debugRunConfig.steps = minimum_steps;
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

static bool debug_setup_pbf_container(DebugRunState *state)
{
    debug_reset_world();
    activePlayers = 0;
    const Color container_color = { 125, 135, 150, 255 };
    const Color fluid_color = { 45, 145, 235, 210 };
    int container_created = 0;
    int fluid_cells_created = 0;

    for (int z = -4; z <= 4; ++z) {
        for (int x = -4; x <= 4; ++x) {
            if (debug_add_cell(x, 0, z, true, false,
                               container_color, DEBUG_CONTAINER_TAG) >= 0) {
                ++container_created;
            }
        }
    }
    for (int y = 1; y <= 20; ++y) {
        for (int z = -4; z <= 4; ++z) {
            if (debug_add_cell(-4, y, z, true, false,
                               container_color, DEBUG_CONTAINER_TAG) >= 0) {
                ++container_created;
            }
            if (debug_add_cell(4, y, z, true, false,
                               container_color, DEBUG_CONTAINER_TAG) >= 0) {
                ++container_created;
            }
        }
        for (int x = -3; x <= 3; ++x) {
            if (debug_add_cell(x, y, -4, true, false,
                               container_color, DEBUG_CONTAINER_TAG) >= 0) {
                ++container_created;
            }
            if (debug_add_cell(x, y, 4, true, false,
                               container_color, DEBUG_CONTAINER_TAG) >= 0) {
                ++container_created;
            }
        }
    }

    for (int y = 8; y <= 11; ++y) {
        for (int z = -2; z <= 1; ++z) {
            for (int x = -2; x <= 1; ++x) {
                float px = ((float)x + 0.5f) * VOXEL_SIZE;
                float py = ((float)y + 0.5f) * VOXEL_SIZE;
                float pz = ((float)z + 0.5f) * VOXEL_SIZE;
                int index = add_fluid_cell(px, py, pz, fluid_color);
                if (index >= 0) {
                    voxels[index].debugClusterTag = DEBUG_FLUID_TAG;
                    ++fluid_cells_created;
                }
            }
        }
    }

    debug_rebuild_world_state();
    state->pbf_container = true;
    state->activation_succeeded = true;
    state->expected_container_voxels = 721;
    state->expected_fluid_particles = 4 * 4 * 4 * VOXEL_CORNER_COUNT;
    state->container_min_x = -1.5f + PBF_PARTICLE_RADIUS;
    state->container_max_x = 2.0f - PBF_PARTICLE_RADIUS;
    state->container_min_z = -1.5f + PBF_PARTICLE_RADIUS;
    state->container_max_z = 2.0f - PBF_PARTICLE_RADIUS;
    state->container_floor_top = VOXEL_SIZE;
    state->camera.position = (Vector3){ 3.0f, 18.0f, 3.0f };
    state->camera.target = (Vector3){ 0.25f, 1.0f, 0.25f };

    DebugSample initial = debug_collect_sample("initial-fluid", 0, 0.0f);
    state->initial_fluid_centroid_y = initial.fluid_centroid_y;
    if (container_created != state->expected_container_voxels) {
        debug_add_failure(state, "pbf-container created %d/%d immutable container voxels",
                          container_created, state->expected_container_voxels);
    }
    if (fluid_cells_created != 64 ||
        fluid_particle_count != state->expected_fluid_particles) {
        debug_add_failure(state,
                          "pbf-container created cells=%d/64 particles=%d/%d",
                          fluid_cells_created, fluid_particle_count,
                          state->expected_fluid_particles);
    }
    return state->failure_count == 0;
}

static bool debug_setup_pbf_dynamic_container(DebugRunState *state)
{
    debug_reset_world();
    activePlayers = 0;
    const Color floor_color = { 225, 145, 55, 255 };
    const Color wall_color = { 125, 135, 150, 255 };
    const Color fluid_color = { 45, 145, 235, 210 };
    int container_created = 0;
    int fluid_cells_created = 0;

    // A 7x7 active floor is glued around its complete perimeter to an immutable rim.
    for (int z = -3; z <= 3; ++z) {
        for (int x = -3; x <= 3; ++x) {
            if (debug_add_cell(x, 4, z, false, true,
                               floor_color, DEBUG_CONTAINER_TAG) >= 0) {
                ++container_created;
            }
        }
    }
    for (int y = 4; y <= 16; ++y) {
        for (int z = -4; z <= 4; ++z) {
            if (debug_add_cell(-4, y, z, true, false,
                               wall_color, DEBUG_CONTAINER_TAG) >= 0) ++container_created;
            if (debug_add_cell(4, y, z, true, false,
                               wall_color, DEBUG_CONTAINER_TAG) >= 0) ++container_created;
        }
        for (int x = -3; x <= 3; ++x) {
            if (debug_add_cell(x, y, -4, true, false,
                               wall_color, DEBUG_CONTAINER_TAG) >= 0) ++container_created;
            if (debug_add_cell(x, y, 4, true, false,
                               wall_color, DEBUG_CONTAINER_TAG) >= 0) ++container_created;
        }
    }

    for (int y = 10; y <= 13; ++y) {
        for (int z = -2; z <= 1; ++z) {
            for (int x = -2; x <= 1; ++x) {
                float px = ((float)x + 0.5f) * VOXEL_SIZE;
                float py = ((float)y + 0.5f) * VOXEL_SIZE;
                float pz = ((float)z + 0.5f) * VOXEL_SIZE;
                int index = add_fluid_cell(px, py, pz, fluid_color);
                if (index >= 0) {
                    voxels[index].debugClusterTag = DEBUG_FLUID_TAG;
                    ++fluid_cells_created;
                }
            }
        }
    }

    debug_rebuild_world_state();
    state->pbf_container = true;
    state->pbf_dynamic_container = true;
    state->activation_succeeded = true;
    state->expected_container_voxels = 49 + 13 * 32;
    state->expected_fluid_particles = 4 * 4 * 4 * VOXEL_CORNER_COUNT;
    state->container_min_x = -3.0f * VOXEL_SIZE + PBF_PARTICLE_RADIUS;
    state->container_max_x = 4.0f * VOXEL_SIZE - PBF_PARTICLE_RADIUS;
    state->container_min_z = -3.0f * VOXEL_SIZE + PBF_PARTICLE_RADIUS;
    state->container_max_z = 4.0f * VOXEL_SIZE - PBF_PARTICLE_RADIUS;
    state->container_floor_top = 5.0f * VOXEL_SIZE;
    state->camera.position = (Vector3){ 3.0f, 6.0f, 3.0f };
    state->camera.target = (Vector3){ 0.0f, 2.0f, 0.0f };

    DebugSample initial = debug_collect_sample("initial-fluid-dynamic-floor", 0, 0.0f);
    state->initial_fluid_centroid_y = initial.fluid_centroid_y;
    if (container_created != state->expected_container_voxels) {
        debug_add_failure(state, "pbf-dynamic-container created %d/%d container voxels",
                          container_created, state->expected_container_voxels);
    }
    if (fluid_cells_created != 64 ||
        fluid_particle_count != state->expected_fluid_particles) {
        debug_add_failure(state,
                          "pbf-dynamic-container created cells=%d/64 particles=%d/%d",
                          fluid_cells_created, fluid_particle_count,
                          state->expected_fluid_particles);
    }
    return state->failure_count == 0;
}

static bool debug_setup_arena_platform_fluid(DebugRunState *state)
{
    debug_reset_world();
    activePlayers = 1;
    buildTestWorld();

    int M = (int)(2.0f * FLOOR_SIZE / VOXEL_SIZE);
    int platform_size = 10;
    int platform_base_height = 5;
    int platform_height = 1;
    int platform_min = M/2 - platform_size/2;
    int platform_max = M/2 + platform_size/2;
    int platform_lip_y = platform_base_height + platform_height + 1;
    Color platform_lip_color = (Color){ 70, 165, 105, 255 };
    Color platform_fluid_color = (Color){ 45, 145, 235, 220 };
    for (int x = platform_min; x <= platform_max; ++x) {
        for (int z = platform_min; z <= platform_max; ++z) {
            if (x != platform_min && x != platform_max &&
                z != platform_min && z != platform_max) {
                continue;
            }
            float px = (x + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
            float py = (platform_lip_y + 0.5f) * VOXEL_SIZE;
            float pz = (z + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
            addVoxel(px, py, pz, true, false, platform_lip_color, VOXEL_TYPE_SOLID);
        }
    }
    for (int x = platform_min + 2; x <= platform_max - 2; ++x) {
        for (int z = platform_min + 2; z <= platform_max - 2; ++z) {
            float px = (x + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
            float py = (platform_lip_y + 0.5f) * VOXEL_SIZE;
            float pz = (z + 0.5f) * VOXEL_SIZE - FLOOR_SIZE;
            add_fluid_cell(px, py, pz, platform_fluid_color);
        }
    }

    ground_structural_pillars();
    heal_structural_seams();
    rebuild_glue_constraints();

    state->arena_platform_fluid_test = true;
    state->activation_succeeded = true;
    state->camera.position = (Vector3){ 0.0f, 15.0f, 15.0f };
    state->camera.target = (Vector3){ 0.0f, 6.0f, 0.0f };

    int fluid_cells = 0;
    for (int i = 0; i < voxel_count; ++i) {
        if (voxel_is_fluid(&voxels[i])) {
            fluid_cells++;
            deactivate_fluid_cell(i);
        }
    }
    fprintf(stderr, "DEBUG REPRO: before platform activation: fluid_cells=%d active_fluid_particles=%d total_fluid_particles=%d\n",
            fluid_cells, active_fluid_particle_count, fluid_particle_count);

    int plat_gx = 0;
    int plat_gy = platform_lip_y - 1;
    int plat_gz = 0;
    int plat_idx = table_get_static_only(plat_gx, plat_gy, plat_gz);
    fprintf(stderr, "DEBUG REPRO: plat_idx at (%d,%d,%d) = %d\n", plat_gx, plat_gy, plat_gz, plat_idx);

    bool activated = activate_static_voxel_for_tether(plat_idx, 0, 1.0f);
    if (!activated) {
        debug_add_failure(state, "Failed to activate platform voxel at (%d,%d,%d)", plat_gx, plat_gy, plat_gz);
    }

    int fluid_cells_after = 0;
    for (int i = 0; i < voxel_count; ++i) {
        if (voxel_is_fluid(&voxels[i])) {
            fluid_cells_after++;
        }
    }

    fprintf(stderr, "DEBUG REPRO: after platform activation: fluid_cells=%d active_fluid_particles=%d total_fluid_particles=%d\n",
            fluid_cells_after, active_fluid_particle_count, fluid_particle_count);

    if (fluid_cells_after != 49) {
        debug_add_failure(state, "Fluid voxels disappeared upon platform activation: expected 49, found %d", fluid_cells_after);
    }
    if (active_fluid_particle_count != 392) {
        debug_add_failure(state, "Active fluid particles disappeared upon platform activation: expected 392, found %d", active_fluid_particle_count);
    }
    if (fluid_particle_count != 392) {
        debug_add_failure(state, "Total fluid particles disappeared upon platform activation: expected 392, found %d", fluid_particle_count);
    }
    return state->failure_count == 0;
}

static bool debug_setup_pbf_scaling(DebugRunState *state)
{
    debug_reset_world();
    activePlayers = 0;
    fluidRenderMode = FLUID_RENDER_SURFACE;

    const Color fluid_color = { 45, 145, 235, 210 };
    const int side = debugRunConfig.pbf_side;
    const float spacing = VOXEL_SIZE;
    const float origin = -0.5f * (float)(side - 1) * spacing;
    int fluid_cells_created = 0;

    for (int y = 0; y < side; ++y) {
        for (int z = 0; z < side; ++z) {
            for (int x = 0; x < side; ++x) {
                int index = add_fluid_cell(origin + (float)x * spacing,
                                           100.0f + (float)y * spacing,
                                           origin + (float)z * spacing,
                                           fluid_color);
                if (index >= 0) {
                    voxels[index].debugClusterTag = DEBUG_FLUID_TAG;
                    ++fluid_cells_created;
                }
            }
        }
    }

    debug_rebuild_world_state();
    state->pbf_scaling_benchmark = true;
    state->activation_succeeded = true;
    state->expected_fluid_particles = side * side * side * VOXEL_CORNER_COUNT;
    state->benchmark_warmup_steps = debugRunConfig.steps / 4;
    if (state->benchmark_warmup_steps < 10) state->benchmark_warmup_steps = 10;
    if (state->benchmark_warmup_steps >= debugRunConfig.steps) {
        state->benchmark_warmup_steps = debugRunConfig.steps - 1;
    }
    state->camera.position = (Vector3){ (float)side * 0.8f, 102.0f + (float)side * 0.5f,
                                       (float)side * 1.2f };
    state->camera.target = (Vector3){ 0.0f, 100.0f + (float)side * 0.25f, 0.0f };

    if (fluid_cells_created != side * side * side ||
        fluid_particle_count != state->expected_fluid_particles) {
        debug_add_failure(state,
                          "pbf-scaling created cells=%d/%d particles=%d/%d",
                          fluid_cells_created, side * side * side,
                          fluid_particle_count, state->expected_fluid_particles);
    }
    return state->failure_count == 0;
}

static bool debug_setup_pbf_dense_neighbors(DebugRunState *state)
{
    int requested_side = debugRunConfig.pbf_side;
    debugRunConfig.pbf_side = 4;
    bool ok = debug_setup_pbf_scaling(state);
    debugRunConfig.pbf_side = requested_side;
    if (!ok) return false;

    const int row = 8;
    const float spacing = PBF_SMOOTHING_RADIUS / 20.0f;
    const Vector3 center = { 0.0f, 100.0f, 0.0f };
    for (int i = 0; i < fluid_particle_count; ++i) {
        int x = i % row;
        int y = (i / row) % row;
        int z = i / (row * row);
        Vector3 position = {
            center.x + ((float)x - 3.5f) * spacing,
            center.y + ((float)y - 3.5f) * spacing,
            center.z + ((float)z - 3.5f) * spacing
        };
        Particle *particle = fluid_particles[i];
        particle->pos = position;
        particle->prev_pos = position;
        particle->predicted_pos = position;
        particle->vel = (Vector3){ 0.0f, 0.0f, 0.0f };
    }
    state->expect_gpu_neighbor_overflow = true;
    state->camera.position = (Vector3){ 2.5f, 102.0f, 3.0f };
    state->camera.target = center;
    return true;
}

static void debug_render_benchmark_frame(DebugRunState *state)
{
    BeginTextureMode(state->target);
        ClearBackground((Color){ 32, 42, 58, 255 });
        BeginMode3D(state->camera);
            DrawVoxels(state->camera);
        EndMode3D();
    EndTextureMode();
    rlDrawRenderBatchActive();
}

static void debug_init_enemy_player(int index, Vector3 pos, EnemyType type, InputType input)
{
    memset(&players[index], 0, sizeof(players[index]));
    memset(&botStates[index], 0, sizeof(botStates[index]));
    memset(&goliathStates[index], 0, sizeof(goliathStates[index]));
    players[index].pos = pos;
    players[index].death_pos = pos;
    players[index].yaw = 0.0f;
    players[index].pitch = 0.0f;
    players[index].onGround = true;
    players[index].respawn_timer = 0.0f;
    players[index].enemyType = type;
    players[index].matterMax = (type == ENEMY_TYPE_SWARMER) ? (MATTER_MAX_DEFAULT * 0.25f) : MATTER_MAX_DEFAULT;
    players[index].matter = players[index].matterMax;
    players[index].last_damage_time = -1000.0f;
    playerInput[index] = input;
    botStates[index].targetIndex = -1;
}

static bool debug_setup_swarmer_wall(DebugRunState *state)
{
    debug_reset_world();
    state->swarmer_wall_test = true;
    state->enemy_bot_index = 0;
    state->enemy_human_index = 1;
    state->swarmer_inside_wall_frames = 0;
    state->swarmer_min_z = 1000.0f;
    state->activation_succeeded = true;
    state->drop_sampled = true;
    state->expected_voxels = 0;
    gameMode = GAME_MODE_FIREFIGHT;
    activePlayers = 2;

    Color wall_color = (Color){ 110, 120, 140, 255 };
    int wall_created = 0;
    for (int y = 0; y <= 4; ++y) {
        for (int x = -2; x <= 2; ++x) {
            if (debug_add_cell(x, y, 0, true, false, wall_color, DEBUG_SUPPORT_TAG) >= 0) {
                ++wall_created;
            }
        }
    }
    debug_rebuild_world_state();
    debug_init_enemy_player(0, (Vector3){ 0.0f, BASE_EYE_HEIGHT, 3.0f },
                            ENEMY_TYPE_SWARMER, INPUT_TYPE_BOT_MEDIUM);
    debug_init_enemy_player(1, (Vector3){ 0.0f, BASE_EYE_HEIGHT, -3.0f },
                            ENEMY_TYPE_STANDARD, INPUT_TYPE_KEYBOARD);
    players[0].yaw = 0.0f;
    players[1].yaw = 180.0f;
    state->camera.position = (Vector3){ 6.0f, 4.5f, 6.0f };
    state->camera.target = (Vector3){ 0.0f, 1.2f, 0.0f };
    if (wall_created < 20) {
        debug_add_failure(state, "swarmer-wall created only %d wall voxels", wall_created);
        return false;
    }
    return true;
}

static bool debug_setup_goliath_flight(DebugRunState *state)
{
    debug_reset_world();
    state->goliath_flight_test = true;
    state->enemy_bot_index = 0;
    state->enemy_human_index = 1;
    state->activation_succeeded = true;
    state->drop_sampled = true;
    state->expected_voxels = 0;
    gameMode = GAME_MODE_FIREFIGHT;
    activePlayers = 2;
    debug_init_enemy_player(0, (Vector3){ 0.0f, BASE_EYE_HEIGHT, 2.0f },
                            ENEMY_TYPE_GOLIATH, INPUT_TYPE_BOT_MEDIUM);
    debug_init_enemy_player(1, (Vector3){ 0.0f, BASE_EYE_HEIGHT, -4.0f },
                            ENEMY_TYPE_STANDARD, INPUT_TYPE_KEYBOARD);
    spawn_goliath_armor(0);
    debug_rebuild_world_state();
    state->goliath_start_y = players[0].pos.y;
    state->goliath_max_y = players[0].pos.y;
    state->goliath_y_after_landing = players[0].pos.y;
    state->goliath_hover_frames = 0;
    players[0].yaw = 0.0f;
    state->camera.position = (Vector3){ 5.5f, 4.0f, 6.0f };
    state->camera.target = (Vector3){ 0.0f, 1.5f, 0.0f };
    if (goliathStates[0].count <= 0) {
        debug_add_failure(state, "goliath-flight spawned no armor voxels");
        return false;
    }
    return true;
}

static bool debug_setup_goliath_behavior(DebugRunState *state)
{
    debug_reset_world();
    state->goliath_behavior_test = true;
    state->enemy_bot_index = 0;
    state->enemy_human_index = 1;
    state->activation_succeeded = true;
    state->drop_sampled = true;
    state->expected_voxels = 0;
    gameMode = GAME_MODE_FIREFIGHT;
    activePlayers = 2;
    debug_init_enemy_player(0, (Vector3){ 0.0f, BASE_EYE_HEIGHT, 5.0f },
                            ENEMY_TYPE_GOLIATH, INPUT_TYPE_BOT_MEDIUM);
    debug_init_enemy_player(1, (Vector3){ 0.0f, BASE_EYE_HEIGHT, -5.0f },
                            ENEMY_TYPE_STANDARD, INPUT_TYPE_KEYBOARD);
    spawn_goliath_armor(0);
    goliathStates[0].launchTimer = 0.5f;
    players[0].yaw = 0.0f;
    players[1].yaw = 180.0f;
    debug_rebuild_world_state();
    state->camera.position = (Vector3){ 13.0f, 7.0f, 4.0f };
    state->camera.target = (Vector3){ 0.0f, 1.0f, -1.0f };
    if (goliathStates[0].count <= 0) {
        debug_add_failure(state, "goliath-behavior spawned no armor voxels");
        return false;
    }
    state->goliath_min_armor_count = goliathStates[0].count;
    state->goliath_final_armor_count = goliathStates[0].count;
    state->goliath_visible_frames = 0;
    state->goliath_voxels_launched = 0;
    return true;
}

static void debug_physics_step(DebugRunState *state)
{
    float dt = PBD_MAX_STEP_DT;
    if (state && state->coarsening_demo) {
        ++state->refinement_elapsed;
        if (state->refinement_elapsed == state->refinement_step) {
            if (!debug_coarsen_grid(state))
                debug_add_failure(state, "coarsening transfer failed at step %d",
                                  state->refinement_elapsed);
            else {
                state->coarsening_complete = true;
                debugRefinementCoarse = true;
                collisionTopologyDirty = true;
            }
        }
    }
    if (state && state->refinement_demo) {
        ++state->refinement_elapsed;
        if (state->refinement_elapsed == state->refinement_step) {
            if (!debug_sync_refinement_grid(state, true))
                debug_add_failure(state, "refinement transfer failed at step %d",
                                  state->refinement_elapsed);
            else {
                state->refinement_complete = true;
                debugRefinementCoarse = false;
                collisionTopologyDirty = true;
            }
        }
    }
    if (state && (state->swarmer_wall_test || state->goliath_flight_test || state->goliath_behavior_test)) {
        int bot = state->enemy_bot_index;
        int human = state->enemy_human_index;
        if (state->goliath_flight_test && state->goliath_forced_jumps < 1 &&
            players[bot].onGround) {
            players[bot].vel.y = JUMP_SPEED;
            players[bot].onGround = false;
            state->goliath_forced_jumps++;
        }
        if (state->swarmer_wall_test) {
            /* Same straight-line charge UpdateBot uses once a swarmer has a target. */
            Vector3 to_target = v_sub(players[human].pos, players[bot].pos);
            to_target.y = 0.0f;
            float len = v_length(to_target);
            if (len > 0.05f) {
                Vector3 moveDir = v_mul(to_target, 1.0f / len);
                float swarmerSpeed = MOVE_SPEED * 1.35f;
                players[bot].vel.x = moveDir.x * swarmerSpeed;
                players[bot].vel.z = moveDir.z * swarmerSpeed;
            }
        } else if (state->goliath_flight_test) {
            /* Walk toward the target like UpdateBot, but only the scripted jump. */
            Vector3 to_target = v_sub(players[human].pos, players[bot].pos);
            to_target.y = 0.0f;
            float len = v_length(to_target);
            if (len > 0.8f) {
                Vector3 moveDir = v_mul(to_target, 1.0f / len);
                float goliathSpeed = MOVE_SPEED * 0.85f;
                players[bot].vel.x = moveDir.x * goliathSpeed;
                players[bot].vel.z = moveDir.z * goliathSpeed;
            } else {
                players[bot].vel.x = 0.0f;
                players[bot].vel.z = 0.0f;
            }
        } else if (state->goliath_behavior_test) {
            int count_before = goliathStates[bot].count;
            UpdateBot(bot, dt);
            if (goliathStates[bot].count < count_before) {
                state->goliath_voxels_launched += (count_before - goliathStates[bot].count);
            }
            float dist = 0.0f;
            if (check_ai_awareness_of_player(bot, human, &dist)) {
                state->goliath_visible_frames++;
            }
            if (goliathStates[bot].count < state->goliath_min_armor_count) {
                state->goliath_min_armor_count = goliathStates[bot].count;
            }
            state->goliath_final_armor_count = goliathStates[bot].count;
        } else {
            UpdateBot(bot, dt);
        }
        net_integrate_player(&players[bot], dt, true);
        net_integrate_player(&players[human], dt, true);
        if (players[bot].enemyType == ENEMY_TYPE_GOLIATH) {
            update_goliath_armor_positions(bot);
            rebuild_voxel_hash();
        }
    }
    if (!state || (!state->active_cube_drop && !state->refinement_demo && !state->coarsening_demo)) {
        activate_static_voxels_near_dynamic();
        prepare_tether_forces();
        wake_sleeping_voxels_near_awake();
    }
    rebuild_static_hash_if_dirty();
    if (state && state->adaptive_cube_drop) {
        // Commit each GPU substep before sampling curvature. A topology change
        // invalidates the resident constraint list for the next substep.
        const float sub_dt = PBD_MAX_STEP_DT / (float)PBD_SUBSTEPS;
        for (int substep = 0; substep < PBD_SUBSTEPS; ++substep) {
            if (state->rip_test) {
                ++state->rip_substeps;
                debug_apply_rip_grips(state);
            }
            simulate_voxel_pbd_steps_counted(sub_dt, 1, 1);
            vgs_hierarchy_sample_curvature(sub_dt);
            // Curvature adapts after every substep. Strain/shear participates
            // only in the final AMR decision of each fixed physics frame.
            if (!vgs_hierarchy_adapt(substep == PBD_SUBSTEPS - 1)) {
                debug_add_failure(state, "adaptive octree transition failed");
                break;
            }
            if (state->rip_test) debug_apply_rip_grips(state);
        }
    } else {
        simulate_voxel_pbd(PBD_MAX_STEP_DT);
    }
    if (state && ((state->refinement_demo && !state->refinement_complete) ||
                  (state->coarsening_demo && state->coarsening_complete)) &&
        !debug_sync_refinement_grid(state, false))
        debug_add_failure(state, "coarse grid synchronization failed at step %d",
                          state->refinement_elapsed);
    if (state && state->active_cube_drop) {
        debug_sample_vgs_temporal_curvature(PBD_MAX_STEP_DT);
        int curvature_frame = state->adaptive_cube_drop
            ? (vgsHierarchy.sample_count - 1) / PBD_SUBSTEPS
            : vgsHierarchy.sample_count - 1;
        debug_write_curvature_levels(state, curvature_frame);
    }
    if (!state || (!state->active_cube_drop && !state->refinement_demo && !state->coarsening_demo)) recycle_dead_voxels();
    rebuild_voxel_hash();
    if (!state || (!state->active_cube_drop && !state->refinement_demo && !state->coarsening_demo)) deactivate_sleeping_voxels();
    rebuild_static_hash_if_dirty();
}

static void debug_record_rip_step(DebugRunState *state, int step)
{
    if (!state->rip_test || !state->rip_history_file) return;
    float half_extension = 0.5f * (state->rip_right_x - state->rip_left_x) *
                           state->rip_commanded_strain;
    float left_target = state->rip_left_x - half_extension;
    float right_target = state->rip_right_x + half_extension;
    float rest_step = 0.5f * pbdSolidVoxelSize;
    int left_count = 0, right_count = 0;
    int left_simulated = 0, right_simulated = 0;
    float max_error = 0.0f;
    for (int i = 0; i < particle_pool_count; ++i) {
        Particle *particle = &particles_pool[i];
        if (!particle->active) continue;
        float rest_x = (float)particle->rest_cell_x * rest_step;
        bool left = rest_x <= state->rip_left_x + state->rip_grip_depth + 1e-5f;
        bool right = rest_x >= state->rip_right_x - state->rip_grip_depth - 1e-5f;
        if (!left && !right) continue;
        if (left) ++left_count;
        if (right) ++right_count;
        if (particle->base_inv_mass > 0.0f) {
            if (left) ++left_simulated;
            if (right) ++right_simulated;
            float target_x = rest_x + (left ? -half_extension : half_extension);
            float error = fabsf(particle->pos.x - target_x);
            if (error > max_error) max_error = error;
        }
    }
    int broken_nodes = 0, broken_leaves = 0, enabled_vgs = 0, enabled_welds = 0;
    float max_strain = 0.0f;
    for (int i = 0; i < vgsHierarchy.node_count; ++i) {
        VgsHierarchyNode *node = &vgsHierarchy.nodes[i];
        if (!node->vgs_enabled) ++broken_nodes;
        if (!node->active) continue;
        if (node->vgs_enabled) ++enabled_vgs;
        if (vgs_terminal_child_mask(node)) ++enabled_welds;
        Vector3 corners[8];
        for (int c = 0; c < 8; ++c) corners[c] = node->particles[c]->pos;
        VgsDeformation d = measure_vgs_deformation(corners, node->rest_edge);
        if (d.max_abs_strain > max_strain) max_strain = d.max_abs_strain;
    }
    for (int i = 0; i < voxel_count; ++i) {
        Voxel *voxel = &voxels[i];
        if (!voxel->simulate || voxel->debugClusterTag != DEBUG_TEST_TAG) continue;
        if (!voxel->vgs_active) ++broken_leaves;
        if (!voxel->simulate_dofs) continue;
        if (voxel->vgs_active) ++enabled_vgs;
        Vector3 corners[8];
        for (int c = 0; c < 8; ++c) corners[c] = voxel->particles[c]->pos;
        VgsDeformation d = measure_vgs_deformation(corners, voxel->rest_edge);
        if (d.max_abs_strain > max_strain) max_strain = d.max_abs_strain;
    }
    state->rip_broken_nodes = broken_nodes;
    state->rip_broken_leaves = broken_leaves;
    if (step > 0 && state->rip_first_break_step < 0 && broken_nodes + broken_leaves > 0)
        state->rip_first_break_step = step;
    fprintf(state->rip_history_file,
            "%d,%.9g,%.9g,%.9g,%d,%d,%d,%d,%.9g,%.9g,%d,%d,%d,%d,%u,%u\n",
            step, (double)state->rip_commanded_strain,
            (double)left_target, (double)right_target, left_count, right_count,
            left_simulated, right_simulated,
            (double)max_error, (double)max_strain, broken_nodes, broken_leaves,
            enabled_vgs, enabled_welds,
            vgsHierarchy.refine_count, vgsHierarchy.coarsen_count);
    if (left_simulated < 4 || right_simulated < 4 || !isfinite(max_error) ||
        max_error > 1e-4f)
        debug_add_failure(state, "rip step %d simulated grip count/error %d/%d %.6g",
                          step, left_simulated, right_simulated, (double)max_error);
    if (step == debugRunConfig.steps / 2 &&
        fabsf(state->rip_commanded_strain - STRAIN_BREAK_THRESHOLD) > 1e-5f)
        debug_add_failure(state, "rip midpoint strain %.6g differs from threshold %.6g",
                          (double)state->rip_commanded_strain,
                          (double)STRAIN_BREAK_THRESHOLD);
}

static float debug_vgs_max_abs_strain(Particle *const particles[8], float rest_edge)
{
    Vector3 corners[8];
    for (int c = 0; c < 8; ++c) corners[c] = particles[c]->pos;
    return measure_vgs_deformation(corners, rest_edge).max_abs_strain;
}

static void debug_record_rip_parent_child_strain(DebugRunState *state, int step)
{
    if (!state->rip_test || !state->rip_strain_levels_file) return;
    for (int level = 0; level < vgsHierarchy.levels; ++level) {
        int count = 0, ratio_count = 0, intact_count = 0;
        double parent_sum = 0.0, child_sum = 0.0, ratio_sum = 0.0;
        double intact_parent_sum = 0.0, intact_child_sum = 0.0;
        int first = vgsHierarchy.level_start[level];
        int end = first + vgsHierarchy.level_count[level];
        for (int i = first; i < end; ++i) {
            VgsHierarchyNode *node = &vgsHierarchy.nodes[i];
            if (!node->active || !vgs_node_has_active_children(node)) continue;
            float parent_strain = debug_vgs_max_abs_strain(node->particles,
                                                            node->rest_edge);
            float child_mean = 0.0f;
            bool intact = node->vgs_enabled;
            for (int c = 0; c < 8; ++c) {
                int child = node->children[c];
                if (child >= 0) {
                    VgsHierarchyNode *child_node = &vgsHierarchy.nodes[child];
                    child_mean += debug_vgs_max_abs_strain(child_node->particles,
                                                            child_node->rest_edge);
                    intact = intact && child_node->vgs_enabled;
                } else {
                    Voxel *child_voxel = &voxels[-child - 1];
                    child_mean += debug_vgs_max_abs_strain(child_voxel->particles,
                                                            child_voxel->rest_edge);
                    intact = intact && child_voxel->vgs_active;
                }
            }
            child_mean *= 0.125f;
            ++count;
            parent_sum += parent_strain;
            child_sum += child_mean;
            if (child_mean > 1e-8f) {
                ratio_sum += parent_strain / child_mean;
                ++ratio_count;
            }
            if (intact) {
                ++intact_count;
                intact_parent_sum += parent_strain;
                intact_child_sum += child_mean;
            }
        }
        if (!count) continue;
        double parent_mean = parent_sum / count;
        double child_mean = child_sum / count;
        double intact_parent_mean = intact_count ? intact_parent_sum / intact_count : 0.0;
        double intact_child_mean = intact_count ? intact_child_sum / intact_count : 0.0;
        fprintf(state->rip_strain_levels_file,
                "%d,%d,%d,%.9g,%.9g,%.9g,%.9g,%d,%.9g,%.9g,%.9g\n",
                step, level, count, parent_mean, child_mean,
                child_mean > 1e-12 ? parent_mean / child_mean : 0.0,
                ratio_count ? ratio_sum / ratio_count : 0.0,
                intact_count, intact_parent_mean, intact_child_mean,
                intact_child_mean > 1e-12 ? intact_parent_mean / intact_child_mean : 0.0);
    }
}

static uint64_t debug_tagged_static_identity_signature(int tag)
{
    uint64_t signature = 0;
    for (int i = 0; i < voxel_count; ++i) {
        const Voxel *voxel = &voxels[i];
        if (voxel->simulate || voxel->debugClusterTag != tag) continue;
        uint64_t value = voxel->identity;
        value ^= (uint64_t)(uint32_t)voxel->gx * UINT64_C(0x9e3779b185ebca87);
        value ^= (uint64_t)(uint32_t)voxel->gy * UINT64_C(0xc2b2ae3d27d4eb4f);
        value ^= (uint64_t)(uint32_t)voxel->gz * UINT64_C(0x165667b19e3779f9);
        value ^= value >> 30;
        value *= UINT64_C(0xbf58476d1ce4e5b9);
        value ^= value >> 27;
        signature += value;
    }
    return signature;
}

static void debug_validate_step(DebugRunState *state, int step)
{
    debug_validate_all_sim_particles(state, step);
    if (state->swarmer_wall_test) {
        Player *bot = &players[state->enemy_bot_index];
        if (bot->pos.z < state->swarmer_min_z) state->swarmer_min_z = bot->pos.z;
        float half = player_collision_half(bot);
        if (player_aabb_hits_world(bot->pos, half)) {
            state->swarmer_inside_wall_frames++;
        }
    }
    if (state->goliath_flight_test) {
        Player *bot = &players[state->enemy_bot_index];
        if (bot->pos.y > state->goliath_max_y) state->goliath_max_y = bot->pos.y;
        if (step == debugRunConfig.steps) state->goliath_y_after_landing = bot->pos.y;
        if (bot->onGround && bot->pos.y > state->goliath_start_y + 0.45f) {
            state->goliath_hover_frames++;
        }
    }
    // Isolate the post-impact lifecycle under test: wait for the real swept
    // tether collision, then wake the struck pillar through the normal static
    // activation path.  Ground support otherwise makes this compact fixture's
    // freeze belief much higher than a damaged gameplay structure's.
    if (state->sleep_snap_original_test &&
        !state->sleep_snap_impact_activation_requested &&
        tetherThrowCcdHitCount > state->tether_thin_wall_hits_before) {
        int seed = debug_find_tagged_voxel(DEBUG_TEST_TAG, true);
        if (seed >= 0) {
            state->sleep_snap_impact_activation_requested =
                activate_static_voxel_for_tether(seed, 0, ACTIVATION_TETHER_BELIEF);
        }
    }
    DebugSample sample = debug_collect_sample("validation", step, state->initial_centroid_y);
    if (state->refinement_demo &&
        (step == 0 || step == state->refinement_step - 1 ||
         step == state->refinement_step || step == debugRunConfig.steps)) {
        bool expect_refined = step >= state->refinement_step;
        if (step == state->refinement_step && state->gpu_readback_performed) {
            // Recheck the uploaded particle masses after enabling the children.
            state->gpu_readback_performed = false;
            state->gpu_zero_mass_corners = 0;
        }
        int mobile_points = 0, enabled_children = 0;
        double total_mass = 0.0;
        for (int r = 0; r < 27; ++r) {
            Particle *particle = state->refine_grid[r];
            if (!particle || !particle->active) {
                debug_add_failure(state, "step %d: refinement grid point %d is invalid", step, r);
                continue;
            }
            if (particle->inv_mass > 0.0f) {
                ++mobile_points;
                total_mass += 1.0 / particle->inv_mass;
            }
        }
        for (int c = 0; c < 8; ++c)
            if (voxels[state->refine_children[c]].simulate_dofs) ++enabled_children;
        if (state->refinement_complete != expect_refined ||
            mobile_points != (expect_refined ? 27 : 8) ||
            enabled_children != (expect_refined ? 8 : 0) ||
            fabs(total_mass - 8000.0) > 0.1) {
            debug_add_failure(state,
                "step %d: refinement state mismatch (done=%d mobile=%d children=%d mass=%.1f)",
                step, state->refinement_complete, mobile_points, enabled_children, total_mass);
        }
        if (sample.invalid_particle_refs > 0)
            debug_add_failure(state, "step %d: %d invalid refinement particle references",
                              step, sample.invalid_particle_refs);
        if (expect_refined) debug_validate_tagged_dynamic(state, &sample, "after refinement");
    }
    if (state->coarsening_demo &&
        (step == 0 || step == state->refinement_step - 1 ||
         step == state->refinement_step || step == debugRunConfig.steps)) {
        bool expect_coarse = step >= state->refinement_step;
        if (step == state->refinement_step && state->gpu_readback_performed) {
            state->gpu_readback_performed = false;
            state->gpu_zero_mass_corners = 0;
        }
        int mobile_points = 0, enabled_children = 0;
        double total_mass = 0.0;
        for (int r = 0; r < 27; ++r) {
            Particle *particle = state->refine_grid[r];
            if (!particle || !particle->active) {
                debug_add_failure(state, "step %d: coarsening grid point %d is invalid", step, r);
                continue;
            }
            if (particle->inv_mass > 0.0f) {
                ++mobile_points;
                total_mass += 1.0 / particle->inv_mass;
            }
        }
        for (int c = 0; c < 8; ++c)
            if (voxels[state->refine_children[c]].simulate_dofs) ++enabled_children;
        if (state->coarsening_complete != expect_coarse ||
            vgsHierarchy.node_count != 1 ||
            vgsHierarchy.nodes[0].active != expect_coarse ||
            mobile_points != (expect_coarse ? 8 : 27) ||
            enabled_children != (expect_coarse ? 0 : 8) ||
            fabs(total_mass - 8000.0) > 0.1 ||
            !isfinite(state->coarsening_projection_rms)) {
            debug_add_failure(state,
                "step %d: coarsening state mismatch (done=%d mobile=%d children=%d mass=%.1f rms=%.4f)",
                step, state->coarsening_complete, mobile_points, enabled_children,
                total_mass, state->coarsening_projection_rms);
        }
        if (sample.invalid_particle_refs > 0)
            debug_add_failure(state, "step %d: %d invalid coarsening particle references",
                              step, sample.invalid_particle_refs);
        debug_validate_tagged_dynamic(state, &sample,
                                      expect_coarse ? "after coarsening" : "before coarsening");
    }
    if (state->active_cube_drop && !state->adaptive_cube_drop &&
        (step == 60 || step == debugRunConfig.steps)) {
        debug_validate_tagged_dynamic(state, &sample, "cube drop");
        if (step == 60 && sample.tagged_rendered != state->expected_voxels) {
            debug_add_failure(state, "cube drop lost rendered PBD cells before impact (%d/%d)",
                              sample.tagged_rendered, state->expected_voxels);
        }
        if (PBD_DISABLE_VGS_DEACTIVATION && step == debugRunConfig.steps && step >= 120 &&
            sample.tagged_rendered != state->expected_voxels) {
            debug_add_failure(state, "cube drop retained only %d/%d rendered cells after impact",
                              sample.tagged_rendered, state->expected_voxels);
        }
    }
    if (state->adaptive_cube_drop) {
        if (getenv("FPS_AMR_TEST_BREAK_PARENT") && step == 1) {
            VgsHierarchyNode *root = &vgsHierarchy.nodes[0];
            int enabled_children = 0;
            for (int c = 0; c < 8; ++c)
                if (vgs_child_enabled(root->children[c])) ++enabled_children;
            if (root->vgs_enabled || enabled_children != 8 ||
                vgs_enabled_child_mask(root) != 0 ||
                vgs_terminal_child_mask(root) != 0)
                debug_add_failure(state,
                    "parent fracture enabled=%d children=%d childMask=0x%x terminalMask=0x%x",
                    root->vgs_enabled, enabled_children,
                    vgs_enabled_child_mask(root), vgs_terminal_child_mask(root));
            if (physics_backend_is_gpu(physicsBackend.active)) {
                uint32_t root_voxel_id = (uint32_t)gpuPhysics.dynamic_voxel_count;
                for (int i = 0; i < gpuPhysics.constraint_count; ++i)
                    if (gpuPhysics.constraint_stage[i].voxel_id == root_voxel_id)
                        debug_add_failure(state,
                            "disabled parent retained GPU constraint kind=%u",
                            gpuPhysics.constraint_stage[i].kind);
            }
        }
        if (getenv("FPS_AMR_TEST_BREAK_CHILD") && step == 1) {
            int broken = 0;
            for (int c = 0; c < 8; ++c)
                if (!vgs_child_enabled(vgsHierarchy.nodes[0].children[c])) ++broken;
            if (broken == 0 || !vgsHierarchy.nodes[0].vgs_enabled ||
                vgs_enabled_child_mask(&vgsHierarchy.nodes[0]) == 0xffu)
                debug_add_failure(state, "child fracture parent=%d broken=%d weldMask=0x%x",
                                  vgsHierarchy.nodes[0].vgs_enabled, broken,
                                  vgs_enabled_child_mask(&vgsHierarchy.nodes[0]));
        }
        if (getenv("FPS_AMR_TEST_BREAK_CHILD") && step > 0 &&
            !vgs_node_has_active_children(&vgsHierarchy.nodes[0]))
            debug_add_failure(state, "root coarsened across a broken child");
        int active_nodes = 0, terminal_nodes = 0, active_leaves = 0;
        double total_mass = 0.0;
        double terminal_volume = 0.0;
        for (int i = 0; i < sim_particle_count; ++i)
            if (sim_particles[i]->base_inv_mass > 0.0f)
                total_mass += 1.0 / sim_particles[i]->base_inv_mass;
        for (int i = 0; i < vgsHierarchy.node_count; ++i) {
            VgsHierarchyNode *node = &vgsHierarchy.nodes[i];
            for (int c = 0; c < 8; ++c) {
                int child = node->children[c];
                bool child_active = child >= 0 ? vgsHierarchy.nodes[child].active :
                    voxels[-child - 1].simulate_dofs;
                if (child_active && !node->active)
                    debug_add_failure(state, "step %d: active child beneath inactive node %d", step, i);
            }
            if (!node->active) continue;
            ++active_nodes;
            int children = 0;
            for (int c = 0; c < 8; ++c) {
                int child = node->children[c];
                children += child >= 0 ? vgsHierarchy.nodes[child].active :
                    voxels[-child - 1].simulate_dofs;
            }
            if (children != 0 && children != 8)
                debug_add_failure(state, "step %d: partial adaptive octet at node %d", step, i);
            if (children == 0) {
                if (node->vgs_enabled) ++terminal_nodes;
                terminal_volume += node->rest_volume;
            }
        }
        for (int i = 0; i < voxel_count; ++i)
            if (voxels[i].simulate && voxels[i].simulate_dofs) {
                if (voxels[i].vgs_active) ++active_leaves;
                terminal_volume += voxels[i].rest_volume;
            }
        if (!vgsHierarchy.nodes[0].active || active_nodes == 0 ||
            terminal_nodes + active_leaves != sample.tagged_rendered ||
            fabs(total_mass - 8000.0) > 0.1 ||
            fabs(terminal_volume - 125.0) > 0.01 || sample.invalid_particle_refs > 0)
            debug_add_failure(state,
                "step %d: adaptive state nodes=%d terminals=%d leaves=%d rendered=%d mass=%.3f volume=%.3f invalid=%d",
                step, active_nodes, terminal_nodes, active_leaves,
                sample.tagged_rendered, total_mass, terminal_volume, sample.invalid_particle_refs);
    }
    if (state->sleep_snap_original_test) {
        float max_error = 0.0f;
        float max_speed = 0.0f;
        float max_strain = 0.0f;
        float max_shear = 0.0f;
        for (int i = 0; i < voxel_count; ++i) {
            Voxel *voxel = &voxels[i];
            if (voxel->debugClusterTag != DEBUG_TEST_TAG || !voxel->simulate) continue;
            max_speed = fmaxf(max_speed, v_length(voxel->vel));
            float strain = 0.0f;
            float shear = 0.0f;
            voxel_measure_strain(voxel, &strain, &shear);
            max_strain = fmaxf(max_strain, strain);
            max_shear = fmaxf(max_shear, shear);
            Vector3 original_center = {
                ((float)voxel->orig_min_gx + (float)voxel->orig_max_gx + 1.0f) *
                    (VOXEL_SIZE * 0.5f),
                ((float)voxel->orig_min_gy + (float)voxel->orig_max_gy + 1.0f) *
                    (VOXEL_SIZE * 0.5f),
                ((float)voxel->orig_min_gz + (float)voxel->orig_max_gz + 1.0f) *
                    (VOXEL_SIZE * 0.5f)
            };
            Vector3 d = v_sub(voxel->pos, original_center);
            float err = v_length(d);
            max_error = fmaxf(max_error, err);
        }
        state->sleep_snap_max_original_error = max_error;
        state->sleep_snap_last_speed = max_speed;
        state->sleep_snap_last_strain = max_strain;
        state->sleep_snap_last_shear = max_shear;

        if (sample.tagged_dynamic == state->expected_voxels) {
            state->activation_succeeded = true;
            if (!state->sleep_snap_active_observed) {
                state->sleep_snap_active_observed = true;
                state->activation_step = step;
                debug_capture(state, "phase-02-impact-activated.png", "impact-activated", step);
            }
        }
        if (sample.tagged_sleeping > 0) {
            state->sleep_snap_sleeping_observed = true;
            state->sleeping_conversion_succeeded = true;
        }
        if (state->sleep_snap_active_observed &&
            sample.tagged_static == state->expected_voxels &&
            sample.tagged_dynamic == 0 && !state->sleep_snap_original_succeeded) {
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
            state->sleep_snap_original_succeeded = original_cells_static;
            state->sleeping_conversion_succeeded = original_cells_static;
            state->sleep_snap_static_step = original_cells_static ? step : -1;
            if (original_cells_static) {
                state->sleep_snap_queue_at_static = recycleQueueCount;
                state->sleep_snap_static_identity_signature =
                    debug_tagged_static_identity_signature(DEBUG_TEST_TAG);
            }
            debug_capture(state, "phase-03-snapped-static.png", "snapped-original-static", step);
        }
        if (state->sleep_snap_rebuild_test && state->sleep_snap_original_succeeded &&
            step > state->sleep_snap_static_step &&
            (sample.tagged_static != state->expected_voxels || sample.tagged_dynamic != 0 ||
             debug_tagged_static_identity_signature(DEBUG_TEST_TAG) !=
                state->sleep_snap_static_identity_signature)) {
            if (!state->sleep_snap_rebuilt_after_static) {
                state->sleep_snap_rebuilt_after_static = true;
                debug_capture(state, "phase-04-unnecessary-rebuild.png",
                              "unnecessary-static-rebuild", step);
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
    if (state->pbf_scaling_benchmark) {
        if (sample.fluid_particles != state->expected_fluid_particles) {
            debug_add_failure(state, "step %d: benchmark fluid count changed to %d (expected %d)",
                              step, sample.fluid_particles, state->expected_fluid_particles);
        }
        if (sample.invalid_fluid_particles > 0) {
            debug_add_failure(state, "step %d: %d benchmark fluid particles are invalid",
                              step, sample.invalid_fluid_particles);
        }
        return;
    }
    if (state->pbf_container) {
        if (!state->pbf_dynamic_container && sample.sim_particles != 0) {
            debug_add_failure(state,
                              "step %d: %d fluid particles leaked into the solid sim list",
                              step, sample.sim_particles);
        }
        if (sample.fluid_particles != state->expected_fluid_particles) {
            debug_add_failure(state, "step %d: fluid particle count changed to %d (expected %d)",
                              step, sample.fluid_particles, state->expected_fluid_particles);
        }
        if (sample.invalid_fluid_particles > 0) {
            debug_add_failure(state, "step %d: %d fluid particles are invalid or mis-indexed",
                              step, sample.invalid_fluid_particles);
        }
        int container_count = 0;
        int dynamic_container_count = 0;
        for (int i = 0; i < voxel_count; ++i) {
            const Voxel *voxel = &voxels[i];
            if (voxel->debugClusterTag != DEBUG_CONTAINER_TAG) continue;
            ++container_count;
            if (voxel->simulate) ++dynamic_container_count;
            if ((!state->pbf_dynamic_container && voxel->simulate) || voxel->pendingActivation) {
                debug_add_failure(state,
                                  "step %d: fluid activated an immutable container voxel",
                                  step);
                break;
            }
        }
        if (container_count != state->expected_container_voxels) {
            debug_add_failure(state, "step %d: container voxel count changed to %d (expected %d)",
                              step, container_count, state->expected_container_voxels);
        }
        if (state->pbf_dynamic_container && step <= 60 && dynamic_container_count < 49) {
            debug_add_failure(state,
                              "step %d: active floor lost dynamic voxels (%d/49 remain)",
                              step, dynamic_container_count);
        }
        if (step == 120) {
            state->fluid_drop_at_step_120 =
                state->initial_fluid_centroid_y - sample.fluid_centroid_y;
        }
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
               !state->sleep_snap_original_test && !state->arena_platform_fluid_test) {
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
            debug_add_failure(state,
                              "tether-activated pillar did not restore static "
                              "(sleeping=%s error=%.4f speed=%.4f strain=%.4f shear=%.4f)",
                              state->sleep_snap_sleeping_observed ? "true" : "false",
                              state->sleep_snap_max_original_error,
                              state->sleep_snap_last_speed,
                              state->sleep_snap_last_strain,
                              state->sleep_snap_last_shear);
        }
        if (state->sleep_snap_rebuild_test) {
            if (state->sleep_snap_queue_at_static != 0) {
                debug_add_failure(state,
                                  "static restoration left %d obsolete queued snapshots",
                                  state->sleep_snap_queue_at_static);
            }
            if (state->sleep_snap_rebuilt_after_static) {
                debug_add_failure(state,
                                  "successfully restored static pillar was rebuilt again");
            }
            if (recycleQueueCount != 0) {
                debug_add_failure(state,
                                  "restored pillar left %d stale recycle snapshots",
                                  recycleQueueCount);
            }
        }
        return;
    }
    if (state->swarmer_wall_test) {
        if (state->swarmer_inside_wall_frames > 3) {
            debug_add_failure(state, "swarmer spent %d frames inside the wall",
                              state->swarmer_inside_wall_frames);
        }
        if (state->swarmer_min_z < 0.20f) {
            debug_add_failure(state, "swarmer walked through the wall to z=%.3f",
                              state->swarmer_min_z);
        }
        return;
    }
    if (state->goliath_flight_test) {
        if (state->goliath_forced_jumps < 1) {
            debug_add_failure(state, "goliath never jumped");
        }
        if (state->goliath_hover_frames > 8) {
            debug_add_failure(state,
                              "goliath hovered on its armor for %d frames (max_y=%.3f final_y=%.3f)",
                              state->goliath_hover_frames, state->goliath_max_y,
                              state->goliath_y_after_landing);
        }
        if (state->goliath_y_after_landing > state->goliath_start_y + 0.55f) {
            Player *bot = &players[state->enemy_bot_index];
            float half = player_collision_half(bot);
            bool hits = player_aabb_hits_world(bot->pos, half);
            Vector3 probe = { bot->pos.x, bot->pos.y - 0.02f, bot->pos.z };
            bool grounded = player_aabb_hits_world(probe, half);
            debug_add_failure(state,
                              "goliath finished airborne at y=%.3f (start %.3f) onGround=%d vel.y=%.3f aabb=%d probe=%d armor=%d",
                              state->goliath_y_after_landing, state->goliath_start_y,
                              bot->onGround ? 1 : 0, bot->vel.y, hits ? 1 : 0, grounded ? 1 : 0,
                              goliathStates[state->enemy_bot_index].count);
        }
        return;
    }
    if (state->goliath_behavior_test) {
        fprintf(stderr,
                "[GoliathBehavior] visible_frames=%d launched=%d final_armor=%d min_armor=%d\n",
                state->goliath_visible_frames, state->goliath_voxels_launched,
                state->goliath_final_armor_count, state->goliath_min_armor_count);
        if (state->goliath_visible_frames < 20) {
            debug_add_failure(state,
                              "goliath was directly aware of visible player for only %d frames",
                              state->goliath_visible_frames);
        }
        if (state->goliath_voxels_launched < 1) {
            debug_add_failure(state,
                              "goliath never launched any voxels (launched=%d visible_frames=%d)",
                              state->goliath_voxels_launched, state->goliath_visible_frames);
        }
        if (state->goliath_final_armor_count <= 0) {
            debug_add_failure(state,
                              "goliath lost all armor voxels (final_armor=%d min_armor=%d)",
                              state->goliath_final_armor_count, state->goliath_min_armor_count);
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
        if (PBD_DISABLE_VGS_DEACTIVATION) {
            if (state->hinge_broken)
                debug_add_failure(state, "single-corner glue broke while VGS deactivation is disabled");
            return;
        }
        if (!state->hinge_broken) {
            debug_add_failure(state,
                              "single-corner glue survived %.1f degrees (limit %.1f)",
                              state->max_hinge_angle_deg, (double)GLUE_BREAK_HINGE_ANGLE_DEG);
        }
        return;
    }
    if (state->voxel_fracture_test) {
        int idx_strain = state->fracture_strain_voxel;
        int idx_shear = state->fracture_shear_voxel;
        int idx_ctrl = state->fracture_control_voxel;
        if (PBD_DISABLE_VGS_DEACTIVATION) {
            if ((idx_strain >= 0 && idx_strain < voxel_count && !voxels[idx_strain].vgs_active) ||
                (idx_shear >= 0 && idx_shear < voxel_count && !voxels[idx_shear].vgs_active) ||
                (idx_ctrl >= 0 && idx_ctrl < voxel_count && !voxels[idx_ctrl].vgs_active)) {
                debug_add_failure(state, "VGS deactivated while the temporary global switch is enabled");
            }
            return;
        }
        if (idx_strain >= 0 && idx_strain < voxel_count && voxels[idx_strain].vgs_active) {
            debug_add_failure(state, "high-strain voxel (strain > 0.20) did not break (vgs_active remains true)");
        }
        if (idx_shear >= 0 && idx_shear < voxel_count && voxels[idx_shear].vgs_active) {
            debug_add_failure(state, "high-shear voxel (shear > 0.20) did not break (vgs_active remains true)");
        }
        if (idx_ctrl >= 0 && idx_ctrl < voxel_count && !voxels[idx_ctrl].vgs_active) {
            debug_add_failure(state, "control voxel broke unexpectedly (vgs_active became false)");
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
    if (state->pbf_scaling_benchmark) {
        if (state->benchmark_measured_steps <= 0 ||
            !isfinite(state->benchmark_frame_total_ms) ||
            state->benchmark_frame_total_ms <= 0.0) {
            debug_add_failure(state, "benchmark did not record valid timed frames");
        }
        if (state->expect_gpu_neighbor_overflow &&
            physicsBackend.active == PHYSICS_BACKEND_GPU_GL43 &&
            state->benchmark_gpu_neighbor_overflow == 0) {
            debug_add_failure(state, "dense PBF case did not exercise the GPU neighbor overflow fallback");
        }
        return;
    }
    if (state->pbf_container) {
        DebugSample final_sample = debug_collect_sample("final", debugRunConfig.steps, 0.0f);
        state->escaped_fluid_particles = 0;
        state->below_container_particles = 0;
        for (int i = 0; i < fluid_particle_count; ++i) {
            const Particle *particle = fluid_particles[i];
            if (!particle || !particle->active || !v_isfinite(particle->pos)) continue;
            if (particle->pos.x < state->container_min_x ||
                particle->pos.x > state->container_max_x ||
                particle->pos.z < state->container_min_z ||
                particle->pos.z > state->container_max_z) {
                ++state->escaped_fluid_particles;
            }
            if (state->pbf_dynamic_container) {
                bool beneath_deformed_floor = false;
                bool inside_container_xz =
                    particle->pos.x >= state->container_min_x &&
                    particle->pos.x <= state->container_max_x &&
                    particle->pos.z >= state->container_min_z &&
                    particle->pos.z <= state->container_max_z;
                for (int voxel_index = 0; voxel_index < voxel_count; ++voxel_index) {
                    const Voxel *floor_voxel = &voxels[voxel_index];
                    if (floor_voxel->debugClusterTag != DEBUG_CONTAINER_TAG ||
                        floor_voxel->orig_min_gy != 4 ||
                        floor_voxel->orig_min_gx < -3 || floor_voxel->orig_max_gx > 3 ||
                        floor_voxel->orig_min_gz < -3 || floor_voxel->orig_max_gz > 3) continue;
                    Vector3 floor_min = { FLT_MAX, FLT_MAX, FLT_MAX };
                    Vector3 floor_max = { -FLT_MAX, -FLT_MAX, -FLT_MAX };
                    for (int corner = 0; corner < 8; ++corner) {
                        Vector3 p = floor_voxel->particles[corner]->pos;
                        if (p.x < floor_min.x) floor_min.x = p.x;
                        if (p.y < floor_min.y) floor_min.y = p.y;
                        if (p.z < floor_min.z) floor_min.z = p.z;
                        if (p.x > floor_max.x) floor_max.x = p.x;
                        if (p.y > floor_max.y) floor_max.y = p.y;
                        if (p.z > floor_max.z) floor_max.z = p.z;
                    }
                    if (inside_container_xz &&
                        particle->pos.x >= floor_min.x &&
                        particle->pos.x <= floor_max.x &&
                        particle->pos.z >= floor_min.z &&
                        particle->pos.z <= floor_max.z &&
                        particle->pos.y < floor_min.y - particle->radius - 0.01f) {
                        beneath_deformed_floor = true;
                        break;
                    }
                }
                if (beneath_deformed_floor) ++state->below_container_particles;
            } else if (particle->pos.y < state->container_floor_top - 0.01f) {
                ++state->below_container_particles;
            }
        }
        int maximum_escaped = state->pbf_dynamic_container
            ? (state->expected_fluid_particles + 9) / 10
            : (state->expected_fluid_particles + 19) / 20;
        float minimum_drop = state->pbf_dynamic_container ? 1.0f : 2.0f;
        if (!isfinite(state->fluid_drop_at_step_120) || state->fluid_drop_at_step_120 < minimum_drop) {
            debug_add_failure(state,
                              "fluid dropped only %.3f metres after 120 steps (minimum %.3f)",
                              state->fluid_drop_at_step_120, minimum_drop);
        }
        if (state->escaped_fluid_particles > maximum_escaped) {
            debug_add_failure(state, "%d fluid particles escaped the container (maximum %d)",
                              state->escaped_fluid_particles, maximum_escaped);
        }
        if (state->below_container_particles > 0) {
            debug_add_failure(state, "%d fluid particles penetrated below the container floor",
                              state->below_container_particles);
        }
        if (!isfinite(final_sample.mean_density_ratio) ||
            final_sample.mean_density_ratio < 0.70f ||
            final_sample.mean_density_ratio > 1.15f) {
            debug_add_failure(state,
                              "final mean density/rest ratio %.3f is outside 0.70..1.15",
                              final_sample.mean_density_ratio);
        }
        if (!isfinite(final_sample.max_density_ratio) ||
            final_sample.max_density_ratio > 1.35f) {
            debug_add_failure(state,
                              "final maximum density/rest ratio %.3f exceeds 1.35",
                              final_sample.max_density_ratio);
        }
        float maximum_mean_speed = state->pbf_dynamic_container ? 1.25f : 0.75f;
        if (!isfinite(final_sample.mean_fluid_speed) ||
            final_sample.mean_fluid_speed > maximum_mean_speed) {
            debug_add_failure(state,
                              "final mean fluid speed %.3f m/s exceeds %.2f",
                              final_sample.mean_fluid_speed, maximum_mean_speed);
        }
        return;
    }
    if (state->arena_platform_fluid_test) {
        int fluid_cells_final = 0;
        for (int i = 0; i < voxel_count; ++i) {
            if (voxel_is_fluid(&voxels[i])) fluid_cells_final++;
        }
        if (fluid_cells_final != 49) {
            debug_add_failure(state, "Fluid voxels lost during simulation: expected 49, found %d", fluid_cells_final);
        }
        if (fluid_particle_count != 392) {
            debug_add_failure(state, "Fluid particles lost during simulation: expected 392, found %d", fluid_particle_count);
        }
        return;
    }
    if (!state->activation_succeeded) {
        debug_add_failure(state, "scenario did not complete its activation stage");
    }
    if (state->refinement_demo && !state->refinement_complete) {
        debug_add_failure(state, "refinement demo never enabled its eight children");
    }
    if (state->coarsening_demo && !state->coarsening_complete) {
        debug_add_failure(state, "coarsening demo never projected its fine grid");
    }
    if (!state->drop_sampled) {
        debug_add_failure(state, "scenario did not reach the required 60-step drop sample");
    } else if (!state->rip_test && !state->overhang_impact && !state->platform_lifecycle_test &&
               state->expected_static_leftovers == 0 &&
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
    if (state->rip_test) {
        if (fabsf(state->rip_commanded_strain - 2.0f * STRAIN_BREAK_THRESHOLD) > 1e-5f)
            debug_add_failure(state, "rip final strain %.6g is not twice break threshold",
                              (double)state->rip_commanded_strain);
        if (state->rip_first_break_step < 0)
            debug_add_failure(state, "rip test did not fracture any VGS constraint");
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
    if (state->active_cube_drop) {
        int active_nodes = 0;
        int enabled_nodes = 0;
        for (int i = 0; i < vgsHierarchy.node_count; ++i)
            if (vgsHierarchy.nodes[i].active) {
                ++active_nodes;
                if (vgsHierarchy.nodes[i].vgs_enabled) ++enabled_nodes;
            }
        fprintf(file, "  \"pbdVoxelSize\": %.9g,\n  \"effectivePbdVoxelSize\": %.9g,\n  \"dynamicVoxels\": %d,\n  \"staticLeftovers\": %d,\n",
                (double)DEBUG_PBD_VOXEL_SIZE, (double)state->effective_pbd_size,
                state->expected_voxels,
                state->expected_static_leftovers);
        fprintf(file, "  \"vgsHierarchy\": {\"levels\": %d, \"internalNodes\": %d, \"activeInternalNodes\": %d, \"enabledInternalNodes\": %d, \"deactivationDisabled\": %s},\n",
                vgsHierarchy.levels + (vgsHierarchy.leaf_count ? 1 : 0),
                vgsHierarchy.node_count, active_nodes, enabled_nodes,
                PBD_DISABLE_VGS_DEACTIVATION ? "true" : "false");
        if (state->adaptive_cube_drop) {
            int enabled_vgs = enabled_nodes, enabled_welds = 0;
            int terminal_nodes = 0, terminal_leaves = 0, visible_terminal_cells = 0;
            int broken_nodes = 0, broken_leaves = 0;
            for (int i = 0; i < voxel_count; ++i)
                if (voxels[i].simulate) {
                    if (!voxels[i].vgs_active) ++broken_leaves;
                    if (voxels[i].simulate_dofs) {
                        ++terminal_leaves;
                        if (voxels[i].vgs_active) {
                            ++enabled_vgs;
                            ++visible_terminal_cells;
                        }
                    }
                }
            for (int i = 0; i < vgsHierarchy.node_count; ++i) {
                VgsHierarchyNode *node = &vgsHierarchy.nodes[i];
                if (!node->vgs_enabled) ++broken_nodes;
                if (!node->active) continue;
                if (vgs_terminal_child_mask(node)) ++enabled_welds;
                if (vgs_node_has_active_children(node)) continue;
                ++terminal_nodes;
                if (node->vgs_enabled) ++visible_terminal_cells;
            }
            fprintf(file, "  \"adaptiveHierarchy\": {\"thresholdPerSecondSquared\": %.9g, \"thresholdLevelFactor\": %d, \"deformationFractionCoarsest\": %.9g, \"deformationFractionFinest\": %.9g, \"cornerDrop\": %s, \"refinements\": %u, \"coarsenings\": %u, \"terminalCells\": %d, \"visibleTerminalCells\": %d, \"terminalNodes\": %d, \"terminalLeaves\": %d, \"brokenInternalNodes\": %d, \"brokenLeaves\": %d, \"enabledVgsConstraints\": %d, \"enabledWeldConstraints\": %d, \"gpuConstraints\": %d, \"simulatedParticles\": %d},\n",
                    (double)VGS_ADAPTIVE_CURVATURE_THRESHOLD,
                    VGS_ADAPTIVE_THRESHOLD_LEVEL_FACTOR,
                    (double)vgsAdaptiveDeformationCoarseFraction,
                    (double)vgsAdaptiveDeformationFineFraction,
                    getenv("FPS_AMR_CORNER_DROP") ? "true" : "false",
                    vgsHierarchy.refine_count, vgsHierarchy.coarsen_count,
                    terminal_nodes + terminal_leaves, visible_terminal_cells,
                    terminal_nodes, terminal_leaves,
                    broken_nodes, broken_leaves,
                    enabled_vgs, enabled_welds,
                    gpuPhysics.constraint_count, sim_particle_count);
        }
        fprintf(file, "  \"vgsTemporalCurvature\": {\"samples\": %d, \"redAtPerSecondSquared\": %.9g, \"lastStepMax\": %.9g, \"lastStepMean\": %.9g},\n",
                debugVgsTemporalSampleCount, (double)VGS_TEMPORAL_CURVATURE_RED_AT,
                (double)debugVgsTemporalLastMax, (double)debugVgsTemporalLastMean);
        fprintf(file, "  \"curvatureByLevelCsv\": \"%s\",\n", state->curvature_levels_path);
        if (state->rip_test)
            fprintf(file, "  \"ripTest\": {\"restEdge\": %.9g, \"gripDepth\": %.9g, \"breakStrain\": %.9g, \"commandedStrain\": %.9g, \"firstBreakStep\": %d, \"brokenNodes\": %d, \"brokenLeaves\": %d, \"historyCsv\": \"%s\"},\n",
                    (double)(state->rip_right_x - state->rip_left_x),
                    (double)state->rip_grip_depth,
                    (double)STRAIN_BREAK_THRESHOLD,
                    (double)state->rip_commanded_strain,
                    state->rip_first_break_step,
                    state->rip_broken_nodes, state->rip_broken_leaves,
                    state->rip_history_path);
    }
    if (state->refinement_demo) {
        fprintf(file, "  \"refinement\": {\"step\": %d, \"completed\": %s, \"sharedGridPoints\": 27, \"children\": 8},\n",
                state->refinement_step, state->refinement_complete ? "true" : "false");
    }
    if (state->coarsening_demo) {
        fprintf(file, "  \"coarsening\": {\"step\": %d, \"completed\": %s, \"fineGridPoints\": 27, \"rootCorners\": 8, \"projectionRms\": %.9g},\n",
                state->refinement_step, state->coarsening_complete ? "true" : "false",
                (double)state->coarsening_projection_rms);
    }
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
            "  \"sleepSnapSleepingObserved\": %s,\n"
            "  \"sleepSnapOriginalSucceeded\": %s,\n  \"sleepSnapStaticStep\": %d,\n"
            "  \"sleepSnapQueueAtStatic\": %d,\n  \"sleepSnapRebuiltAfterStatic\": %s,\n"
            "  \"sleepSnapMaxOriginalError\": %.9g,\n"
            "  \"sleepSnapLastSpeed\": %.9g,\n  \"sleepSnapLastStrain\": %.9g,\n"
            "  \"sleepSnapLastShear\": %.9g,\n"
            "  \"supportRemovalSucceeded\": %s,\n  \"activationStep\": %d,"
            "\n  \"dropAtStep60\": %.9g,\n",
            state->activation_succeeded ? "true" : "false",
            state->sleeping_conversion_succeeded ? "true" : "false",
            state->sleep_snap_original_test ? "true" : "false",
            state->sleep_snap_active_observed ? "true" : "false",
            state->sleep_snap_sleeping_observed ? "true" : "false",
            state->sleep_snap_original_succeeded ? "true" : "false",
            state->sleep_snap_static_step,
            state->sleep_snap_queue_at_static,
            state->sleep_snap_rebuilt_after_static ? "true" : "false",
            (double)state->sleep_snap_max_original_error,
            (double)state->sleep_snap_last_speed,
            (double)state->sleep_snap_last_strain,
            (double)state->sleep_snap_last_shear,
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

    fprintf(file,
            "  \"pbfDiagnostics\": {\"enabled\": %s, \"expectedParticles\": %d,"
            " \"dropAtStep120\": %.9g, \"escapedParticles\": %d,"
            " \"belowContainerParticles\": %d},\n",
            state->pbf_container ? "true" : "false",
            state->expected_fluid_particles,
            (double)state->fluid_drop_at_step_120,
            state->escaped_fluid_particles,
             state->below_container_particles);

    double benchmark_physics_ms = state->benchmark_measured_steps > 0
        ? state->benchmark_physics_total_ms / (double)state->benchmark_measured_steps : 0.0;
    double benchmark_render_ms = state->benchmark_measured_steps > 0
        ? state->benchmark_render_total_ms / (double)state->benchmark_measured_steps : 0.0;
    double benchmark_frame_ms = state->benchmark_measured_steps > 0
        ? state->benchmark_frame_total_ms / (double)state->benchmark_measured_steps : 0.0;
    double benchmark_fps = benchmark_frame_ms > 0.0 ? 1000.0 / benchmark_frame_ms : 0.0;
    double benchmark_gpu_compute_ms = state->benchmark_measured_steps > 0
        ? state->benchmark_gpu_compute_total_ms / (double)state->benchmark_measured_steps : 0.0;
    double benchmark_gpu_hull_ms = state->benchmark_measured_steps > 0
        ? state->benchmark_gpu_hull_total_ms / (double)state->benchmark_measured_steps : 0.0;
    fprintf(file,
            "  \"benchmark\": {\"enabled\": %s, \"fluidParticles\": %d,"
            " \"hullRender\": true, \"hullBackend\": \"%s\", \"warmupSteps\": %d, \"measuredSteps\": %d,"
            " \"averagePhysicsMs\": %.9g, \"averageRenderMs\": %.9g,"
            " \"averageFrameMs\": %.9g, \"framesPerSecond\": %.9g,"
            " \"gpuComputeMs\": %.9g, \"gpuHullMs\": %.9g,"
            " \"gpuPackMs\": %.9g, \"gpuUploadMs\": %.9g, \"gpuReadbackMs\": %.9g,"
            " \"gpuUploadBytesPerFrame\": %.9g, \"gpuReadbackBytesPerFrame\": %.9g,"
            " \"gpuDispatchesPerFrame\": %.9g, \"gpuHashCapacity\": %d,"
            " \"gpuNeighborOverflowMax\": %u},\n",
            state->pbf_scaling_benchmark ? "true" : "false",
            state->expected_fluid_particles,
            state->benchmark_gpu_hull ? "gpu-compute" : "cpu",
            state->benchmark_warmup_steps, state->benchmark_measured_steps,
            benchmark_physics_ms, benchmark_render_ms, benchmark_frame_ms, benchmark_fps,
            benchmark_gpu_compute_ms, benchmark_gpu_hull_ms,
            state->benchmark_measured_steps > 0 ? state->benchmark_gpu_pack_total_ms / state->benchmark_measured_steps : 0.0,
            state->benchmark_measured_steps > 0 ? state->benchmark_gpu_upload_total_ms / state->benchmark_measured_steps : 0.0,
            state->benchmark_measured_steps > 0 ? state->benchmark_gpu_readback_total_ms / state->benchmark_measured_steps : 0.0,
            state->benchmark_measured_steps > 0 ? (double)state->benchmark_gpu_upload_bytes / state->benchmark_measured_steps : 0.0,
            state->benchmark_measured_steps > 0 ? (double)state->benchmark_gpu_readback_bytes / state->benchmark_measured_steps : 0.0,
            state->benchmark_measured_steps > 0 ? (double)state->benchmark_gpu_dispatches / state->benchmark_measured_steps : 0.0,
            state->benchmark_gpu_hash_capacity,
            state->benchmark_gpu_neighbor_overflow);

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
                " \"taggedRendered\": %d,"
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
                sample->tagged_rendered,
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
    if (!debug_make_directory_recursive(state->output_dir)) {
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
    activePlayers = 1;
    state.expected_voxels = DEBUG_EXPECTED_VOXELS;
    state.camera.position = (Vector3){ 9.0f, 8.5f, 12.0f };
    state.camera.target = (Vector3){ 0.0f, 4.0f, 0.0f };
    state.camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    state.camera.fovy = 45.0f;
    state.camera.projection = CAMERA_PERSPECTIVE;

    if (!debug_prepare_output_directory(&state)) return 4;
    FILE *timing_file = NULL;
    if (getenv("FPS_AMR_RECORD_TIMINGS")) {
        char timing_path[DEBUG_PATH_CAPACITY];
        if (!debug_join_path(timing_path, sizeof(timing_path),
                             state.output_dir, "frame-timings.csv") ||
            !(timing_file = fopen(timing_path, "w"))) {
            fprintf(stderr, "Unable to create frame timing CSV\n");
            return 4;
        }
        fputs("step,physicsMs,debugFrameMs\n", timing_file);
    }
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
        debug_record_rip_step(&state, 0);
        debug_record_rip_parent_child_strain(&state, 0);
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
            double frame_started = GetTime();
            debug_physics_step(&state);
            double physics_finished = GetTime();
            if (physicsBackend.active == PHYSICS_BACKEND_GPU_GL43 &&
                gpuPhysics.last_neighbor_overflow > state.benchmark_gpu_neighbor_overflow) {
                state.benchmark_gpu_neighbor_overflow = gpuPhysics.last_neighbor_overflow;
            }
            if (state.pbf_scaling_benchmark) {
                debug_render_benchmark_frame(&state);
            }
            double frame_finished = GetTime();
            if (state.pbf_scaling_benchmark && step > state.benchmark_warmup_steps) {
                state.benchmark_physics_total_ms += (physics_finished - frame_started) * 1000.0;
                state.benchmark_render_total_ms += (frame_finished - physics_finished) * 1000.0;
                state.benchmark_frame_total_ms += (frame_finished - frame_started) * 1000.0;
                if (physicsBackend.active == PHYSICS_BACKEND_GPU_GL43) {
                    state.benchmark_gpu_compute_total_ms += gpuPhysics.last_gpu_compute_ms;
                    state.benchmark_gpu_hull_total_ms += gpuPhysics.last_hull_gpu_ms;
                    state.benchmark_gpu_pack_total_ms += gpuPhysics.last_pack_ms;
                    state.benchmark_gpu_upload_total_ms += gpuPhysics.last_upload_ms;
                    state.benchmark_gpu_readback_total_ms += gpuPhysics.last_readback_ms;
                    state.benchmark_gpu_upload_bytes += gpuPhysics.last_upload_bytes;
                    state.benchmark_gpu_readback_bytes += gpuPhysics.last_readback_bytes;
                    state.benchmark_gpu_dispatches += gpuPhysics.last_dispatches;
                    state.benchmark_gpu_hash_capacity = gpuPhysics.hash_capacity;
                    state.benchmark_gpu_hull = state.benchmark_gpu_hull || gpuPhysics.last_hull_used_gpu;
                }
                state.benchmark_measured_steps++;
            }
            debug_validate_step(&state, step);
            debug_record_rip_step(&state, step);
            debug_record_rip_parent_child_strain(&state, step);
            debug_gif_capture(&state, step);
            if (debug_should_capture_step(step)) {
                char file_name[64];
                snprintf(file_name, sizeof(file_name), "step-%04d.png", step);
                debug_capture(&state, file_name, "simulation", step);
            }
            if (timing_file)
                fprintf(timing_file, "%d,%.9g,%.9g\n", step,
                        (physics_finished - frame_started) * 1000.0,
                        (GetTime() - frame_started) * 1000.0);
        }
        debug_finalize_validation(&state);
    }

    if (timing_file && fclose(timing_file) != 0)
        debug_add_failure(&state, "failed to close frame timing CSV");
    debug_gif_end(&state);
    if (state.rip_history_file) {
        if (fclose(state.rip_history_file) != 0)
            debug_add_failure(&state, "failed to close rip history CSV");
        state.rip_history_file = NULL;
    }
    if (state.rip_strain_levels_file) {
        if (fclose(state.rip_strain_levels_file) != 0)
            debug_add_failure(&state, "failed to close rip parent/child strain CSV");
        state.rip_strain_levels_file = NULL;
    }
    if (state.curvature_levels_file) {
        if (fclose(state.curvature_levels_file) != 0)
            debug_add_failure(&state, "failed to close curvature distribution CSV");
        state.curvature_levels_file = NULL;
    }
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
    if (state.pbf_scaling_benchmark && state.benchmark_measured_steps > 0) {
        double frame_ms = state.benchmark_frame_total_ms / (double)state.benchmark_measured_steps;
        fprintf(stderr,
                "debug-benchmark particles=%d backend=%s hull=surface physicsMs=%.6f renderMs=%.6f frameMs=%.6f fps=%.3f measured=%d warmup=%d\n",
                state.expected_fluid_particles, physics_backend_name(physicsBackend.active),
                state.benchmark_physics_total_ms / (double)state.benchmark_measured_steps,
                state.benchmark_render_total_ms / (double)state.benchmark_measured_steps,
                frame_ms, frame_ms > 0.0 ? 1000.0 / frame_ms : 0.0,
                state.benchmark_measured_steps, state.benchmark_warmup_steps);
    }
    for (int i = 0; i < state.failure_count; ++i) {
        fprintf(stderr, "debug-failure: %s\n", state.failures[i]);
    }
    UnloadRenderTexture(state.target);
    return passed ? 0 : 4;
}
