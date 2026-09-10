#include "raylib.h"
#include "../forerunner_grammar.h"
#include "../forerunner_rasterizer.h"
#include "../forerunner_grammar.c"
#include <stdio.h>
#include <math.h>

#define MAX_PREVIEW_VOXELS 524288

typedef struct {
    Vector3 pos;
    Color color;
} SimpleVoxel;

static SimpleVoxel preview_voxels[MAX_PREVIEW_VOXELS];
static int preview_voxel_count = 0;

static void plot_forerunner_preview(int gx, int gy, int gz, Color c) {
    if (preview_voxel_count >= MAX_PREVIEW_VOXELS) return;
    float voxel_size = 0.5f;
    float floor_size = 25.0f; // 50m arena
    float px = (gx + 0.5f) * voxel_size - floor_size;
    float py = (gy + 0.5f) * voxel_size;
    float pz = (gz + 0.5f) * voxel_size - floor_size;
    preview_voxels[preview_voxel_count++] = (SimpleVoxel){ (Vector3){ px, py, pz }, c };
}

static void render_forerunner_shot(ForerunnerStage stage, uint32_t seed,
                                  Vector3 cam_pos, Vector3 cam_target, float fovy,
                                  const char *out_path, const char *title, const char *subtitle) {
    preview_voxel_count = 0;
    ForerunnerPlan plan = generate_forerunner_structure(seed, stage);
    int center_g = 50; // 25m / 0.5m = 50 voxels
    rasterize_forerunner_plan(&plan, center_g, center_g, 10, plot_forerunner_preview);

    Camera3D camera = { 0 };
    camera.position = cam_pos;
    camera.target = cam_target;
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera.fovy = fovy;
    camera.projection = CAMERA_PERSPECTIVE;

    RenderTexture2D target = LoadRenderTexture(1280, 720);
    BeginTextureMode(target);
        // Halo Ring-World Atmosphere: Deep sci-fi slate/azure sky
        ClearBackground((Color){ 45, 60, 75, 255 });

        BeginMode3D(camera);
            // Installation ground plane
            DrawCube((Vector3){ 0.0f, 4.5f, 0.0f }, 60.0f, 1.0f, 60.0f, (Color){ 55, 62, 70, 255 });
            DrawCubeWires((Vector3){ 0.0f, 4.5f, 0.0f }, 60.0f, 1.0f, 60.0f, (Color){ 40, 48, 55, 255 });

            float vs = 0.5f;
            for (int i = 0; i < preview_voxel_count; ++i) {
                DrawCube(preview_voxels[i].pos, vs, vs, vs, preview_voxels[i].color);
                DrawCubeWires(preview_voxels[i].pos, vs, vs, vs, ColorBrightness(preview_voxels[i].color, -0.20f));
            }
        EndMode3D();

        // Brutalist UI telemetry banner
        DrawRectangle(20, 20, 700, 85, Fade((Color){ 20, 28, 38, 255 }, 0.85f));
        DrawRectangleLines(20, 20, 700, 85, (Color){ 45, 225, 255, 255 });
        DrawText(title, 35, 30, 22, (Color){ 230, 245, 255, 255 });
        DrawText(subtitle, 35, 56, 16, (Color){ 120, 200, 235, 255 });
        DrawText(TextFormat("Nodes: %d | Voxels: %d | Seed: %u | Stage: %d",
                            plan.node_count, preview_voxel_count, seed, stage),
                 35, 78, 15, (Color){ 255, 175, 50, 255 });

    EndTextureMode();

    Image img = LoadImageFromTexture(target.texture);
    ImageFlipVertical(&img);
    ExportImage(img, out_path);
    printf("Rendered [%s] -> %s (%d voxels)\n", title, out_path, preview_voxel_count);

    UnloadImage(img);
    UnloadRenderTexture(target);
}

int main(void) {
    SetConfigFlags(FLAG_WINDOW_HIDDEN);
    InitWindow(1280, 720, "Forerunner Gallery Renderer");
    if (!IsWindowReady()) return 1;

    const char *out_dir = "/Users/migr4479/.gemini/antigravity/brain/86453e72-1dcc-4b9c-ba5a-2222f2dbe7ff";

    // 1. Apex Cartographer Installation (Aerial Overview)
    render_forerunner_shot(FORERUNNER_STAGE_CARTOGRAPHER, 1337,
                           (Vector3){ 24.0f, 22.0f, 24.0f },
                           (Vector3){ 0.0f, 6.0f, 0.0f },
                           48.0f,
                           TextFormat("%s/forerunner_cartographer_overview.png", out_dir),
                           "Installation 04: The Apex Cartographer",
                           "Canted Chevron Pylons, Suspended Hard-Light Bridge & Hovering Gravity Core");

    // 2. First-Person Skybridge & Chasm Crossing Vista
    render_forerunner_shot(FORERUNNER_STAGE_CARTOGRAPHER, 1337,
                           (Vector3){ 0.0f, 9.2f, 8.5f },
                           (Vector3){ 0.0f, 8.8f, -10.0f },
                           54.0f,
                           TextFormat("%s/forerunner_bridge_chasm_vista.png", out_dir),
                           "Abyssal Transit: Suspended Hard-Light Bridge",
                           "Looking North across the Hard-Light Span toward the Core & Telemetry Spires");

    // 3. Upward Monumental Chevron Gateway Portal
    render_forerunner_shot(FORERUNNER_STAGE_GATEWAY, 1337,
                           (Vector3){ 0.0f, 10.0f, 32.0f },
                           (Vector3){ 0.0f, 11.0f, 16.0f },
                           46.0f,
                           TextFormat("%s/forerunner_chevron_portal.png", out_dir),
                           "Monumental Chevron Gateway: Canted Pylons",
                           "25-Degree Inward Rake, Heavy Lintel & Recessed Cyan Energy Grooves");

    // 4. Hovering Gravity Matrix Core & Twin Telemetry Spires
    render_forerunner_shot(FORERUNNER_STAGE_CARTOGRAPHER, 1337,
                           (Vector3){ 0.0f, 9.5f, -2.5f },
                           (Vector3){ 0.0f, 10.0f, -11.0f },
                           50.0f,
                           TextFormat("%s/forerunner_core_spire.png", out_dir),
                           "Apex Matrix: Hovering Gravity Core & Telemetry Spires",
                           "Octahedral Energy Matrix, Containment Ring & Solar Amber Beam Emitters");

    CloseWindow();
    return 0;
}
