#include "raylib.h"
#include "../hyperborean_grammar.h"
#include "../hyperborean_rasterizer.h"
#include "../hyperborean_grammar.c"
#include <stdio.h>
#include <math.h>

#define MAX_PREVIEW_VOXELS 524288

typedef struct {
    Vector3 pos;
    Color color;
} SimpleVoxel;

static SimpleVoxel preview_voxels[MAX_PREVIEW_VOXELS];
static int preview_voxel_count = 0;

static void plot_hyper_preview(int gx, int gy, int gz, Color c) {
    if (preview_voxel_count >= MAX_PREVIEW_VOXELS) return;
    float voxel_size = 0.5f;
    float floor_size = 30.0f; // 60m arena
    float px = (gx + 0.5f) * voxel_size - floor_size;
    float py = (gy + 0.5f) * voxel_size;
    float pz = (gz + 0.5f) * voxel_size - floor_size;
    preview_voxels[preview_voxel_count++] = (SimpleVoxel){ (Vector3){ px, py, pz }, c };
}

static void render_hyper_shot(HyperStage stage, uint32_t seed,
                              Vector3 cam_pos, Vector3 cam_target, float fovy,
                              const char *out_path, const char *title, const char *subtitle) {
    preview_voxel_count = 0;
    HyperPlan plan = generate_hyperborean_structure(seed, stage);
    int center_g = 60; // 30m / 0.5m = 60 voxels
    rasterize_hyperborean_plan(&plan, center_g, center_g, 2, plot_hyper_preview);

    Camera3D camera = { 0 };
    camera.position = cam_pos;
    camera.target = cam_target;
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera.fovy = fovy;
    camera.projection = CAMERA_PERSPECTIVE;

    RenderTexture2D target = LoadRenderTexture(1280, 720);
    BeginTextureMode(target);
        // Nordic/Arctic clear dawn sky (Hyperborean light)
        ClearBackground((Color){ 170, 195, 215, 255 });

        BeginMode3D(camera);
            // Grassy moorland terrain
            DrawCube((Vector3){ 0.0f, 0.5f, 0.0f }, 75.0f, 1.0f, 75.0f, (Color){ 65, 95, 55, 255 });
            DrawCubeWires((Vector3){ 0.0f, 0.5f, 0.0f }, 75.0f, 1.0f, 75.0f, (Color){ 45, 75, 40, 255 });

            float vs = 0.5f;
            for (int i = 0; i < preview_voxel_count; ++i) {
                DrawCube(preview_voxels[i].pos, vs, vs, vs, preview_voxels[i].color);
                DrawCubeWires(preview_voxels[i].pos, vs, vs, vs, ColorBrightness(preview_voxels[i].color, -0.22f));
            }
        EndMode3D();

        // Architectural vignette banner
        DrawRectangle(20, 20, 680, 85, Fade(BLACK, 0.75f));
        DrawRectangleLines(20, 20, 680, 85, (Color){ 225, 215, 180, 255 });
        DrawText(title, 35, 30, 22, RAYWHITE);
        DrawText(subtitle, 35, 56, 16, (Color){ 200, 225, 190, 255 });
        DrawText(TextFormat("Nodes: %d | Voxels: %d | Seed: %u | Stage: %d",
                            plan.node_count, preview_voxel_count, seed, stage),
                 35, 78, 15, (Color){ 240, 200, 70, 255 });

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
    InitWindow(1280, 720, "Hyperborean Gallery Renderer");
    if (!IsWindowReady()) return 1;

    const char *out_dir = "/Users/migr4479/.gemini/antigravity/brain/86453e72-1dcc-4b9c-ba5a-2222f2dbe7ff";

    // 1. Complete Hyperborean Sun-Henge Complex (Aerial Overview)
    render_hyper_shot(HYPER_STAGE_FULL_SANCTUM, 1337,
                      (Vector3){ 26.0f, 22.0f, 32.0f },
                      (Vector3){ 0.0f, 3.0f, 2.0f },
                      48.0f,
                      TextFormat("%s/hyperborean_full_complex.png", out_dir),
                      "The Hyperborean Sun-Henge (Apollo Hyperboreus)",
                      "Outer Sarsen Menhirs, Fluted Marble Peristyle, Gold Pediment & PBF Moat");

    // 2. Solstice Processional Avenue Vista (Looking past Heel Stone down dromos)
    render_hyper_shot(HYPER_STAGE_FULL_SANCTUM, 1337,
                      (Vector3){ 2.8f, 3.8f, 26.5f },
                      (Vector3){ 0.0f, 3.0f, 0.0f },
                      48.0f,
                      TextFormat("%s/hyperborean_avenue_vista.png", out_dir),
                      "Solstice Processional Avenue & Heel Stone",
                      "Looking Past the Heel Stone toward the Concentric Peristyle & Great Trilithons");

    // 3. Close-up of the Great Hybrid Trilithons & Doric Pediment
    render_hyper_shot(HYPER_STAGE_FULL_SANCTUM, 1337,
                      (Vector3){ 0.0f, 4.2f, 6.5f },
                      (Vector3){ 0.0f, 4.5f, -6.0f },
                      54.0f,
                      TextFormat("%s/hyperborean_trilithon_pediment.png", out_dir),
                      "The Great Hybrid Trilithon: Megalithic Parthenon",
                      "Rough Cyclopean Orthostats Supporting Classical Doric Pediment & Gold Tympanum");

    // 4. Close-up of the Sunken Marble Tholos & PBF Fluid Moat
    render_hyper_shot(HYPER_STAGE_FULL_SANCTUM, 1337,
                      (Vector3){ 5.0f, 4.0f, 5.0f },
                      (Vector3){ 0.0f, 1.2f, 0.0f },
                      46.0f,
                      TextFormat("%s/hyperborean_tholos_moat.png", out_dir),
                      "Sanctum Core: Sunken Tholos & Sacred PBF Spring",
                      "Tiered Marble Podium, Concentric Fluid Moat & Eternal Flame Brazier");

    CloseWindow();
    return 0;
}
