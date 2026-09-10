#include "raylib.h"
#include "../greek_grammar.h"
#include "../greek_rasterizer.h"
#include "../greek_grammar.c"
#include <stdio.h>

#define MAX_PREVIEW_VOXELS 262144
typedef struct {
    Vector3 pos;
    Color color;
} SimpleVoxel;

static SimpleVoxel preview_voxels[MAX_PREVIEW_VOXELS];
static int preview_voxel_count = 0;

static void plot_preview_voxel(int gx, int gy, int gz, Color c) {
    if (preview_voxel_count >= MAX_PREVIEW_VOXELS) return;
    float voxel_size = 0.5f;
    float floor_size = 20.0f;
    float px = (gx + 0.5f) * voxel_size - floor_size;
    float py = (gy + 0.5f) * voxel_size;
    float pz = (gz + 0.5f) * voxel_size - floor_size;
    preview_voxels[preview_voxel_count++] = (SimpleVoxel){ (Vector3){ px, py, pz }, c };
}

static void render_stage_image(TempleStage stage, uint32_t seed, const char *out_path, const char *title) {
    preview_voxel_count = 0;
    TemplePlan plan = generate_greek_temple(seed, stage, 2);
    int center_g = 40;
    rasterize_temple_plan(&plan, center_g, center_g, 3, plot_preview_voxel);

    // Adaptive camera position based on stage size
    float cam_dist = (stage == TEMPLE_STAGE_SANCTUARY) ? 38.0f : (stage == TEMPLE_STAGE_PERIPTERAL) ? 26.0f : 20.0f;
    float cam_y    = (stage == TEMPLE_STAGE_SANCTUARY) ? 29.0f : (stage == TEMPLE_STAGE_PERIPTERAL) ? 17.0f : 14.0f;
    float target_z = 0.0f;
    float fovy     = (stage == TEMPLE_STAGE_SANCTUARY) ? 52.0f : 45.0f;

    Camera3D camera = { 0 };
    camera.position = (Vector3){ cam_dist, cam_y, cam_dist };
    camera.target = (Vector3){ 0.0f, 3.5f, target_z };
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera.fovy = fovy;
    camera.projection = CAMERA_PERSPECTIVE;

    RenderTexture2D target = LoadRenderTexture(1280, 720);
    BeginTextureMode(target);
        ClearBackground((Color){ 140, 190, 230, 255 }); // Mediterranean sky
        BeginMode3D(camera);
            // Arena Floor
            DrawCube((Vector3){ 0.0f, 0.5f, 0.0f }, 40.0f, 1.0f, 40.0f, (Color){ 75, 105, 70, 255 });
            DrawCubeWires((Vector3){ 0.0f, 0.5f, 0.0f }, 40.0f, 1.0f, 40.0f, (Color){ 55, 80, 50, 255 });

            float vs = 0.5f;
            for (int i = 0; i < preview_voxel_count; ++i) {
                DrawCube(preview_voxels[i].pos, vs, vs, vs, preview_voxels[i].color);
                DrawCubeWires(preview_voxels[i].pos, vs, vs, vs, ColorBrightness(preview_voxels[i].color, -0.2f));
            }
        EndMode3D();

        // On-screen watermark / title banner
        DrawRectangle(20, 20, 500, 75, Fade(BLACK, 0.65f));
        DrawRectangleLines(20, 20, 500, 75, RAYWHITE);
        DrawText(title, 35, 32, 22, RAYWHITE);
        DrawText(TextFormat("Voxels: %d | Seed: %u", preview_voxel_count, seed), 35, 62, 18, GOLD);

    EndTextureMode();

    Image img = LoadImageFromTexture(target.texture);
    ImageFlipVertical(&img);
    ExportImage(img, out_path);
    printf("Rendered [%s] -> %s (%d voxels)\n", title, out_path, preview_voxel_count);

    UnloadImage(img);
    UnloadRenderTexture(target);
}

static void render_closeup_view(TempleStage stage, uint32_t seed, Vector3 cam_pos, Vector3 cam_target, float fovy, const char *out_path, const char *title) {
    preview_voxel_count = 0;
    TemplePlan plan = generate_greek_temple(seed, stage, 2);
    int center_g = 40;
    rasterize_temple_plan(&plan, center_g, center_g, 3, plot_preview_voxel);

    Camera3D camera = { 0 };
    camera.position = cam_pos;
    camera.target = cam_target;
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera.fovy = fovy;
    camera.projection = CAMERA_PERSPECTIVE;

    RenderTexture2D target = LoadRenderTexture(1280, 720);
    BeginTextureMode(target);
        ClearBackground((Color){ 140, 190, 230, 255 });
        BeginMode3D(camera);
            DrawCube((Vector3){ 0.0f, 0.5f, 0.0f }, 40.0f, 1.0f, 40.0f, (Color){ 75, 105, 70, 255 });
            DrawCubeWires((Vector3){ 0.0f, 0.5f, 0.0f }, 40.0f, 1.0f, 40.0f, (Color){ 55, 80, 50, 255 });

            float vs = 0.5f;
            for (int i = 0; i < preview_voxel_count; ++i) {
                DrawCube(preview_voxels[i].pos, vs, vs, vs, preview_voxels[i].color);
                DrawCubeWires(preview_voxels[i].pos, vs, vs, vs, ColorBrightness(preview_voxels[i].color, -0.2f));
            }
        EndMode3D();

        DrawRectangle(20, 20, 520, 75, Fade(BLACK, 0.65f));
        DrawRectangleLines(20, 20, 520, 75, RAYWHITE);
        DrawText(title, 35, 32, 22, RAYWHITE);
        DrawText(TextFormat("Voxels: %d | Seed: %u", preview_voxel_count, seed), 35, 62, 18, GOLD);

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
    InitWindow(1280, 720, "Gallery Renderer");
    if (!IsWindowReady()) return 1;

    const char *out_dir = "/Users/migr4479/.gemini/antigravity/brain/86453e72-1dcc-4b9c-ba5a-2222f2dbe7ff";

    // 1. All 5 Maturation Stages
    render_stage_image(TEMPLE_STAGE_CELLA, 1337,
                       TextFormat("%s/stage0_cella.png", out_dir),
                       "Stage 0: Sacred Cella (Naos)");

    render_stage_image(TEMPLE_STAGE_PROSTYLE, 1337,
                       TextFormat("%s/stage1_prostyle.png", out_dir),
                       "Stage 1: Prostyle (Front Porch)");

    render_stage_image(TEMPLE_STAGE_AMPHIPROSTYLE, 1337,
                       TextFormat("%s/stage2_amphi.png", out_dir),
                       "Stage 2: Amphiprostyle (Front+Rear)");

    render_stage_image(TEMPLE_STAGE_PERIPTERAL, 1337,
                       TextFormat("%s/stage3_peripteral.png", out_dir),
                       "Stage 3: Peripteral (Colonnade Wrap)");

    render_stage_image(TEMPLE_STAGE_SANCTUARY, 1337,
                       TextFormat("%s/stage4_sanctuary.png", out_dir),
                       "Stage 4: Sanctuary Precinct & Sacred Garden");

    // 2. Close-up of Courtyard, PBF Pool, and Olive Trees (Axial vista from Propylaea toward Temple)
    render_closeup_view(TEMPLE_STAGE_SANCTUARY, 1337,
                        (Vector3){ 0.0f, 6.5f, 24.0f },
                        (Vector3){ 0.0f, 3.2f, 4.0f },
                        50.0f,
                        TextFormat("%s/sanctuary_courtyard_closeup.png", out_dir),
                        "Sanctuary: Temenos Pool & Sacred Grove Axial Vista");

    // 3. Close-up of The Tholos (Round Monopteros Rotunda, Cypress Grove & Braziers)
    render_closeup_view(TEMPLE_STAGE_SANCTUARY, 1337,
                        (Vector3){ 9.5f, 6.5f, -24.5f },
                        (Vector3){ 0.0f, 3.2f, -15.0f },
                        50.0f,
                        TextFormat("%s/tholos_rotunda_view.png", out_dir),
                        "Sacred Glade: The Tholos Rotunda & Cypress Grove");

    // 4. Close-up of Garden Flank (PBF Fountain, Marble Exedra & Olive Trees)
    render_closeup_view(TEMPLE_STAGE_SANCTUARY, 1337,
                        (Vector3){ -21.0f, 8.0f, -1.0f },
                        (Vector3){ -14.5f, 2.0f, -5.0f },
                        46.0f,
                        TextFormat("%s/garden_exedra_fountain_view.png", out_dir),
                        "Garden Precinct: PBF Fountain, Exedra & Olive Grove");

    // 5. Different Random Seeds
    render_stage_image(TEMPLE_STAGE_PERIPTERAL, 42,
                       TextFormat("%s/seed_42.png", out_dir),
                       "Peripteral Temple [Seed 42]");

    render_stage_image(TEMPLE_STAGE_PERIPTERAL, 7777,
                       TextFormat("%s/seed_7777.png", out_dir),
                       "Peripteral Temple [Seed 7777]");

    CloseWindow();
    return 0;
}
