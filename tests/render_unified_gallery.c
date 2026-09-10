#include "raylib.h"
#include "../unified_sanctum_grammar.h"
#include "../unified_sanctum_rasterizer.h"
#include "../unified_sanctum_grammar.c"
#include <stdio.h>
#include <math.h>

#define MAX_PREVIEW_VOXELS 1048576

typedef struct {
    Vector3 pos;
    Color color;
} SimpleVoxel;

static SimpleVoxel preview_voxels[MAX_PREVIEW_VOXELS];
static int preview_voxel_count = 0;

static void plot_sanctum_preview(int gx, int gy, int gz, Color c) {
    if (preview_voxel_count >= MAX_PREVIEW_VOXELS) return;
    float voxel_size = 0.5f;
    float floor_size = 40.0f; // 80m arena
    float px = (gx + 0.5f) * voxel_size - floor_size;
    float py = (gy + 0.5f) * voxel_size;
    float pz = (gz + 0.5f) * voxel_size - floor_size;
    preview_voxels[preview_voxel_count++] = (SimpleVoxel){ (Vector3){ px, py, pz }, c };
}

static void render_sanctum_shot(uint32_t seed, int growth_steps,
                                Vector3 cam_pos, Vector3 cam_target, float fovy,
                                const char *out_path, const char *title, const char *subtitle) {
    preview_voxel_count = 0;
    SanctumCitadelPlan plan = generate_unified_sanctum(seed, growth_steps);
    int center_g = 80; // 40m / 0.5m = 80 voxels
    rasterize_sanctum_plan(&plan, center_g, center_g, 10, plot_sanctum_preview, NULL);

    Camera3D camera = { 0 };
    camera.position = cam_pos;
    camera.target = cam_target;
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera.fovy = fovy;
    camera.projection = CAMERA_PERSPECTIVE;

    RenderTexture2D target = LoadRenderTexture(1280, 720);
    BeginTextureMode(target);
        // Atmospheric Twilight Horizon: Deep indigo twilight with subtle cosmic aurora
        ClearBackground((Color){ 28, 36, 48, 255 });

        BeginMode3D(camera);
            // Sacred terrain ground plane
            DrawCube((Vector3){ 0.0f, 4.5f, 0.0f }, 90.0f, 1.0f, 90.0f, (Color){ 44, 48, 54, 255 });
            DrawCubeWires((Vector3){ 0.0f, 4.5f, 0.0f }, 90.0f, 1.0f, 90.0f, (Color){ 32, 36, 42, 255 });

            float vs = 0.5f;
            for (int i = 0; i < preview_voxel_count; ++i) {
                DrawCube(preview_voxels[i].pos, vs, vs, vs, preview_voxels[i].color);
                DrawCubeWires(preview_voxels[i].pos, vs, vs, vs, ColorBrightness(preview_voxels[i].color, -0.18f));
            }
        EndMode3D();

        // Syncretic UI telemetry banner
        const char *topo_names[] = {
            "Stronghold (360-Def Ziggurat)",
            "Abyssal Rift (Chasm & 3 Bridges)",
            "Sunken Crucible (Inverted Colosseum)",
            "Asymmetric Outpost (Bunker & Chokepoint)"
        };
        DrawRectangle(20, 20, 840, 85, Fade((Color){ 16, 22, 30, 255 }, 0.88f));
        DrawRectangleLines(20, 20, 840, 85, (Color){ 45, 225, 255, 255 });
        DrawText(title, 35, 28, 22, (Color){ 240, 248, 255, 255 });
        DrawText(subtitle, 35, 54, 15, (Color){ 255, 205, 90, 255 });
        DrawText(TextFormat("Topology: %s | Spawns: %d | Supplies: %d | JumpPads: %d | Nodes: %d | Voxels: %d",
                            topo_names[plan.topology], plan.spawn_point_count, plan.supply_count,
                            plan.jump_pad_count, plan.node_count, preview_voxel_count),
                 35, 78, 14, (Color){ 120, 220, 255, 255 });

    EndTextureMode();

    Image img = LoadImageFromTexture(target.texture);
    ImageFlipVertical(&img);
    ExportImage(img, out_path);
    printf("Rendered [%s] -> %s (%d voxels, %d nodes, %d spawns, %d supplies)\n",
           title, out_path, preview_voxel_count, plan.node_count, plan.spawn_point_count, plan.supply_count);

    UnloadImage(img);
    UnloadRenderTexture(target);
}

int main(void) {
    SetConfigFlags(FLAG_WINDOW_HIDDEN);
    InitWindow(1280, 720, "Unified Syncretic Firefight Gallery Renderer");
    if (!IsWindowReady()) return 1;

    const char *out_dir = "/Users/migr4479/.gemini/antigravity/brain/86453e72-1dcc-4b9c-ba5a-2222f2dbe7ff";

    // 1. Topology 0: Stronghold (Seed 100, Step 12)
    render_sanctum_shot(100, 12,
                        (Vector3){ 42.0f, 38.0f, 42.0f },
                        (Vector3){ 0.0f, 10.0f, 0.0f },
                        54.0f,
                        TextFormat("%s/firefight_stronghold.png", out_dir),
                        "Firefight Arena 1: The Stronghold Ziggurat",
                        "Concentric 360-Def Hilltop Holdout, 4 Wave Ingress Crypts, Moat Canals & Battlements");

    // 2. Topology 1: Abyssal Rift (Seed 101, Step 12)
    render_sanctum_shot(101, 12,
                        (Vector3){ 38.0f, 26.0f, 0.0f },
                        (Vector3){ 0.0f, 8.0f, 0.0f },
                        55.0f,
                        TextFormat("%s/firefight_abyssal_rift.png", out_dir),
                        "Firefight Arena 2: The Abyssal Rift",
                        "Dual-Plateau Fortress split by Deep Void Chasm, 3 Crossing Bridges & Suspended Altar");

    // 3. Topology 2: Sunken Crucible (Seed 102, Step 12)
    render_sanctum_shot(102, 12,
                        (Vector3){ 38.0f, 32.0f, 38.0f },
                        (Vector3){ 0.0f, 6.0f, 0.0f },
                        52.0f,
                        TextFormat("%s/firefight_sunken_crucible.png", out_dir),
                        "Firefight Arena 3: The Sunken Crucible",
                        "Inverted Colosseum Battle Bowl, Central Weapon Island, Perimeter Promenade & 4 Jump Pads");

    // 4. Topology 3: Asymmetric Outpost (Seed 103, Step 12)
    render_sanctum_shot(103, 12,
                        (Vector3){ -36.0f, 30.0f, 28.0f },
                        (Vector3){ 2.0f, 8.0f, 2.0f },
                        54.0f,
                        TextFormat("%s/firefight_asymmetric_outpost.png", out_dir),
                        "Firefight Arena 4: The Asymmetric Outpost",
                        "Urban Megaron Bunker, Y=28 Watchtower, Zigzagging Curtain Wall, Chokepoint & Staging Grounds");

    CloseWindow();
    return 0;
}
