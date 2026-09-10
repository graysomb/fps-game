#include "raylib.h"
#include "../greek_grammar.h"
#include "../greek_rasterizer.h"
#include "../greek_grammar.c"
#include <stdio.h>

#define MAX_PREVIEW_VOXELS 32768
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

int main(int argc, char **argv) {
    SetConfigFlags(FLAG_WINDOW_HIDDEN);
    InitWindow(1280, 720, "Temple Renderer");
    if (!IsWindowReady()) {
        fprintf(stderr, "Failed to initialize window\n");
        return 1;
    }

    uint32_t seed = 1337;
    TempleStage stage = TEMPLE_STAGE_PERIPTERAL;
    if (argc > 2) seed = (uint32_t)atoi(argv[2]);
    if (argc > 3) stage = (TempleStage)atoi(argv[3]);

    TemplePlan plan = generate_greek_temple(seed, stage, 2);
    int center_g = 40;
    rasterize_temple_plan(&plan, center_g, center_g, 3, plot_preview_voxel);
    printf("Rasterized %d voxels for temple (stage %d, seed %u)\n", preview_voxel_count, stage, seed);

    Camera3D camera = { 0 };
    // 3/4 isometric perspective looking at entrance and flank
    camera.position = (Vector3){ 26.0f, 18.0f, 26.0f };
    camera.target = (Vector3){ 0.0f, 4.0f, 0.0f };
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera.fovy = 40.0f;
    camera.projection = CAMERA_PERSPECTIVE;

    RenderTexture2D target = LoadRenderTexture(1280, 720);
    BeginTextureMode(target);
        ClearBackground((Color){ 140, 190, 230, 255 }); // Sky blue
        BeginMode3D(camera);
            // Green arena floor
            DrawCube((Vector3){ 0.0f, 0.5f, 0.0f }, 40.0f, 1.0f, 40.0f, (Color){ 75, 105, 70, 255 });
            DrawCubeWires((Vector3){ 0.0f, 0.5f, 0.0f }, 40.0f, 1.0f, 40.0f, (Color){ 55, 80, 50, 255 });

            float vs = 0.5f;
            for (int i = 0; i < preview_voxel_count; ++i) {
                DrawCube(preview_voxels[i].pos, vs, vs, vs, preview_voxels[i].color);
                DrawCubeWires(preview_voxels[i].pos, vs, vs, vs, ColorBrightness(preview_voxels[i].color, -0.2f));
            }
        EndMode3D();
    EndTextureMode();

    Image img = LoadImageFromTexture(target.texture);
    ImageFlipVertical(&img);
    const char *out_path = (argc > 1) ? argv[1] : "temple_engine_actual.png";
    ExportImage(img, out_path);
    printf("Exported actual in-engine screenshot to %s\n", out_path);

    UnloadImage(img);
    UnloadRenderTexture(target);
    CloseWindow();
    return 0;
}
