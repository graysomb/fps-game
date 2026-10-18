/* Deterministic menu captures: no gameplay world or physics required. */
#define main game_main
#include "../fps_ray.c"
#undef main
#include <assert.h>

int main(void) {
    SetConfigFlags(FLAG_WINDOW_HIDDEN);
    InitWindow(960,720,"Menu preview regression");
    character_renderer_init();
    const int roster[]={-1,ENEMY_TYPE_SWARMER,ENEMY_TYPE_STANDARD,
        ENEMY_TYPE_GOLIATH,ENEMY_TYPE_SHIELDBEARER,ENEMY_TYPE_SAPPER,
        ENEMY_TYPE_LOBBER,ENEMY_TYPE_BUILDER,ENEMY_TYPE_ANCHOR,ENEMY_TYPE_FOUNDRY};
    for(int angle=0;angle<360;angle+=90)for(int n=0;n<10;n++) {
        BeginDrawing();
        clarity_menu_render(roster[n],angle,0);
        EndDrawing();
        Image capture=LoadImageFromTexture(clarityShowcase.texture);
        ImageFlipVertical(&capture);
        assert(capture.data && capture.width==720 && capture.height==720);
        assert(ExportImage(capture,TextFormat("artifacts/menu-preview-%d-%d.png",n,angle)));
        UnloadImage(capture);
    }
    /* Exercise cached resource teardown/recreation as when returning to menus. */
    clarity_menu_shutdown();
    assert(!clarityShowcaseReady && !clarityShowcase.id);
    BeginDrawing();clarity_menu_render(-1,90,0);EndDrawing();
    Image side=LoadImageFromTexture(clarityShowcase.texture);ImageFlipVertical(&side);
    ExportImage(side,"artifacts/menu-preview-rat-side.png");UnloadImage(side);
    UnloadShader(clarityShowcaseShader);clarityShowcaseShader=(Shader){0};
    clarityShowcaseCube.materials[0].shader=(Shader){rlGetShaderIdDefault(),rlGetShaderLocsDefault()};
    BeginDrawing();clarity_menu_render(ENEMY_TYPE_FOUNDRY,180,0);EndDrawing();
    Image fallback=LoadImageFromTexture(clarityShowcase.texture);ImageFlipVertical(&fallback);ExportImage(fallback,"artifacts/menu-preview-fallback.png");UnloadImage(fallback);
    clarity_menu_shutdown();character_renderer_shutdown();weapon_renderer_shutdown();CloseWindow();
    return 0;
}
