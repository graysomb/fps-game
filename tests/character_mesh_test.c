/* Inspect generated data directly without adding test-only public APIs. */
#include "../character_renderer.c"
#include <assert.h>
#include <stdio.h>
int main(int argc,char **argv) {
    (void)argv;
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(160,120,"Rat mesh checks");
    bool fallback=argc>1;
    if(fallback) assert(ChangeDirectory(".build"));
    assert(character_renderer_init());assert(shaded!=fallback);
    for(int k=0;k<CR_COUNT;k++) {
        Mesh m=parts[k].meshes[0];assert(m.vertexCount>0 && m.vertexCount<=2048);
        for(int i=0;i<m.vertexCount*3;i++) {assert(isfinite(m.vertices[i]));assert(isfinite(m.normals[i]));}
        BoundingBox box=GetMeshBoundingBox(m);
        assert(box.min.y>=-.25f && box.max.y<=.25f);
    }
    CharacterPose pose=character_pose((Vector3){0},0,90,4,true,false,SKYBLUE,1);
    assert(pose.pitch==25 && isfinite(pose.stride));
    BeginDrawing();ClearBackground(BLACK);BeginMode3D((Camera3D){.position={1,1,-2},.up={0,1,0},.fovy=45});
    character_draw(pose,false);EndMode3D();EndDrawing();
    character_renderer_shutdown();character_renderer_shutdown();CloseWindow();
    printf("Character mesh checks passed (%s)\n",fallback?"fallback":"shader");return 0;
}
