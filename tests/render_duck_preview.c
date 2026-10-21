/* Standalone visual harness; deliberately does not change matches or players. */
#include "../duck_renderer.c"
#include <stdio.h>
#include <stdlib.h>

static DuckPose preview_pose(float time,bool moving) {
    float angle=time*2*PI/8;
    return (DuckPose){.position=moving?(Vector3){1.05f*cosf(angle),0,.65f*sinf(angle)}:(Vector3){0},
        .yaw=moving?atan2f(1.05f*sinf(angle),-.65f*cosf(angle))*RAD2DEG:-28,
        .scale=1,.time=time,.speed=moving?.7f:0};
}
static void scene(float time,bool moving,bool wake,int view) {
    Camera3D camera={{2.8f,2.1f,-3.7f},{0,.28f,0},{0,1,0},39,CAMERA_PERSPECTIVE};
    if(view>=0) {
        float angle=view*PI/2;
        if(view==3)angle=-PI/4;
        camera.position=(Vector3){sinf(angle)*2.8f,1.25f,-cosf(angle)*2.8f};
    }
    BeginMode3D(camera);
    DrawPlane((Vector3){0,-.005f,0},(Vector2){200,200},(Color){173,183,178,255});
    for(int n=-5;n<=5;n++) {
        DrawLine3D((Vector3){n,0,-5},(Vector3){n,0,5},(Color){153,164,159,255});
        DrawLine3D((Vector3){-5,0,n},(Vector3){5,0,n},(Color){153,164,159,255});
    }
    DuckPose pose=preview_pose(time,moving);
    if(view>=0)pose.yaw=0;
    duck_draw_floor(pose,wake);duck_draw(pose);
    EndMode3D();
    DrawText("MALLARD",24,21,25,(Color){26,57,49,255});
    DrawText(moving?"Swimming on solid ground":"Folded wings / tucked feet",25,52,15,(Color){57,79,69,255});
}
static void save(RenderTexture2D rt,const char *out,const char *name) {
    char path[1024];snprintf(path,sizeof(path),"%s/%s.png",out,name);
    Image image=LoadImageFromTexture(rt.texture);ImageFlipVertical(&image);
    assert(ExportImage(image,path));UnloadImage(image);
}
static void check_geometry(void) {
    assert(duck_renderer_triangle_count()>500 && duck_renderer_triangle_count()<2000);
    for(int part=0;part<PART_COUNT;part++) {
        Mesh *mesh=&duckParts[part].meshes[0];
        for(int v=0;v<mesh->vertexCount;v++) {
            Vector3 p={mesh->vertices[v*3],mesh->vertices[v*3+1],mesh->vertices[v*3+2]};
            Vector3 n={mesh->normals[v*3],mesh->normals[v*3+1],mesh->normals[v*3+2]};
            assert(isfinite(p.x)&&isfinite(p.y)&&isfinite(p.z));
            assert(isfinite(n.x)&&isfinite(n.y)&&isfinite(n.z));
            assert(Vector3Length(n)>.99f && Vector3Length(n)<1.01f);
        }
    }
    assert(!valid_pose((DuckPose){.scale=NAN}));
    assert(valid_pose(preview_pose(0,true)));
}
static void check_occlusion(RenderTexture2D rt) {
    Camera3D camera={{0,1,-4},{0,.4f,0},{0,1,0},40,CAMERA_PERSPECTIVE};
    Image captures[2];
    for(int i=0;i<2;i++) {
        BeginTextureMode(rt);ClearBackground(BLACK);BeginMode3D(camera);
        DrawCube((Vector3){0,1,-1.8f},6,4,.2f,GRAY);
        if(i) {DuckPose p=preview_pose(0,false);duck_draw_floor(p,true);duck_draw(p);}
        EndMode3D();EndTextureMode();captures[i]=LoadImageFromTexture(rt.texture);
    }
    Color *a=LoadImageColors(captures[0]),*b=LoadImageColors(captures[1]);
    assert(!memcmp(a,b,rt.texture.width*rt.texture.height*sizeof(Color)));
    UnloadImageColors(a);UnloadImageColors(b);UnloadImage(captures[0]);UnloadImage(captures[1]);
}
int main(int argc,char **argv) {
    bool capture=argc>1 && !strcmp(argv[1],"--capture");
    bool checks=argc>1 && !strcmp(argv[1],"--checks");
    bool fallback=argc>1 && !strcmp(argv[1],"--fallback");
    const char *out=argc>2?argv[2]:"artifacts/duck";
    SetConfigFlags((capture||checks||fallback?FLAG_WINDOW_HIDDEN:0)|FLAG_MSAA_4X_HINT);
    InitWindow(800,600,"Mallard swimming preview");assert(IsWindowReady());
    assert(duck_renderer_init(fallback?"missing-duck-shaders":NULL));
    assert(duck_renderer_has_shader()!=fallback);
    if(!fallback)assert(floorShaded);
    check_geometry();
    unsigned int vao=duckParts[BODY].meshes[0].vaoId;
    assert(duck_renderer_init(NULL));assert(duckParts[BODY].meshes[0].vaoId==vao);
    RenderTexture2D rt=LoadRenderTexture(800,600);
    check_occlusion(rt);
    if(capture||fallback) {
        for(int view=0;view<4;view++) {
            BeginTextureMode(rt);ClearBackground(RAYWHITE);scene(.8f,false,false,view);EndTextureMode();
            char name[64];snprintf(name,sizeof(name),"%sview-%d",fallback?"fallback-":"",view);save(rt,out,name);
        }
        if(capture)for(int frame=0;frame<160;frame++) {
            BeginTextureMode(rt);ClearBackground(RAYWHITE);scene(frame/20.0f,true,true,-1);EndTextureMode();
            char name[64];snprintf(name,sizeof(name),"swim-%03d",frame);save(rt,out,name);
        }
    } else if(!checks) {
        SetTargetFPS(60);float time=0;bool paused=false,moving=true,wake=true;
        while(!WindowShouldClose()) {
            if(IsKeyPressed(KEY_SPACE))paused=!paused;
            if(IsKeyPressed(KEY_M))moving=!moving;
            if(IsKeyPressed(KEY_W))wake=!wake;
            if(IsKeyPressed(KEY_R))time=0;
            if(IsKeyDown(KEY_RIGHT))time+=GetFrameTime();
            if(IsKeyDown(KEY_LEFT))time=fmaxf(0,time-GetFrameTime());
            if(!paused)time+=GetFrameTime();
            BeginDrawing();ClearBackground(RAYWHITE);scene(time,moving,wake,-1);
            DrawText("SPACE pause   arrows scrub   M swim/idle   W wake   R reset",24,572,14,DARKGRAY);
            EndDrawing();
        }
    }
    printf("Mallard: %d cached triangles; geometry, occlusion, initialization passed (%s).\n",triangles,fallback?"fallback":"shader");
    UnloadRenderTexture(rt);duck_renderer_shutdown();duck_renderer_shutdown();
    assert(!ready && !triangles);
    assert(duck_renderer_init(NULL));duck_renderer_shutdown();CloseWindow();return 0;
}
