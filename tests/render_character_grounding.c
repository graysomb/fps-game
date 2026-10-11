#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>

/* Deterministic capture using the production movement and foot solver.
 * Create artifacts/rat/grounded/frames before running this executable. */
int main(void) {
    SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);
    InitWindow(640,640,"Rat ground contact");assert(IsWindowReady());
    init_static_hash();assert(character_renderer_init());
    RenderTexture2D rt=LoadRenderTexture(640,640);
    const char *names[]={"walk-stop","turn","jump-land","platform-edge"};
    for(int scene=0;scene<4;scene++) {
        clear_world_voxels();
        if(scene==3)for(int x=-3;x<3;x++)for(int z=-3;z<3;z++)for(int y=0;y<2;y++)
            assert(add_static_voxel_at_grid(x,y,z,GRAY,0)>=0);
        rebuild_static_hash_if_dirty();
        netTransport.role=NET_ROLE_OFFLINE;activePlayers=1;playerInput[0]=INPUT_TYPE_KEYBOARD;
        players[0]=(Player){.pos={0,scene==3?2:1,0},.onGround=true,.tetherVoxel=-1};
        CharacterAnimation animation={0};CharacterPose pose={0};
        for(int frame=0;frame<120;frame++) {
            for(int sub=0;sub<2;sub++) {
                int tick=frame*2+sub;float t=tick/60.0f;
                players[0].vel.x=(scene==0 && t>.4f && t<2.8f)?1.2f:(scene==3 && t>.4f && t<3)?1:0;
                if(scene==2 && tick==40){players[0].vel.y=JUMP_SPEED;players[0].onGround=false;}
                if(scene==1)players[0].yaw=t<.5f?0:fminf(180,(t-.5f)*90);
                net_integrate_player(&players[0],1.0f/60,true);
                pose=character_pose(players[0].pos,players[0].yaw,0,fabsf(players[0].vel.x),players[0].onGround,false,SKYBLUE,t);
                character_animate(&animation,&pose,1.0f/60,character_ground_query,0);
            }
            float x=players[0].pos.x;
            Camera3D camera={.position={x-2.5f,2.5f,-3.5f},.target={x,scene==3?1.15f:.85f,0},.up={0,1,0},.fovy=43};
            if(scene==2){camera.position=(Vector3){x-5,4,-7};camera.target.y=2;}
            BeginTextureMode(rt);ClearBackground((Color){214,216,218,255});BeginMode3D(camera);
            DrawPlane((Vector3){0,-.003f,0},(Vector2){40,40},(Color){190,193,196,255});
            for(int line=-20;line<=20;line++) {
                DrawLine3D((Vector3){line*.5f,.005f,-10},(Vector3){line*.5f,.005f,10},(Color){160,164,168,255});
                DrawLine3D((Vector3){-10,.005f,line*.5f},(Vector3){10,.005f,line*.5f},(Color){160,164,168,255});
            }
            if(scene==3){DrawCube((Vector3){0,.5f,0},3,1,3,(Color){125,133,140,255});DrawCubeWires((Vector3){0,.5f,0},3,1,3,DARKGRAY);}
            character_draw(pose,false);
            for(int side=-1;side<=1;side+=2)character_draw_arm(
                Vector3Transform((Vector3){side*.17f,.11f,0},character_body_transform(pose)),
                Vector3Transform((Vector3){side*.18f,.03f,-.16f},character_body_transform(pose)),pose.gear);
            character_draw_visor(pose);EndMode3D();
            DrawText(names[scene],24,20,26,DARKGRAY);
            DrawText(players[0].onGround?"GROUNDED":"AIRBORNE",24,54,18,DARKGRAY);EndTextureMode();
            char path[256];snprintf(path,sizeof(path),"artifacts/rat/grounded/frames/%s-%03d.png",names[scene],frame);
            Image image=LoadImageFromTexture(rt.texture);ImageFlipVertical(&image);assert(ExportImage(image,path));UnloadImage(image);
        }
    }
    UnloadRenderTexture(rt);character_renderer_shutdown();CloseWindow();return 0;
}
