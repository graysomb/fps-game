#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static void empty(void){
    enemy_cleanup(true);clear_all_chunks();clear_world_voxels();init_static_hash();worldPeriodic=false;
    memset(players,0,sizeof(players));memset(enemySpecial,0,sizeof(enemySpecial));memset(ffTestReviveHeld,0,sizeof(ffTestReviveHeld));
    combat_ui_reset();gameMode=GAME_MODE_FIREFIGHT;gameState=GAME_STATE_PLAYING;firefightHumans=2;activePlayers=5;enemyLab=false;firefightPractice=false;
    netTransport.role=NET_ROLE_OFFLINE;worldFloorHalfSize=40;firefightWaveState=FIREFIGHT_STATE_ACTIVE;
    for(int i=0;i<MAX_COMBATANTS;i++){players[i].pos=(Vector3){20+i,BASE_EYE_HEIGHT,20};players[i].matter=100;players[i].matterMax=100;players[i].kd_ratio=1;players[i].last_damage_time=-1000;combatantLife[i]=i<2?LIFE_ALIVE:LIFE_UNUSED;combatantTeam[i]=i<2?0:1;playerInput[i]=i<2?INPUT_TYPE_KEYBOARD:INPUT_TYPE_BOT_HARD;reviveOwner[i<4?i:0]=-1;}
}
static void free_revive(void){
    empty();players[0].pos=(Vector3){0,BASE_EYE_HEIGHT,0};players[1].pos=(Vector3){1,BASE_EYE_HEIGHT,0};
    players[1].isExposed=true;kill_player(1,4,false,false);players[0].matter=0;ffTestReviveHeld[0]=true;
    ff_update_coop(1);assert(reviveProgress[1]==1 && reviveOwner[1]==0);
    ffTestReviveHeld[0]=false;ff_update_coop(.1f);assert(reviveProgress[1]==0);
    ffTestReviveHeld[0]=true;ff_update_coop(3.1f);assert(combatantLife[1]==LIFE_ALIVE && players[0].matter==0 && players[1].matter==35 && reviveComplete[1]>0);
}
static void bombs(void){
    empty();enemy_lab_spawn(4,ENEMY_TYPE_SAPPER,(Vector3){0,2,0});
    add_static_box_at_grid(0,0,0,2,3,3,GRAY);init_static_hash();
    int support=table_get_static_only(0,2,3);assert(support>=0);
    assert(enemy_create_charge(4,(Vector3){0,2,1},voxels[support].identity));
    players[4].isExposed=true;kill_player(4,0,true,false);int count=0,drop=-1;
    for(int n=0;n<ENEMY_CHARGE_MAX;n++)if(enemyCharges[n].active){count++;assert(enemyCharges[n].fuse==4);if(!enemyCharges[n].support)drop=n;}
    assert(count==2 && drop>=0);uint64_t id=enemyCharges[drop].id;
    enemy_lab_spawn(4,ENEMY_TYPE_STANDARD,(Vector3){25,2,25});enemy_tick(.5f);
    assert(enemyCharges[drop].active && enemyCharges[drop].id==id && enemyCharges[drop].fuse==3.5f);
    enemy_detach(enemy_voxel(id));enemy_tick(.1f);assert(!enemyCharges[drop].active);
    enemy_cleanup(true);for(int n=0;n<ENEMY_CHARGE_MAX;n++)assert(!enemyCharges[n].active);
}
static void lobbers(void){
    empty();players[0].pos=(Vector3){0,2,-35};players[1].respawn_timer=9999;
    enemy_lab_spawn(4,ENEMY_TYPE_LOBBER,(Vector3){0,2,35});
    add_static_box_at_grid(-5,5,0,8,0,0,GRAY);init_static_hash();ai_tick(4,.1f);
    assert(aiCombat[4].target==0 && aiCombat[4].visible);
    clear_world_voxels();init_static_hash();
    Vector3 velocity;assert(enemy_lob_arc((Vector3){0,3,35},players[0].pos,(Vector3){1,1,1},&velocity));assert(v_isfinite(velocity));
    add_static_box_at_grid(0,7,0,3,0,3,GRAY);init_static_hash();rebuild_voxel_hash();
    int seed=table_get_static_only(2,1,1);assert(seed>=0);
    double begin=GetTime();int root=enemy_extract_payload(4,seed);double ms=(GetTime()-begin)*1000;
    assert(root>=0);int count;Vector3 center,extent;assert(enemy_payload(4,root,&count,&center,&extent));
    assert(count==64 && extent.x<=1.01f && extent.y<=1.01f && extent.z<=1.01f);
    assert(table_get_static_only(7,3,3)>=0 && table_get_static_only(2,1,1)<0);
    players[4].pos=v_add(center,(Vector3){0,0,-4});players[4].goldTetherCharged=true;acquire_tether_hit(4,root,(Vector3){0,0,1});assert(players[4].tetherHolding);
    float initialY=center.y;
    for(int frame=0;frame<120;frame++){prepare_tether_forces();simulate_voxel_pbd_steps(1.f/60,2);}
    assert(resolve_tether_voxel(&players[4]));
    bool valid=enemy_payload(4,players[4].tetherVoxel,&count,&center,&extent);
    fprintf(stderr,"Lift valid=%d group=%d count=%d extent=%.2f %.2f %.2f centerY=%.2f\n",valid,build_glue_cluster_indices(players[4].tetherVoxel,glueClusterIndices),count,extent.x,extent.y,extent.z,center.y);
    assert(valid);
    assert(count==64 && center.y>initialY+.5f);
    for(int v=0;v<voxel_count;v++)assert(v_isfinite(voxels[v].pos));
    fprintf(stderr,"64-voxel extraction %.2f ms\n",ms);
    empty();players[0].pos=(Vector3){0,2,-30};players[1].respawn_timer=9999;
    enemy_lab_spawn(4,ENEMY_TYPE_LOBBER,(Vector3){1,2,-3});
    add_static_box_at_grid(0,7,0,3,0,3,GRAY);init_static_hash();rebuild_voxel_hash();
    enemySpecial[4].cooldown=0;enemy_update(4,.1f);assert(players[4].tetherHolding && enemySpecial[4].state==3);
    for(int frame=0;frame<120;frame++){
        enemy_update(4,1.f/60);prepare_tether_forces();simulate_voxel_pbd_steps(1.f/60,2);
    }
    assert(enemySpecial[4].actions>0 && !players[4].tetherHolding);
    bool launched=false;for(int n=0;n<MAX_COMBATANTS*2;n++)if(enemyLobs[n].id){assert(enemyLobs[n].count>=32);launched=true;}
    assert(launched);
}
static void feedback(void){
    empty();players[0].pos=(Vector3){0,2,0};players[4].pos=(Vector3){-5,2,0};
    players[0].invuln_timer=1;apply_matter_damage(0,4,10);assert(damageCues[0][0].until==0);
    players[0].invuln_timer=0;apply_matter_damage(0,4,10);assert(damageCues[0][0].known && damageCues[0][0].origin.x==-5);
    Camera3D camera={{0,2,0},{0,2,-1},{0,1,0},60,CAMERA_PERSPECTIVE};
    for(int n=0;n<8;n++){float a=n*PI/4;Vector2 d=combat_ui_direction((Vector3){sinf(a),0,-cosf(a)},camera);assert(fabsf(d.x-sinf(a))<.001f && fabsf(d.y+cosf(a))<.001f);}
    uint8_t packet[128];NetWriter w;NetReader r;netTransport.role=NET_ROLE_CLIENT;
    net_writer_init(&w,packet,sizeof(packet));net_write_u8(&w,1);net_write_u8(&w,0);net_write_u8(&w,1);net_write_f32(&w,0);net_write_f32(&w,2);net_write_f32(&w,6);
    net_reader_init(&r,packet,w.length);assert(combat_ui_receive(NET_MSG_COMBAT_UI,&r,1));
    bool rear=false;for(int n=0;n<8;n++)rear|=damageCues[0][n].origin.z==6;assert(rear);
    netTransport.role=NET_ROLE_OFFLINE;
}
static void preview(int w,int h,const char *path,bool fallback){
    empty();character_renderer_init();players[0].pos=(Vector3){-1,BASE_EYE_HEIGHT,-2};players[1].pos=(Vector3){1,BASE_EYE_HEIGHT,0};players[1].death_pos=players[1].pos;
    combatantLife[1]=LIFE_DOWNED;players[1].respawn_timer=9999;reviveOwner[1]=0;reviveProgress[1]=1.8f;
    grenadeInventory[0]=2;holograms[0].cooldown=0;grenadeCooldown[0]=0;
    Camera3D camera={{4,3,-6},{0,1,0},{0,1,0},60,CAMERA_PERSPECTIVE};
    combat_ui_add_damage(0,(Vector3){-8,2,-8},true);
    if(fallback){combat_ui_shutdown();reviveShaderAttempted=true;}
    RenderTexture2D rt=LoadRenderTexture(w,h);BeginTextureMode(rt);ClearBackground((Color){22,28,36,255});BeginMode3D(camera);
    DrawPlane((Vector3){0},(Vector2){20,20},DARKGRAY);weaponFrameSampled=false;drawViewPlayerIndex=-1;draw_rat_world(0,false);draw_rat_world(1,false);combat_ui_world();EndMode3D();
    combat_ui_hud(0,camera,w,h);combat_ui_abilities(0,w,h);EndTextureMode();
    if(!fallback)assert(reviveShader.id && reviveShader.id!=rlGetShaderIdDefault());
    Image image=LoadImageFromTexture(rt.texture);ImageFlipVertical(&image);ExportImage(image,path);UnloadImage(image);UnloadRenderTexture(rt);character_renderer_shutdown();combat_ui_shutdown();
}
int main(int argc,char **argv){
    assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,540,"Combat feedback");init_pbd_thread_pool();assert(initialize_physics_backend());
    free_revive();bombs();lobbers();feedback();
    preview(960,540,"artifacts/combat-feedback.png",false);preview(480,270,"artifacts/combat-feedback-split.png",false);preview(480,270,"artifacts/combat-feedback-fallback.png",true);
    clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Combat feedback checks passed");return 0;
}

