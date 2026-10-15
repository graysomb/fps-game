#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
static void fresh(void) {
    clear_all_chunks();clear_world_voxels();activePlayers=0;worldPeriodic=false;worldFloorHalfSize=40;
    gameMode=GAME_MODE_FIREFIGHT;firefightHumans=1;firefightPractice=false;firefightFoundry=true;
}
static void release_checks(void) {
    fresh();buildFoundryHoldout();assert(foundryTrapCount==6);
    fprintf(stderr,"Foundry cells=%d, fragile=%d vs %d (%.1f%% reduction)\n",voxel_count,foundryFragileCells,foundryFragileReferenceCells,100.f*(1-(float)foundryFragileCells/foundryFragileReferenceCells));
    assert(foundryFragileCells*2<foundryFragileReferenceCells);
    for(int x=-18;x<18;x++)for(int y=0;y<4;y++)for(int z=-14;z<14;z++)assert(table_get_static_only(x,y,z)>=0);
    for(int n=0;n<6;n++) {
        FoundryTrap *t=&foundryTraps[n];assert(foundry_trap_complete(t));
        assert(foundry_trap_hit(table_get_static_only(t->rx,t->ry,t->rz),0,false));
    }
    double start=GetTime();foundry_traps_prepare();fprintf(stderr,"Six-trap release %.2f ms\n",1000*(GetTime()-start));
    for(int n=0;n<6;n++) {
        FoundryTrap *t=&foundryTraps[n];int count=0;
        for(int i=0;i<voxel_count;i++)if(voxels[i].foundryTrapId==n+1){count++;assert(voxels[i].restorationQueued);assert(voxels[i].owner==0);}
        assert(count==(t->x1-t->x0+1)*(t->y1-t->y0+1)*(t->z1-t->z0+1));assert(count<=256);
        assert(table_get_static_only(t->rx,t->ry,t->rz)<0);
    }
    float before=0,after=0;int count=0;
    for(int i=0;i<voxel_count;i++)if(voxels[i].foundryTrapId==1){before+=voxels[i].pos.y;count++;}
    double worst=0;
    for(int f=0;f<45;f++) {start=GetTime();simulate_voxel_pbd_steps(1.f/60,1);foundry_traps_hits();worst=fmax(worst,1000*(GetTime()-start));}
    for(int i=0;i<voxel_count;i++)if(voxels[i].foundryTrapId==1){after+=voxels[i].pos.y;assert(isfinite(voxels[i].pos.y));}
    fprintf(stderr,"Drop %.3f -> %.3f, worst simulation %.2f ms backend=%s\n",before/count,after/count,worst,physics_backend_name(physicsBackend.active));assert(after<before-count*.2f);
    // No rearm while the payload remains; then normal chunk regeneration.
    foundry_traps_prepare();assert(foundryTraps[0].active);
    for(int i=voxel_count-1;i>=0;i--)if(voxels[i].foundryTrapId)remove_voxel_index(i);
    recycleFrameCounter+=RECYCLE_STATIC_RESTORE_DELAY+1;
    for(int f=0;f<3000 && recycleQueueCount;f++){recycleFrameCounter++;regeneration_step(1.f/60);}
    foundry_traps_prepare();for(int n=0;n<6;n++){assert(foundry_trap_complete(&foundryTraps[n]));assert(!foundryTraps[n].active);}
    assert(foundry_trap_hit(table_get_static_only(foundryTraps[0].rx,foundryTraps[0].ry,foundryTraps[0].rz),0,false));foundry_traps_prepare();assert(foundryTraps[0].activation==2);
}
static void damage_checks(void) {
    fresh();foundry_trap_define(FOUNDRY_DROP,0,3,8,8,0,3,-1,8,0);
    foundry_trap_hit(table_get_static_only(-1,8,0),0,false);foundry_traps_prepare();
    activePlayers=2;memset(players,0,sizeof(players));
    for(int j=0;j<2;j++){players[j].pos=(Vector3){.75f,2,.75f};players[j].matter=100;combatantLife[j]=LIFE_ALIVE;combatantTeam[j]=j;}
    players[0].pos.x=10;players[1].enemyType=ENEMY_TYPE_STANDARD;playerInput[1]=INPUT_TYPE_BOT_HARD;
    foundry_traps_hits();assert(players[1].respawn_timer==0);
    for(int i=0;i<voxel_count;i++)if(voxels[i].foundryTrapId){voxels[i].foundryTrapPrevious=voxels[i].pos;voxels[i].pos.y=.25f;voxels[i].vel.y=-5;}
    foundry_traps_hits();assert(players[1].respawn_timer>0);int kills=players[0].debrisKills;assert(kills==1);
    foundry_traps_hits();assert(players[0].debrisKills==kills);
    // Invulnerability and a clear miss do not consume the activation hit mask.
    players[1].respawn_timer=0;players[1].invuln_timer=1;foundryTraps[0].victims=0;
    foundry_traps_hits();assert(players[1].respawn_timer==0 && !foundryTraps[0].victims);
    players[1].invuln_timer=0;players[1].pos.x=12;
    foundry_traps_hits();assert(players[1].respawn_timer==0 && !foundryTraps[0].victims);
    players[1].pos.x=.75f;firefightHumans=2;playerInput[1]=INPUT_TYPE_GAMEPAD;combatantTeam[1]=0;
    foundry_traps_hits();assert(combatantLife[1]==LIFE_DOWNED);
    // A protected Walking Foundry cannot be bypassed; exposed core takes one hit.
    players[1].respawn_timer=0;players[1].enemyType=ENEMY_TYPE_FOUNDRY;foundryTraps[0].victims=0;
    enemySpecial[1]=(EnemySpecialState){.coreHealth=300};foundry_traps_hits();assert(enemySpecial[1].coreHealth==300);
    enemySpecial[1].coreOpen=true;foundry_traps_hits();assert(enemySpecial[1].coreHealth==200);
    foundry_traps_hits();assert(enemySpecial[1].coreHealth==200);
    puts("Trap stationary safety, swept impact, kill credit and dedup passed");
}
static void trigger_checks(void) {
    fresh();foundry_trap_define(FOUNDRY_DROP,0,3,8,8,0,3,-1,8,0);
    int shot=addVoxel(-2,4.25f,.25f,false,true,WHITE,0);voxels[shot].isBullet=true;voxels[shot].owner=0;voxels[shot].vel=(Vector3){30,0,0};
    update_projectiles(.1f);assert(foundryTraps[0].pending);foundry_traps_prepare();assert(foundryTraps[0].active);
    fresh();foundry_trap_define(FOUNDRY_TIP,0,3,1,6,0,3,0,0,0);
    enemy_blast_terrain((EnemyCharge){.pos={0,.25f,0},.owner=0});foundry_traps_prepare();assert(foundryTraps[0].active && foundryTraps[0].owner==0);
    assert(recycleQueueCount>=foundryTraps[0].originalCount);
    for(int i=voxel_count-1;i>=0;i--)if(voxels[i].simulate)remove_voxel_index(i);
    recycleFrameCounter+=RECYCLE_STATIC_RESTORE_DELAY+1;
    for(int f=0;f<3000 && recycleQueueCount;f++){recycleFrameCounter++;regeneration_step(1.f/60);}
    foundry_traps_prepare();assert(foundry_trap_complete(&foundryTraps[0]) && !foundryTraps[0].active);
    puts("Projectile release and blast-destroyed assembly regeneration passed");
}
static void capture(void) {
    fresh();buildFoundryHoldout();voxel_ca_init();
    Camera3D cameras[3]={{{29,22,29},{0,0,0},{0,1,0},65,CAMERA_PERSPECTIVE},{{28,5,9},{21,2,0},{0,1,0},65,CAMERA_PERSPECTIVE},{{13,4,14},{21,2,21},{0,1,0},65,CAMERA_PERSPECTIVE}};
    RenderTexture2D target=LoadRenderTexture(960,540);
    for(int n=0;n<3;n++) {
        voxel_render_frame++;BeginTextureMode(target);ClearBackground((Color){30,37,46,255});BeginMode3D(cameras[n]);
        DrawPlane((Vector3){0,-.02f,0},(Vector2){80,80},(Color){55,61,68,255});DrawVoxels(cameras[n]);foundry_traps_draw();EndMode3D();EndTextureMode();
        Image im=LoadImageFromTexture(target.texture);ImageFlipVertical(&im);ExportImage(im,TextFormat("artifacts/foundry-traps-%d.png",n));UnloadImage(im);
    }
    UnloadRenderTexture(target);voxel_ca_shutdown();
}
int main(int argc,char **argv) {
    assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,540,"Foundry traps");init_pbd_thread_pool();assert(initialize_physics_backend());
    release_checks();damage_checks();trigger_checks();capture();clear_all_chunks();clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Foundry traps passed");
}

