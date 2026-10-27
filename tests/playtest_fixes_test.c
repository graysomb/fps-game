#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static void fresh(void){
    enemy_cleanup(true);clear_all_chunks();clear_world_voxels();init_static_hash();worldPeriodic=false;
    memset(players,0,sizeof(players));memset(holograms,0,sizeof(holograms));memset(aiCombat,0,sizeof(aiCombat));
    gameMode=GAME_MODE_DEATHMATCH;gameState=GAME_STATE_PLAYING;activePlayers=2;enemyLab=false;creativeModeActive=false;
    netTransport.role=NET_ROLE_OFFLINE;worldFloorHalfSize=40;
    for(int i=0;i<MAX_COMBATANTS;i++){players[i].pos=(Vector3){20+i,2,20};players[i].matter=100;players[i].matterMax=100;players[i].kd_ratio=1;combatantLife[i]=LIFE_ALIVE;combatantTeam[i]=i;playerInput[i]=i<2?INPUT_TYPE_KEYBOARD:INPUT_TYPE_BOT_HARD;}
}
static bool held_compaction(bool old){
    fresh();players[0].pos=(Vector3){0,2,.5f};
    int root=addVoxel(0,2,-.5f,false,true,GRAY,0);
    int dummy=addVoxel(10,2,0,false,true,GRAY,0);
    int member=addVoxel(0,2,0,false,true,GRAY,0);
    rebuild_voxel_hash();glue_neighbor_faces_for_voxel(root);glue_neighbor_faces_for_voxel(member);
    players[0].tetherHolding=true;players[0].tetherVoxel=root;players[0].tetherVoxelIdentity=voxels[root].identity;
    prepare_tether_forces();assert(tetherHoldClusterTag[member]==1);
    uint64_t id=voxels[member].identity;remove_voxel_index(dummy);
    if(old)tetherHoldClusterTag[dummy]=0;
    handle_pbd_projectile_hits();
    bool ok=enemy_voxel(id)>=0 && players[0].matter==100;
    fprintf(stderr,"Held cluster after compaction: %s, matter %.0f\n",ok?"protected":"SELF HIT / DELETED",players[0].matter);
    return ok;
}
static void armor_and_release(void){
    fresh();players[1].enemyType=ENEMY_TYPE_GOLIATH;
    int armor=addVoxel(0,2,0,false,true,GRAY,0);uint64_t id=voxels[armor].identity;
    goliathStates[1].count=1;goliathStates[1].voxels[0].voxelIdentity=id;
    int bullet=addVoxel(0,2,1,false,true,WHITE,0);voxels[bullet].isBullet=true;voxels[bullet].owner=0;voxels[bullet].vel=(Vector3){0,0,-20};
    rebuild_voxel_hash();update_projectiles(.1f);
    assert(enemy_voxel(id)>=0 && get_goliath_owner_of_voxel(enemy_voxel(id))==1 && goliathStates[1].count==1);
    goliathStates[1].count=0;
    fresh();players[0].pos=(Vector3){0,2,0};int v=addVoxel(0,2,0,false,true,GRAY,0);
    voxels[v].tetherReleaseOverlap=1;handle_pbd_projectile_hits();assert(voxel_count==1 && players[0].matter==100);
    voxels[v].pos.x=5;handle_pbd_projectile_hits();assert(!voxels[v].tetherReleaseOverlap);
    voxels[v].pos.x=0;handle_pbd_projectile_hits();assert(players[0].matter<100);
    puts("Attached armor and initial release overlap passed");
}
static void upgrades(void){
    fresh();gameMode=GAME_MODE_FIREFIGHT;firefightHumans=1;
    memset(teamUpgradeLevels,0,sizeof(teamUpgradeLevels));assert(ff_build_cost(0)==10);
    teamUpgradeLevels[0]=1;teamUpgradeLevels[1]=2;teamUpgradeLevels[2]=3;
    assert(ff_build_cost(0)==5 && ff_harvest_amount(0)==20 && ff_throw_multiplier(0)==1.75f);
    firefightWave=2;firefightEndless=true;firefightFoundry=false;firefightPractice=false;ff_wave_complete();
    assert(teamUpgradeLevels[0]==1 && teamUpgradeLevels[1]==2 && teamUpgradeLevels[2]==3);
    teamUpgradeLevels[0]=9;teamUpgradeLevels[1]=teamUpgradeLevels[2]=8;
    assert(ff_upgrades_maxed() && ff_upgrade_choice(0,1)==-1 && ff_build_cost(0)==1);
    reset_firefight_match();assert(teamUpgradeLevels[0]==0 && teamUpgradeLevels[1]==0 && teamUpgradeLevels[2]==0);
}
static void decoys(void){
    // Aim at a floor point four meters ahead; stop precisely, then remain idle.
    fresh();players[0].pos=(Vector3){0,BASE_EYE_HEIGHT,0};
    players[0].pitch=-atan2f(BASE_EYE_HEIGHT,4)*RAD2DEG;
    assert(hologram_deploy(0));assert(fabsf(holograms[0].travelRemaining-4)<.001f);
    players[0].yaw=90; // Turning after deployment must not steer the decoy.
    for(int frame=0;frame<180;frame++)hologram_tick(1.f/60);
    assert(fabsf(holograms[0].body.pos.z+4)<.001f && fabsf(holograms[0].body.pos.x)<.001f);
    assert(holograms[0].travelRemaining==0 && holograms[0].body.vel.z==0 && holograms[0].life>0);
    Vector3 stopped=holograms[0].body.pos;
    hologram_tick(.25f);assert(Vector3Equals(stopped,holograms[0].body.pos));
    fresh();players[0].pos=(Vector3){0,BASE_EYE_HEIGHT,0};players[0].pitch=45;
    assert(hologram_deploy(0));assert(holograms[0].travelRemaining>24 && holograms[0].travelRemaining<=hologram_max_range());
    // Cross nearly the entire map, then idle for the remaining eight seconds.
    fresh();players[0].pos=(Vector3){0,BASE_EYE_HEIGHT,35};
    players[0].pitch=-atan2f(BASE_EYE_HEIGHT,70)*RAD2DEG;
    assert(hologram_deploy(0));assert(fabsf(holograms[0].travelRemaining-70)<.01f);
    assert(holograms[0].life>14 && holograms[0].life<=hologram_max_life());
    for(int frame=0;frame<15*60;frame++)hologram_tick(1.f/60);
    assert(fabsf(holograms[0].body.pos.z+35)<.01f && holograms[0].travelRemaining==0);
    assert(holograms[0].life>6 && holograms[0].life<=8);
    fresh();activePlayers=5;players[0].pos=(Vector3){0,BASE_EYE_HEIGHT,0};
    assert(hologram_deploy(0));holograms[0].direction=(Vector3){0};
    players[4].enemyType=ENEMY_TYPE_SWARMER;players[4].pos=players[0].pos;
    for(int hit=0;hit<5;hit++){
        players[4].contactDamageTimer=0;
        hologram_tick(.01f);
        assert(players[4].respawn_timer==0 && combatantLife[4]==LIFE_ALIVE);
        assert(holograms[0].health==40-hit*10 && holograms[0].damageFlicker>0);
        assert(players[4].contactDamageTimer==.5f);
        hologram_tick(.01f);assert(holograms[0].health==40-hit*10);
    }
    assert(holograms[0].life==0);
    fresh();players[0].pos=(Vector3){0,BASE_EYE_HEIGHT,0};players[0].yaw=0;
    assert(hologram_deploy(0));assert(!hologram_deploy(0));
    for(int i=0;i<60;i++)hologram_tick(1.f/60);
    assert(holograms[0].body.pos.z<-1 && holograms[0].body.onGround);
    hologram_damage(MAX_COMBATANTS,49);assert(holograms[0].life>0 && holograms[0].health==1);
    hologram_damage(MAX_COMBATANTS,1);assert(!holograms[0].life && players[0].matter==100);
    fresh();players[0].pos=(Vector3){0,BASE_EYE_HEIGHT,0};players[0].yaw=0;
    add_static_box_at_grid(-4,4,0,5,-5,-5,GRAY);rebuild_voxel_hash();assert(hologram_deploy(0));
    for(int i=0;i<180;i++)hologram_tick(1.f/60);
    assert(holograms[0].body.pos.z>-2 && v_length(holograms[0].direction)==0);
    fresh();activePlayers=5;players[0].pos=(Vector3){0,BASE_EYE_HEIGHT,-4};players[4].pos=(Vector3){0,BASE_EYE_HEIGHT,0};players[4].yaw=0;
    assert(hologram_deploy(0));players[0].pos.x=30;players[1].pos.x=30;
    ai_tick(4,.1f);assert(aiCombat[4].target==MAX_COMBATANTS && aiCombat[4].visible);
    for(int type=ENEMY_TYPE_SCRAPPER;type<ENEMY_TYPE_COUNT;type++){players[4].enemyType=type;enemySpecial[4]=(EnemySpecialState){0};enemy_update(4,.1f);}
    hologram_tick(.1f);assert(holograms[0].life>0);
    players[0].respawn_timer=3;hologram_tick(.1f);assert(!holograms[0].life);
    fresh();players[0].pos=(Vector3){0,2,-3};assert(hologram_deploy(0));
    int bullet=addVoxel(0,2,0,false,true,WHITE,0);voxels[bullet].owner=1;voxels[bullet].isBullet=true;
    hologram_projectile(bullet,(Vector3){0,2,0},(Vector3){0,2,-6},10);
    assert(holograms[0].health==40);hologram_projectile(bullet,(Vector3){0,2,0},(Vector3){0,2,-6},10);assert(holograms[0].health==40);
    assert(voxel_count==1); /* A decoy never consumes a projectile. */
    hologram_blast((Vector3){0,2,-3},4,60);assert(!holograms[0].life);
    puts("Hologram movement, damage, targeting and lifecycle passed");
}
static void appearance(void){
    fresh();character_renderer_init();weapon_renderer_init();players[0].pos=(Vector3){1,BASE_EYE_HEIGHT,0};players[0].yaw=0;players[0].onGround=true;
    assert(hologram_deploy(0));holograms[0].direction=(Vector3){0};hologram_tick(.01f);
    players[0].pos.x=-1;weaponFrameSampled=false;
    Camera3D camera={{3,2.8f,-6},{0,.7f,0},{0,1,0},45,CAMERA_PERSPECTIVE};
    RenderTexture2D image=LoadRenderTexture(960,540);BeginTextureMode(image);ClearBackground((Color){20,25,35,255});BeginMode3D(camera);
    DrawPlane((Vector3){0,0,0},(Vector2){20,20},DARKGRAY);draw_rat_world(0,false);draw_player_weapon(0,camera.position,false,WEAPON_COLOR);character_draw_visor(characterSamples[0]);hologram_draw(false);
    EndMode3D();EndTextureMode();Image capture=LoadImageFromTexture(image.texture);ImageFlipVertical(&capture);ExportImage(capture,"artifacts/player-and-hologram.png");UnloadImage(capture);UnloadRenderTexture(image);
    weapon_renderer_shutdown();character_renderer_shutdown();
}
static void trap_map(bool crosscut){
    fresh();activePlayers=0;if(crosscut)buildCrosscut();else buildFoundryHoldout();
    assert(foundryTrapCount==(crosscut?13:4));rebuild_voxel_hash();rebuild_static_hash_if_dirty();
    activePlayers=3;players[0].pos=(Vector3){38,2,38};players[2].pos=(Vector3){-38,2,38};
    playerInput[1]=INPUT_TYPE_BOT_HARD;players[1].enemyType=ENEMY_TYPE_STANDARD;
    for(int round=0;round<2;round++){
    for(int n=0;n<foundryTrapCount;n++){
        FoundryTrap *t=&foundryTraps[n];players[1].pos=(Vector3){(t->x0+t->x1+1)*.25f,1,(t->z0+t->z1+1)*.25f};
        if(t->kind==FOUNDRY_TIP)players[1].pos.z=(t->z1+1)*.5f+1;
        players[1].matter=100;players[1].respawn_timer=players[1].invuln_timer=0;combatantLife[1]=LIFE_ALIVE;
        assert(foundry_trap_hit(table_get_static_only(t->rx,t->ry,t->rz),0,false));foundry_traps_prepare();assert(t->active);
        for(int frame=0;frame<180 && players[1].respawn_timer==0;frame++){
            foundry_traps_prepare();simulate_voxel_pbd_steps(1.f/60,1);foundry_traps_hits();
        }
        fprintf(stderr,"%s trap %d: %s\n",crosscut?"Crosscut":"Foundry",n,players[1].respawn_timer>0?"KILL":"MISS");
        if(players[1].respawn_timer<=0){for(int v=0;v<voxel_count;v++)if(voxels[v].foundryTrapId==n+1){fprintf(stderr,"payload pos %.2f %.2f %.2f speed %.2f active %d heavy %d\n",voxels[v].pos.x,voxels[v].pos.y,voxels[v].pos.z,v_length(voxels[v].vel),voxels[v].vgs_active,voxels[v].foundryTrapHeavy);break;}}
        assert(players[1].respawn_timer>0);assert(players[2].matter==100 && players[2].respawn_timer==0);
        for(int v=voxel_count-1;v>=0;v--)if(voxels[v].foundryTrapId==n+1)remove_voxel_index(v);
    }
    if(!round){
        recycleFrameCounter+=RECYCLE_STATIC_RESTORE_DELAY+1;
        for(int step=0;step<4000 && recycleQueueCount;step++){recycleFrameCounter++;regeneration_step(1.f/60);}
        foundry_traps_prepare();
        for(int n=0;n<foundryTrapCount;n++)assert(foundry_trap_complete(&foundryTraps[n]) && !foundryTraps[n].active);
        fprintf(stderr,"%s all traps restored; repeating physical kills\n",crosscut?"Crosscut":"Foundry");
    }
    }
}
int main(int argc,char **argv){
    assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);
    SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,540,"Playtest fixes");init_pbd_thread_pool();assert(initialize_physics_backend());
    bool old=false;for(int i=1;i<argc;i++)if(!strcmp(argv[i],"--repro-old"))old=true;
    if(old){bool ok=held_compaction(true);return ok?1:0;}
    assert(held_compaction(false));armor_and_release();upgrades();decoys();appearance();trap_map(true);trap_map(false);
    if(physicsBackend.requested==PHYSICS_BACKEND_GPU_GL43)assert(physics_backend_is_gpu(physicsBackend.active));
    fprintf(stderr,"Validated backend: %s\n",physics_backend_name(physicsBackend.active));
    clear_all_chunks();clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Playtest fixes passed");return 0;
}

