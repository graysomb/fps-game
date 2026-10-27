#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
static void fresh(void) {
    clear_all_chunks();clear_world_voxels();activePlayers=0;worldPeriodic=false;worldFloorHalfSize=40;
    gameMode=GAME_MODE_FIREFIGHT;firefightHumans=1;firefightPractice=false;firefightFoundry=true;
}
static void release_checks(void) {
    fresh();buildFoundryHoldout();assert(foundryTrapCount==4);
    fprintf(stderr,"Foundry cells=%d, fragile=%d vs %d (%.1f%% reduction)\n",voxel_count,foundryFragileCells,foundryFragileReferenceCells,100.f*(1-(float)foundryFragileCells/foundryFragileReferenceCells));
    assert(foundryFragileCells*2<foundryFragileReferenceCells);
    for(int x=-18;x<18;x++)for(int y=0;y<4;y++)for(int z=-14;z<14;z++)assert(table_get_static_only(x,y,z)>=0);
    for(int n=0;n<foundryTrapCount;n++) {
        FoundryTrap *t=&foundryTraps[n];assert(foundry_trap_complete(t));
        assert(foundry_trap_hit(table_get_static_only(t->rx,t->ry,t->rz),0,false));
    }
    double start=GetTime();foundry_traps_prepare();fprintf(stderr,"All-trap release %.2f ms\n",1000*(GetTime()-start));
    for(int n=0;n<foundryTrapCount;n++) {
        FoundryTrap *t=&foundryTraps[n];int count=0;
        for(int i=0;i<voxel_count;i++)if(voxels[i].foundryTrapId==n+1){count++;assert(voxels[i].restorationQueued);assert(voxels[i].owner==0);}
        int expected=0;for(int k=0;k<t->originalCount;k++){Voxel *v=&t->original[k];expected+=foundry_trap_contains(t,v->gx,v->gy,v->gz);}
        assert(count==expected);assert(count<=FOUNDRY_TRAP_PAYLOAD);
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
    foundry_traps_prepare();for(int n=0;n<foundryTrapCount;n++){assert(foundry_trap_complete(&foundryTraps[n]));assert(!foundryTraps[n].active);}
    assert(foundry_trap_hit(table_get_static_only(foundryTraps[0].rx,foundryTraps[0].ry,foundryTraps[0].rz),0,false));foundry_traps_prepare();assert(foundryTraps[0].activation==2);
}
static void post_checks(void) {
    for(int corner=0;corner<4;corner++) {
        fresh();buildFoundryHoldout();FoundryTrap *t=&foundryTraps[0];assert(t->fourPosts);
        int x=corner&1?t->x1:t->x0,z=corner&2?t->z1:t->z0;
        if(corner==3){
            int id=table_get_static_only(x,t->ry+1,z);Vector3 at=voxels[id].pos;
            remove_voxel_index(id);assert(addVoxel(at.x,at.y,at.z,false,true,ORANGE,0)>=0);
        }

        assert(foundry_trap_hit(table_get_static_only(x,t->ry,z),0,false));foundry_traps_prepare();
        for(int a=0;a<4;a++){x=a&1?t->x1:t->x0;z=a&2?t->z1:t->z0;
            assert(table_get_static_only(x,t->ry,z)<0);
            assert(table_get_static_only(x,t->ry-1,z)>=0);
            assert(table_get_static_only(x,t->ry+1,z)<0);
        }
        for(int v=0;v<voxel_count;v++)if(voxels[v].foundryTrapId==1){
            assert(voxels[v].orig_min_gy>=t->y0);
            assert(ColorToInt(voxels[v].color)!=ColorToInt(ORANGE));
        }
        for(int v=0;v<voxel_count;v++){
            Voxel *cell=&voxels[v];
            assert(!(cell->simulate && !cell->isBullet && cell->orig_min_gy>=t->ry && cell->orig_min_gy<t->y0 &&
                (cell->orig_min_gx==t->x0 || cell->orig_min_gx==t->x1) && (cell->orig_min_gz==t->z0 || cell->orig_min_gz==t->z1)));
        }
        assert(!foundryTraps[1].active);
        if(corner==0){
            activePlayers=1;playerInput[0]=INPUT_TYPE_KEYBOARD;combatantTeam[0]=0;combatantLife[0]=LIFE_ALIVE;
            players[0]=(Player){.pos={(t->x0+t->x1+1)*.25f,1.5f,(t->z0+t->z1+1)*.25f},.matter=100};
            for(int frame=0;frame<90 && players[0].respawn_timer<=0;frame++){simulate_voxel_pbd_steps(1.f/60,1);foundry_traps_hits();}
            assert(players[0].respawn_timer>0);
        }
    }
    puts("Each corner trigger releases all four posts; other slab stays supported");
}
/* Seed the public random map construction, without a production test override. */
static unsigned selection_seed(int corner,bool pillar) {
    for(unsigned seed=1;seed<10000;seed++){
        SetRandomSeed(seed);int bay=GetRandomValue(0,3),post=GetRandomValue(0,3);
        if((pillar?post:bay)==corner)return seed;
    }
    assert(false);return 1;
}
static void bay_checks(void) {
    for(int corner=0;corner<4;corner++){
        int bay=2;unsigned seed=selection_seed(corner,false);SetRandomSeed(seed);
        fresh();buildFoundryHoldout();FoundryTrap *t=&foundryTraps[bay];assert(t->kind==FOUNDRY_BAY);
        int orange=0,roof=0,wall=0;unsigned signature=2166136261u;
        for(int k=0;k<t->originalCount;k++){
            Voxel *v=&t->original[k];bool crack=foundry_bay_crack(t,v->gx,v->gy,v->gz);
            signature=(signature^(unsigned)crack)*16777619u;
            if(crack){orange++;roof+=v->gy==11;wall+=v->gy<11;assert(v->gy>0);}
        }
        assert(roof>10 && wall>20 && orange<t->originalCount/2);
        /* Rebuilding the same map must reproduce the exact crack network. */
        fresh();SetRandomSeed(seed);buildFoundryHoldout();t=&foundryTraps[bay];unsigned rebuilt=2166136261u;
        for(int k=0;k<t->originalCount;k++){Voxel *v=&t->original[k];rebuilt=(rebuilt^(unsigned)foundry_bay_crack(t,v->gx,v->gy,v->gz))*16777619u;}
        assert(signature==rebuilt);
        int trigger=table_get_static_only(t->rx,t->ry,t->rz);assert(trigger>=0 && ColorToInt(voxels[trigger].color)==ColorToInt(ORANGE));
        assert(foundry_trap_hit(trigger,-1,false));foundry_traps_prepare();assert(t->active);
        int payload=0;
        for(int k=0;k<t->originalCount;k++){Voxel *v=&t->original[k];int id=table_get_static_only(v->gx,v->gy,v->gz);
            if(v->gy>=1)assert(id<0);else assert(id>=0);
            payload+=foundry_trap_contains(t,v->gx,v->gy,v->gz);
        }
        int actual=0;for(int v=0;v<voxel_count;v++)if(voxels[v].foundryTrapId==bay+1){actual++;assert(ColorToInt(voxels[v].color)!=ColorToInt(ORANGE));}
        assert(actual==payload);
        rebuild_glue_topology_neighbors();build_glue_cluster_ids(collisionVoxelClusterIds);
        int fragments=0;bool seen[MAX_VOXELS]={0};
        for(int v=0;v<voxel_count;v++)if(voxels[v].foundryTrapId==bay+1){
            int group=collisionVoxelClusterIds[v];assert(group>=0);
            if(!seen[group]){seen[group]=true;fragments++;}
        }
        fprintf(stderr,"Bay %d: %d orange crack cells, %d physical fragments\n",bay,orange,fragments);
        assert(fragments>=4);
        for(int n=0;n<foundryTrapCount;n++)if(n!=bay)assert(!foundryTraps[n].active && foundry_trap_complete(&foundryTraps[n]));
        activePlayers=1;playerInput[0]=INPUT_TYPE_KEYBOARD;combatantTeam[0]=0;combatantLife[0]=LIFE_ALIVE;
        players[0]=(Player){.pos={(t->x0+t->x1+1)*.25f,1.5f,(t->z0+t->z1+1)*.25f},.matter=100};
        for(int frame=0;frame<180 && players[0].respawn_timer<=0;frame++){simulate_voxel_pbd_steps(1.f/60,1);foundry_traps_hits();}
        fprintf(stderr,"Bay %d: %d dynamic cells, center crush=%d\n",bay,payload,players[0].respawn_timer>0);
        assert(players[0].respawn_timer>0);
    }
    puts("Four sparse factory bays: branching orange fractures, isolated activation and lethal roof collapse passed");
}
static void pillar_checks(void) {
    for(int corner=0;corner<4;corner++){
        int n=3;SetRandomSeed(selection_seed(corner,true));
        fresh();buildFoundryHoldout();FoundryTrap *t=&foundryTraps[n];
        assert(t->kind==FOUNDRY_PILLAR && t->originalCount==1280);
        assert(t->x1-t->x0+1==8 && t->z1-t->z0+1==8 && t->y1==19);
        assert(foundry_trap_hit(table_get_static_only(t->rx,t->ry,t->rz),-1,false));foundry_traps_prepare();
        int expected=0,actual=0;float before=0,after=0;
        for(int k=0;k<t->originalCount;k++){Voxel *v=&t->original[k];
            expected+=foundry_trap_contains(t,v->gx,v->gy,v->gz);
            assert((table_get_static_only(v->gx,v->gy,v->gz)>=0)==(v->gy==0));
        }
        rebuild_glue_topology_neighbors();build_glue_cluster_ids(collisionVoxelClusterIds);
        bool seen[MAX_VOXELS]={0};int fragments=0;
        for(int i=0;i<voxel_count;i++)if(voxels[i].foundryTrapId==n+1){
            actual++;before+=voxels[i].pos.y;assert(ColorToInt(voxels[i].color)!=ColorToInt(ORANGE));
            int g=collisionVoxelClusterIds[i];assert(g>=0);if(!seen[g]){seen[g]=true;fragments++;}
        }
        assert(actual==expected && fragments>=4);
        for(int j=0;j<foundryTrapCount;j++)if(j!=n)assert(!foundryTraps[j].active && foundry_trap_complete(&foundryTraps[j]));
        for(int f=0;f<60;f++)simulate_voxel_pbd_steps(1.f/60,1);
        for(int i=0;i<voxel_count;i++)if(voxels[i].foundryTrapId==n+1){assert(isfinite(voxels[i].pos.y));after+=voxels[i].pos.y;}
        assert(after<before-actual*.5f);
        fprintf(stderr,"Corner pillar %d: %d fragments, %d falling cells\n",n,fragments,actual);
    }
}
static void selection_checks(void) {
    unsigned bayMask=0,pillarMask=0;
    for(unsigned seed=1;seed<=16;seed++){
        fresh();SetRandomSeed(seed);buildFoundryHoldout();assert(foundryTrapCount==4);
        FoundryTrap *bay=&foundryTraps[2],*pillar=&foundryTraps[3];
        assert(bay->kind==FOUNDRY_BAY && pillar->kind==FOUNDRY_PILLAR);
        int b=(bay->x0>0?2:0)+(bay->z0>0?1:0),p=(pillar->x0>0?2:0)+(pillar->z0>0?1:0);
        bayMask|=1u<<b;pillarMask|=1u<<p;
        for(int corner=0;corner<4;corner++){
            int sx=corner&2?1:-1,sz=corner&1?1:-1;
            if(corner!=p)for(int z=sz*70-4;z<=sz*70+3;z++)for(int y=0;y<20;y++)for(int x=sx*70-4;x<=sx*70+3;x++){
                int id=table_get_static_only(x,y,z);assert(id>=0);
                assert(ColorToInt(voxels[id].color)!=ColorToInt(ORANGE));assert(!foundry_trap_hit(id,-1,false));
            }
            if(corner!=b){
                int x0=(int)floorf((sx*23-6)*1.875f),x1=(int)ceilf((sx*23+7)*1.875f)-1;
                int z0=(int)floorf((sz*22-5)*1.875f),z1=(int)ceilf((sz*22+6)*1.875f)-1;
                for(int z=z0;z<=z1;z++)for(int y=0;y<=11;y++)for(int x=x0;x<=x1;x++){
                    int id=table_get_static_only(x,y,z);if(id<0)continue;
                    assert(ColorToInt(voxels[id].color)!=ColorToInt(ORANGE));assert(!foundry_trap_hit(id,-1,false));
                }
            }
        }
    }
    assert(bayMask==15 && pillarMask==15);
    int previousBay=foundryCrackedBay,previousPillar=foundryCrackedPillar;
    fresh();fpRestoring=true;buildFoundryHoldout();fpRestoring=false;
    assert(foundryCrackedBay==previousBay && foundryCrackedPillar==previousPillar && foundryTrapCount==4);
    puts("Random selection: exactly one bay and one pillar, normal unselected structures, all corners covered");
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
    release_checks();post_checks();bay_checks();pillar_checks();selection_checks();damage_checks();trigger_checks();capture();clear_all_chunks();clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Foundry traps passed");
}

