#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
static void step(float dt){
    tower_tick(dt);eq_tick(dt);
    for(int i=4;i<activePlayers;i++)if(players[i].respawn_timer<=0){UpdateBot(i,dt);net_integrate_player(&players[i],dt,true);}
}
static EnemyStructure *built(void){for(int b=0;b<ENEMY_STRUCTURE_MAX;b++)if(enemyStructures[b].tower && enemyStructures[b].built)return &enemyStructures[b];return NULL;}
static void fixture(void){
    enemy_lab_start(ENEMY_TYPE_BUILDER);worldPeriodic=false;enemyLabGod=true;
    players[0].invuln_timer=1000;players[0].pos=(Vector3){0,1,-22};
}
static void build_and_climb(void){
    fixture();EnemyStructure *t=NULL;
    for(int f=0;f<3600;f++){
        step(1.f/60);t=built();
        if(f%300==0)fprintf(stderr,"frame %d builder %.2f %.2f rifle %.2f %.2f %.2f state %d wp %d tower %d reason %s\n",f,players[4].pos.x,players[4].pos.z,players[5].pos.x,players[5].pos.y,players[5].pos.z,towerRifles[5].state,towerRifles[5].waypoint,t!=NULL,enemySpecial[4].reason);
        if(t && tower_stationed(5))break;
    }
    assert(t && tower_stationed(5));assert(players[5].pos.y>4);
    fprintf(stderr,"Built %d cells; rifleman climbed to %.2fm\n",t->count,players[5].pos.y);
    assert(t->occupant==5);int baseline=voxel_count;
    players[4].matter=0;players[4].isExposed=true;
    int victim=enemy_voxel(t->ids[0]);assert(victim>=0);remove_voxel_index(victim);
    for(int f=0;f<30;f++)step(1.f/60);assert(voxel_count==baseline);
    /* Actual activation must retain tower ownership and not enqueue map restoration. */
    int k=t->count/2;victim=enemy_voxel(t->ids[k]);assert(victim>=0);
    static UnitVoxelBuffer buffer;buffer.count=0;Voxel v=voxels[victim];
    assert(unit_voxel_buffer_push(&buffer,v.gx,v.gy,v.gz,v.color,0,false,victim,0,0));
    remove_buffered_static_voxels(&buffer);emit_unit_voxels_from_units(&buffer,false,true,-1);
    int dynamic=0;for(int n=0;n<voxel_count;n++)if(voxels[n].towerTag==t->tag && voxels[n].simulate)dynamic++;
    assert(dynamic>0);for(int f=0;f<30;f++)step(1.f/60);
    assert(voxel_count==baseline);for(int n=0;n<voxel_count;n++)assert(!(voxels[n].towerTag==t->tag&&voxels[n].simulate));
    /* Fully reclaim/rebuild repeatedly, with bounded material count. */
    players[5].pos=tower_point(t,(Vector3){-2,1,8});tower_release(5);players[5].respawn_timer=9999;
    double repairStart=GetTime();
    for(int attempt=0;attempt<100;attempt++){
        for(int n=voxel_count-1;n>=0;n--)if(voxels[n].towerTag==t->tag)remove_voxel_index(n);
        for(int f=0;f<30;f++)step(1.f/60);assert(voxel_count==baseline);
    }
    fprintf(stderr,"100 full collapses/repairs: %.2f ms each, %d voxels retained\n",(GetTime()-repairStart)*10,voxel_count);
    victim=enemy_voxel(t->ids[0]);remove_voxel_index(victim);players[4].respawn_timer=9999;
    for(int f=0;f<30;f++)step(1.f/60);assert(voxel_count==baseline-1);
    puts("Construction, real stair navigation, rapid repairs, activated ownership, total collapse and owner death passed");
}
static EnemyStructure *ready_tower(void){
    fixture();tower_blueprint();
    EnemyStructure *t=&enemyStructures[0];*t=(EnemyStructure){.tower=true,.owner=4,.generation=enemySpecial[4].generation,.origin={0,0,5},.quarter=0,.occupant=-1,.tag=towerNextTag++,.construction=2};
    players[4].pos=tower_point(t,(Vector3){-1,1,6.5f});players[5].pos=tower_point(t,(Vector3){-2,1,2});assert(tower_place(t));return t;
}
static void garrison_checks(void){
    EnemyStructure *t=ready_tower();players[5].pos=tower_waypoint(t,9);players[5].onGround=true;t->occupant=5;towerRifles[5]=(TowerRifle){.tag=t->tag,.state=2};
    players[0].pos=tower_point(t,(Vector3){1.5f,1,-25});players[0].matter=100;
    ai_select(5,0,players[0].pos,AI_ALERT_SEEN);players[5].yaw=players[5].pitch=0;players[5].last_shot_time=-100;
    uint32_t before=weaponVisuals[5].shot_sequence;
    for(int f=0;f<300;f++)step(1.f/60);
    assert(tower_stationed(5) && aiCombat[5].visible);assert(weaponVisuals[5].shot_sequence>before);
    players[0].pos=tower_point(t,(Vector3){5,1,7.5f});ai_select(5,0,players[0].pos,AI_ALERT_SEEN);
    for(int f=0;f<120;f++)step(1.f/60);assert(tower_stationed(5));
    players[0].matter=0;players[0].isExposed=true;bool jumped=false,landed=false;
    for(int f=0;f<360;f++){step(1.f/60);jumped|=towerRifles[5].state==3;landed|=jumped&&players[5].pos.y<1.1f;if(landed)break;}
    fprintf(stderr,"Jump attack: jumped=%d landed=%d rifle=(%.2f %.2f %.2f) state=%d\n",jumped,landed,players[5].pos.x,players[5].pos.y,players[5].pos.z,towerRifles[5].state);
    assert(jumped && landed);
    players[0].matter=100;players[0].isExposed=false;
    for(int f=0;f<120;f++)step(1.f/60);assert(towerRifles[5].state!=3);
    const int guns[]={EQ_DEFAULT,EQ_GLOB,EQ_GLOCK,EQ_DMR,EQ_RAIL};
    const int types[]={0,PICKUP_GLOB,PICKUP_GLOCK,PICKUP_DMR,PICKUP_RAILGUN};
    for(int n=0;n<5;n++){
        tower_release(5);eq_reset(5);fpCarrierType[5]=types[n];equipment[5].gun[0]=guns[n];players[5].matter=100;
        assert(tower_eligible(5)&&eq_gun(5)==guns[n]);
        players[5].pos=tower_waypoint(t,9);players[5].vel=(Vector3){0};players[5].onGround=true;
        players[0].pos=tower_point(t,(Vector3){1.5f,1,-25});players[0].matter=100;players[0].invuln_timer=1000;
        t->occupant=5;towerRifles[5]=(TowerRifle){.tag=t->tag,.state=2};ai_select(5,0,players[0].pos,AI_ALERT_SEEN);players[5].last_shot_time=-100;
        uint32_t sequence=weaponVisuals[5].shot_sequence;
        for(int f=0;f<300 && weaponVisuals[5].shot_sequence==sequence;f++)step(1.f/60);
        fprintf(stderr,"Tower gun %s shot=%d\n",eqGuns[guns[n]].name,weaponVisuals[5].shot_sequence!=sequence);
        assert(weaponVisuals[5].shot_sequence!=sequence);
    }
    equipment[5].gun[0]=EQ_SHOTGUN;assert(!tower_eligible(5));
    puts("Garrison long-range fire, protected-player hold, exposed-player jump, return and weapon eligibility passed");
}
static void safety_checks(void){
    EnemyStructure *t=ready_tower();int n=voxel_count;
    Vector3 room=tower_point(t,(Vector3){1.5f,1,7.5f});
    assert(is_view_occluded_by_voxels(room,tower_point(t,(Vector3){-2,1,7.5f})));
    assert(is_view_occluded_by_voxels(room,tower_point(t,(Vector3){1.5f,1,5})));
    assert(!is_view_occluded_by_voxels(room,tower_point(t,(Vector3){1.5f,1,10})));
    players[5].pos=tower_point(t,(Vector3){-3,1,-3});
    tower_move(5,tower_point(t,(Vector3){-3,1,2}),1.f/60,false);
    assert(hypotf(players[5].vel.x,players[5].vel.z)>MOVE_SPEED*.99f);
    int v=enemy_voxel(t->ids[0]);Vector3 at=voxels[v].pos;remove_voxel_index(v);
    players[0].pos=v_add(at,(Vector3){0,.1f,0});assert(!tower_place(t));assert(voxel_count==n-1);
    players[0].pos=(Vector3){-20,1,-20};int obstruction=add_static_voxel_at_grid((int)floorf(at.x*2),(int)floorf(at.y*2),(int)floorf(at.z*2),RED,0);
    uint64_t identity=voxels[obstruction].identity;assert(!tower_place(t));assert(enemy_voxel(identity)>=0);
    remove_voxel_index(enemy_voxel(identity));assert(tower_place(t));assert(voxel_count==n);
    players[4].pos=(Vector3){-20,1,20};remove_voxel_index(enemy_voxel(t->ids[0]));for(int f=0;f<10;f++)step(1.f/60);assert(voxel_count==n-1);
    enemy_lab_spawn(6,ENEMY_TYPE_STANDARD,tower_waypoint(t,0));players[5].pos=tower_waypoint(t,0);towerRifles[5].retry=0;tower_rifle(5,.01f);tower_rifle(6,.01f);assert(t->occupant==5 && !towerRifles[6].tag);
    players[5].respawn_timer=9999;tower_tick(.01f);towerRifles[6].retry=0;tower_rifle(6,.01f);assert(t->occupant==6);
    puts("Repair actor/unrelated-terrain safety, range limit, single occupancy and replacement passed");
}
#ifndef TOWER_PREVIEW
int main(int argc,char **argv){
    assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,540,"Builder tower tests");init_pbd_thread_pool();assert(initialize_physics_backend());
    build_and_climb();garrison_checks();safety_checks();enemy_cleanup(true);assert(!voxel_count);clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Builder tower tests passed");return 0;
}

#endif
