#define TOWER_PREVIEW
#include "builder_towers_test.c"
int main(int argc,char **argv){
    assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(640,360,"Foundry builder regression");init_pbd_thread_pool();assert(initialize_physics_backend());
    for(int side=-1;side<=1;side+=2)for(int lane=-2;lane<=2;lane++){
        enemyLab=false;gameMode=GAME_MODE_FIREFIGHT;firefightFoundry=true;firefightPractice=false;firefightHumans=1;activePlayers=MAX_COMBATANTS;ResetGame();
        for(int i=4;i<MAX_COMBATANTS;i++){players[i].respawn_timer=9999;combatantLife[i]=LIFE_UNUSED;}
        init_firefight_bot_entity(4,ENEMY_TYPE_BUILDER,INPUT_TYPE_BOT_EASY);
        players[4].pos=(Vector3){side*36,1,lane*3};players[0].invuln_timer=1000;
        int frame=0;for(;frame<3600 && !built();frame++)step(1.f/60);
        fprintf(stderr,"Foundry entrance %d/%d: built=%d time=%.1f pos=%.1f,%.1f reason=%s\n",side,lane,built()!=NULL,frame/60.f,players[4].pos.x,players[4].pos.z,enemySpecial[4].reason);
        assert(built());
        EnemyStructure *t=built();
        for(int f=0;f<1800;f++)step(1.f/60);
        Vector3 local=crosscut_rotate(world_delta(players[4].pos,t->origin),(4-t->quarter)%4);
        fprintf(stderr,"Shelter local %.2f %.2f %.2f stage %d\n",local.x,local.y,local.z,t->shelterStep);
        assert(local.x>.75f && local.x<2.25f && local.z>6.5f && local.z<8.5f && local.y<1.5f);
    }
    enemy_cleanup(true);clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();return 0;
}
