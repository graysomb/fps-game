#define FIREFIGHT_TESTING
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
bool ffTestReviveHeld[4];
static void fixture(void){
 enemy_cleanup(true);clear_all_chunks();clear_world_voxels();clear_pickups();init_static_hash();worldPeriodic=false;
 memset(players,0,sizeof(players));memset(enemySpecial,0,sizeof(enemySpecial));
 gameMode=GAME_MODE_FIREFIGHT;gameState=GAME_STATE_PLAYING;activePlayers=MAX_COMBATANTS;firefightHumans=2;
 firefightPractice=enemyLab=equipmentLab=creativeModeActive=false;firefightCheckpoints=false;firefightWaveState=FIREFIGHT_STATE_ACTIVE;
 netTransport.role=NET_ROLE_OFFLINE;worldFloorHalfSize=40;fp_new_run();eq_reset_all();clarity_reset();ex_reset();
 for(int i=0;i<MAX_COMBATANTS;i++){playerInput[i]=i<4?INPUT_TYPE_KEYBOARD:INPUT_TYPE_BOT_EASY;combatantLife[i]=i<2?LIFE_ALIVE:LIFE_UNUSED;combatantTeam[i]=i<4?0:1;players[i].respawn_timer=i<2?0:9999;players[i].matter=players[i].matterMax=100;players[i].pos=(Vector3){i*3,1,0};reset_player_vitals(&players[i],false);ai_reset(i);}
}
static void recovery(void) {
    const int frameRates[] = {30, 60, 144};
    for (unsigned n = 0; n < sizeof(frameRates) / sizeof(frameRates[0]); ++n) {
        int fps = frameRates[n];
        fixture();
        apply_player_damage(0, 4, 85, false, true, false);
        assert(players[0].shields == 0 && players[0].health == 15);
        for (int i = 0; i < fps * 4; ++i) player_vitals_tick(1.f / fps);
        assert(players[0].shields == 0 && players[0].health == 15);
        for (int i = 0; i < fps * 2; ++i) player_vitals_tick(1.f / fps);
        assert(fabsf(players[0].shields - 25) < .01f && players[0].health == 15);
        for (int i = 0; i < fps * 3; ++i) player_vitals_tick(1.f / fps);
        assert(players[0].shields == 75 && fabsf(players[0].health - 20) < .01f);
        assert(players[0].matter == 100);
    }

    fixture();
    assert(eq_pickup(0, PICKUP_SHIELD_REGEN, false));
    apply_player_damage(0, 4, 85, false, true, false);
    player_vitals_tick(6);
    assert(players[0].shields == 35 && players[0].health == 15);
    apply_matter_damage(0, 4, 10);
    player_vitals_tick(4.9f);
    assert(players[0].shields == 25 && players[0].health == 15);
    player_vitals_tick(.2f);
    assert(fabsf(players[0].shields - 28.5f) < .01f);

    players[0].recoveryDelay = 0;
    exSuppression[0] = 5;
    player_vitals_tick(4);
    assert(fabsf(players[0].shields - 28.5f) < .01f && players[0].health == 15);
    exSuppression[0] = 0;
    player_vitals_tick(3);
    assert(players[0].shields == 75 && players[0].health > 15);

    players[0].recoveryDelay = 0;
    players[0].invuln_timer = 2;
    apply_matter_damage(0, 4, 10);
    players[0].invuln_timer = 0;
    apply_matter_damage(0, 1, 10);
    assert(players[0].shields == 75 && players[0].recoveryDelay == 0);
    players[4].respawn_timer = 0;
    apply_matter_damage(4, 0, 10);
    assert(players[4].matter == 90);
    netTransport.role = NET_ROLE_CLIENT;
    apply_player_damage(0, 4, 100, false, true, false);
    player_vitals_tick(100);
    assert(players[0].shields == 75 && players[0].health < 25);
    netTransport.role = NET_ROLE_OFFLINE;
    assert(!eq_pickup(4, PICKUP_SHIELD_REGEN, false));
    puts("Recovery delays, rate, frame rates, suppression and host authority passed");
}

static void combat(void) {
    fixture();
    gameMode = GAME_MODE_DEATHMATCH;
    for (int n = 0; n < 3; ++n) {
        apply_matter_damage(0, 1, eqGuns[EQ_DEFAULT].damage);
    }
    eq_damage(1, 0, eqGuns[EQ_GLOB].damage, false, DEATH_GLOB, players[1].pos);
    assert(players[0].shields == 0 && players[0].health == 25 && players[0].matter == 100);
    assert(players[0].shieldBreakFlash > 0 && players[0].healthHitFlash == 0);
    clarity_actor_damage(0, 1, 20, true);
    assert(players[0].health == 5 && players[0].healthHitFlash > 0);
    clarity_actor_damage(0, 1, 20, true);
    assert(players[0].health == 0 && players[0].respawn_timer > 0);

    fixture();
    gameMode = GAME_MODE_DEATHMATCH;
    for (int n = 0; n < 4; ++n) {
        eq_damage(1, 0, eqGuns[EQ_DMR].damage, true, DEATH_DMR, players[1].pos);
        assert((players[0].respawn_timer > 0) == (n == 3));
    }

    fixture();
    apply_player_damage(0, 4, 110, false, true, false);
    assert(players[0].health == 0 && combatantLife[0] == LIFE_DOWNED);
    player_vitals_tick(20);
    assert(players[0].health == 0 && players[0].shields == 0);
    players[1].pos = players[0].death_pos;
    players[1].last_damage_time = -1000;
    ffTestReviveHeld[1] = true;
    ff_update_coop(FF_REVIVE_SECONDS);
    ffTestReviveHeld[1] = false;
    assert(combatantLife[0] == LIFE_ALIVE && players[0].health == 12.5f);
    assert(players[0].shields == 0 && players[0].matter == FF_REVIVE_MATTER);

    fixture();
    players[0].health = 10;
    players[0].shields = 12;
    eq_spend(0, 100);
    assert(players[0].matter == 0 && !players[0].isExposed);
    add_player_matter(&players[0], 15);
    assert(players[0].matter == 15 && players[0].shields == 12 && players[0].health == 10);
    pickups[0] = (Pickup){ .pos = players[0].pos, .active = true, .type = PICKUP_HEALTH };
    pickups[1] = (Pickup){ .pos = players[0].pos, .active = true, .type = PICKUP_AMMO };
    update_pickups(0);
    assert(players[0].matter == 40 && players[0].shields == 12 && players[0].health == 25);
    puts("Shield stripping, debris, four-hit DMR, death, revive and independent resources passed");
}
static void drops(void){
 fixture();
 for(int i=4;i<7;i++){
  init_firefight_bot_entity(i,ENEMY_TYPE_SHIELDBEARER,INPUT_TYPE_BOT_MEDIUM);fp_spawn(i,ENEMY_TYPE_SHIELDBEARER);
  assert(!fp_carrier(i) && !(equipment[i].upgrades&EQ_SHIELD));
  players[i].invuln_timer=0;players[i].isExposed=true;kill_player(i,0,false,false);
 }
 assert(fpCount==3);for(int n=0;n<3;n++)assert(fpItems[n].type==PICKUP_SHIELD_REGEN && fpItems[n].owner==-1);
 uint32_t id=fpItems[0].id;assert(fp_collect(0,0,false));assert(!fp_collect(0,1,false));assert(fpItems[1].owner==-1);
 assert(fp_collect(1,1,false));assert(fpBannerWrite==1);assert(eq_gun(0)==EQ_DEFAULT && equipment[0].gun[1]==EQ_EMPTY);
 players[0].invuln_timer=0;players[0].matter=0;players[0].isExposed=true;kill_player(0,4,false,false);
 assert(combatantLife[0]==LIFE_DOWNED && fpItems[0].owner==0 && eq_upgrade(0,EQ_SHIELD));
 ff_restore_human(0,35);assert(eq_upgrade(0,EQ_SHIELD));
 players[0].death_pos=players[0].pos;ff_resolve_death(0);assert(!eq_upgrade(0,EQ_SHIELD) && fpItems[0].owner==-1 && fpItems[0].id==id);
 ff_restore_human(0,100);assert(fp_collect(0,0,false));
 firefightCheckpoints=true;firefightWave=4;fp_capture();assert(fpCheckpoint.equipment[0].upgrades&EQ_SHIELD);
 eq_reset(0);fp_drop(&fpItems[0],players[0].pos);fp_restore();assert(eq_upgrade(0,EQ_SHIELD) && fpItems[0].owner==0 && fpItems[0].id==id && players[0].shields==75 && players[0].health==25);
 fp_new_run();eq_reset_all();assert(fpCount==0 && !eq_upgrade(0,EQ_SHIELD));
 puts("Drops: every bearer, native AI, duplicates, one celebration, downs/revives/death, checkpoint and reset passed");
}
static void networking(void) {
    fixture();
    uint8_t data[FPS_NET_MAX_PACKET];
    NetWriter writer;
    net_writer_init(&writer, data, sizeof(data));
    net_write_header(&writer, NET_MSG_PLAYER_STATE, 0, 17, 1);
    net_write_u8(&writer, 1);
    players[0].health = 17;
    players[0].shields = 31;
    players[0].matter = 42;
    players[0].matterMax = 130;
    players[0].rechargeFlags = RECHARGE_SHIELDS;
    net_write_player_state(&writer, 0);
    assert(!writer.failed);

    reset_player_vitals(&players[0], false);
    netTransport.role = NET_ROLE_CLIENT;
    netTransport.session_id = 17;
    netServerTick = 0;
    net_on_receive(&netTransport, 0, FPS_NET_CHANNEL_SNAPSHOT, data, writer.length, NULL);
    assert(players[0].health == 17 && players[0].shields == 31);
    assert(players[0].matter == 42 && players[0].matterMax == 130);
    assert(players[0].rechargeFlags == RECHARGE_SHIELDS);
    assert(players[0].shieldHitFlash > 0 && players[0].healthHitFlash > 0);

    net_writer_init(&writer, data, sizeof(data));
    players[0].shields = NAN;
    net_write_player_state(&writer, 0);
    NetReader reader;
    NetPlayerWireState state;
    int slot;
    net_reader_init(&reader, data, writer.length);
    assert(!net_read_player_state(&reader, &slot, &state));
    netTransport.role = NET_ROLE_OFFLINE;
    puts("Player snapshot round trip and malformed-vitals rejection passed");
}
static void captures(void){
 fixture();weapon_renderer_init();character_renderer_init();eq_pickup(0,PICKUP_SHIELD_REGEN,false);players[0].matter=20;players[0].shields=42;player_vitals_tick(.1f);
 drawViewPlayerIndex=0;weaponFrameSampled=true;weaponFrameTime=1;
 for(int small=0;small<2;small++)for(int aim=0;aim<2;aim++){
  int w=small?480:960,h=small?270:540;RenderTexture2D rt=LoadRenderTexture(w,h);
  Camera3D camera={{3,2.3f,6},{0,1,0},{0,1,0},aim?30:60,CAMERA_PERSPECTIVE};
  BeginTextureMode(rt);ClearBackground((Color){25,32,43,255});BeginMode3D(camera);
  DrawPlane((Vector3){0,0,0},(Vector2){20,20},(Color){60,67,75,255});
  CharacterPose pose=character_pose((Vector3){-.7f,BASE_EYE_HEIGHT,0},145,0,0,true,false,BLUE,1);pose.effect=SKYBLUE;pose.effectAmount=.2f;character_draw(pose,false);
  eq_shield_pickup((Vector3){.8f,1.2f,0},1,WEAPON_COLOR);eq_shield_pickup((Vector3){.8f,1.2f,0},1,WEAPON_GLASS);
  EndMode3D();draw_hud_bars(0,&players[0],w,h);eq_hud(0,w,h);DrawCircle(w/2,h/2,2,WHITE);EndTextureMode();
  Image im=LoadImageFromTexture(rt.texture);ImageFlipVertical(&im);ExportImage(im,TextFormat(".build/shield-regen-%d-%d.png",small,aim));UnloadImage(im);UnloadRenderTexture(rt);
 }
 RenderTexture2D split=LoadRenderTexture(960,540),tile=LoadRenderTexture(480,270);
 for(int i=0;i<4;i++){
  combatantLife[i]=LIFE_ALIVE;players[i].respawn_timer=0;players[i].matter=20+i*8;eq_pickup(i,PICKUP_SHIELD_REGEN,false);
  players[i].shields=i==1?0:23+i*5;players[i].health=10+i*4;players[i].rechargeFlags=i==1?0:RECHARGE_SHIELDS;
  Camera3D c={{0,1.8f,4},{0,1,0},{0,1,0},60,CAMERA_PERSPECTIVE};
  drawViewPlayerIndex=i;
  BeginTextureMode(tile);ClearBackground((Color){25,32,43,255});BeginMode3D(c);
  DrawPlane((Vector3){0,0,0},(Vector2){20,20},(Color){60,67,75,255});
  CharacterPose pose=character_pose((Vector3){-.7f,BASE_EYE_HEIGHT,0},135+i*20,0,0,true,false,player_palette_color(i),1);
  if(i!=1){pose.effect=SKYBLUE;pose.effectAmount=.2f;}character_draw(pose,false);
  eq_shield_pickup((Vector3){.8f,1.2f,0},1,WEAPON_COLOR);eq_shield_pickup((Vector3){.8f,1.2f,0},1,WEAPON_GLASS);
  if(i==3)DrawCube((Vector3){.8f,1.2f,.7f},1,1.2f,.3f,GRAY);
  EndMode3D();draw_player_radar(i,480,270);draw_hud_bars(i,&players[i],480,270);
  combat_ui_abilities(i,480,270);eq_hud(i,480,270);EndTextureMode();
  BeginTextureMode(split);if(i==0)ClearBackground(BLACK);DrawTextureRec(tile.texture,(Rectangle){0,0,480,-270},(Vector2){(i%2)*480,(i/2)*270},WHITE);EndTextureMode();
 }
 Image four=LoadImageFromTexture(split.texture);ImageFlipVertical(&four);ExportImage(four,".build/shield-regen-four-player.png");UnloadImage(four);UnloadRenderTexture(tile);UnloadRenderTexture(split);
 character_renderer_shutdown();weapon_renderer_shutdown();weaponFrameSampled=false;
}
int main(int argc,char**argv){
 assert(parse_physics_arguments(argc,argv));SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(960,540,"Shield regeneration tests");
 init_pbd_thread_pool();assert(initialize_physics_backend());recovery();combat();drops();networking();captures();
 enemy_cleanup(true);fp_shutdown();clear_world_voxels();shutdown_pbd_thread_pool();gpu_physics_shutdown();CloseWindow();puts("Shield regeneration passed");return 0;
}
