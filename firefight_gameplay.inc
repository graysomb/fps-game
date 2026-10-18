/* Firefight orchestration; included after enemy initialization. */
static int ffRemaining[FF_TYPE_COUNT], ffSpawnOrdinal;
static float ffSpawnTimer, ffCleanupTimer;
static int ffEntrance;
static float ff_build_cost(int i) {return ff_bonus(i,UPGRADE_BUILD)?fmaxf(1,10.f/(1+teamUpgradeLevels[0])):MATTER_BUILD_COST;}
static float ff_throw_multiplier(int i) {return ff_bonus(i,UPGRADE_THROW)?1+.25f*teamUpgradeLevels[2]:1;}
static bool ff_upgrades_maxed(void) {return teamUpgradeLevels[0]>=9 && teamUpgradeLevels[1]>=8 && teamUpgradeLevels[2]>=8;}
static int ff_upgrade_choice(int choice,int direction) {
    for(int n=0;n<3;n++){choice=(choice+direction+3)%3;if(teamUpgradeLevels[choice]<teamUpgradeCaps[choice])return choice;}
    return -1;
}
static const char *ff_upgrade_label(int n) {
    unsigned level=teamUpgradeLevels[n];
    if(level>=teamUpgradeCaps[n])return TextFormat("%s Lv%u MAX",n==0?"Build":n==1?"Harvest":"Throw",level);
    if(n==0)return TextFormat("Build Lv%u: %.2g -> %.2g Matter",level,fmaxf(1,10.f/(1+level)),fmaxf(1,10.f/(2+level)));
    if(n==1)return TextFormat("Harvest Lv%u: %u -> %u Matter",level,10+5*level,15+5*level);
    return TextFormat("Throw Lv%u: %.2gx -> %.2gx",level,1+.25f*level,1.25f+.25f*level);
}

static EnemyType pick_firefight_enemy_type(int wave, int spawn_index) {
    (void)wave; (void)spawn_index;
    for(int t=0;t<FF_TYPE_COUNT;t++) if(ffRemaining[t]>0) return (EnemyType)t;
    return ENEMY_TYPE_STANDARD;
}

static void ff_restore_human(int i, float matter) {
    Player *p=&players[i];
    if(combatantLife[i]!=LIFE_DOWNED) grenadeInventory[i]=1;
    p->pos=pick_player_spawn(i); p->death_pos=p->pos;
    p->vel=(Vector3){0}; p->matter=matter; p->matterMax=MATTER_MAX_DEFAULT;
    p->isExposed=false; p->respawn_timer=0; p->invuln_timer=2;
    p->last_damage_time=(float)GetTime();
    p->tetherHolding=p->goldTetherHolding=false; p->tetherVoxel=-1;
    p->meleeSwingActive=p->meleeSwingHitApplied=false;
    p->meleeKnockbackActive=false; p->onGround=false;
    combatantLife[i]=LIFE_ALIVE;
    bleedTimer[i]=reviveProgress[i]=0; reviveOwner[i]=-1;
}

static void ff_resolve_death(int i) {
    if(enemyLab) {combatantLife[i]=LIFE_WAITING;players[i].respawn_timer=3;return;}
    if(firefightPractice) { ff_restore_human(i,MATTER_MAX_DEFAULT); return; }
    if(firefightTeamLives>0) {
        firefightTeamLives--; combatantLife[i]=LIFE_WAITING;
        players[i].respawn_timer=PLAYER_RESPAWN_TIME;
    } else {
        combatantLife[i]=LIFE_SPECTATING; players[i].respawn_timer=9999;
    }
    reviveOwner[i]=-1; reviveProgress[i]=0;
}

static bool ff_revive_held(int i) {
#ifdef FIREFIGHT_TESTING
    extern bool ffTestReviveHeld[4];
    return ffTestReviveHeld[i];
#else
    return playerInput[i]==INPUT_TYPE_GAMEPAD ?
        IsGamepadButtonDown(seat_device(i), GAMEPAD_BUTTON_RIGHT_THUMB) :
        IsKeyDown(i==0?KEY_X:KEY_COMMA);
#endif
}

static void ff_update_coop(float dt) {
    bool alive=false, pending=false;
    for(int i=0;i<firefightHumans;i++) {
        if(combatantLife[i]==LIFE_WAITING) {
            players[i].respawn_timer-=dt;
            if(players[i].respawn_timer<=0) ff_restore_human(i,MATTER_MAX_DEFAULT);
        }
        if(combatantLife[i]==LIFE_DOWNED) {
            float look=0;
            if(playerInput[i]==INPUT_TYPE_GAMEPAD) look=-GetGamepadAxisMovement(seat_device(i),GAMEPAD_AXIS_RIGHT_X)*localSeats[i].sensitivity;
            else look=(IsKeyDown(i==0?KEY_F:KEY_LEFT)?1.0f:0)-(IsKeyDown(i==0?KEY_H:KEY_RIGHT)?1.0f:0);
            players[i].death_yaw+=look*90*dt;
            bleedTimer[i]-=dt;
            int rescuer=-1;
            for(int j=0;j<firefightHumans;j++) {
                if(j==i || combatantLife[j]!=LIFE_ALIVE || players[j].matter<=FF_REVIVE_COST || !ff_revive_held(j)) continue;
                Vector3 delta=world_delta(players[i].death_pos,players[j].pos);
                float distance=v_length(delta);
                if(distance>2.0f || (float)GetTime()-players[j].last_damage_time<dt*1.5f) continue;
                if(distance>0.05f && first_voxel_hit((Ray){players[j].pos,v_mul(delta,1.0f/distance)},distance,-1)>=0) continue;
                /* One rescuer cannot progress two bodies. */
                bool busy=false;
                for(int k=0;k<i;k++) if(combatantLife[k]==LIFE_DOWNED && reviveOwner[k]==j) busy=true;
                if(!busy) { rescuer=j; break; }
            }
            if(rescuer<0 || reviveOwner[i]!=rescuer) reviveProgress[i]=0;
            reviveOwner[i]=rescuer;
            if(rescuer>=0) {
                reviveProgress[i]+=dt;
                if(reviveProgress[i]>=FF_REVIVE_SECONDS) {
                    Vector3 position=players[i].death_pos;
                    players[rescuer].matter-=FF_REVIVE_COST;
                    reviveCount[rescuer]++;
                    ff_restore_human(i,FF_REVIVE_MATTER);
                    players[i].pos=position;
                    play_sfx(SFX_SHIELD);
                }
            }
            if(combatantLife[i]==LIFE_DOWNED && bleedTimer[i]<=0) {
                players[i].deaths++; ff_resolve_death(i);
            }
        }
        alive |= combatantLife[i]==LIFE_ALIVE;
        pending |= combatantLife[i]==LIFE_WAITING;
    }
    if(!alive && !pending) {
        for(int i=0;i<firefightHumans;i++) if(combatantLife[i]==LIFE_DOWNED) {
            players[i].deaths++; ff_resolve_death(i);
            pending |= combatantLife[i]==LIFE_WAITING;
        }
        if(!pending) { firefightWaveState=FIREFIGHT_STATE_DEFEAT; gameState=GAME_STATE_GAMEOVER; }
    }
}

static void start_firefight_wave(int wave_number) {
    enemy_cleanup(true);
    firefightWave=wave_number<1?1:wave_number;
    firefightWaveState=FIREFIGHT_STATE_ACTIVE;
    firefightWaveBannerTimer=3.5f;
    firefightWaveEnemiesTotal=0; firefightEnemiesSpawned=0;
    for(int t=0;t<FF_TYPE_COUNT;t++) {
        ffRemaining[t]=firefightPractice ? (t==0?1:0) : ff_count(firefightWave,firefightHumans,t);
        firefightWaveEnemiesTotal+=ffRemaining[t];
    }
    firefightWaveEnemiesRemaining=firefightWaveEnemiesTotal;
    ffSpawnOrdinal=0; ffSpawnTimer=2; ffCleanupTimer=0; ffEntrance=0;
    for(int i=MAX_PLAYERS;i<MAX_COMBATANTS;i++) {
        if(goliathStates[i].count>0) explode_goliath_armor(i);
        players[i].respawn_timer=9999; combatantLife[i]=LIFE_UNUSED;
    }
}

static void reset_firefight_match(void) {
    firefightScore=firefightKillsTotal=0; firefightTeamLives=3;
    firefightIntermissionTimer=0; teamUpgrade=UPGRADE_NONE; memset(teamUpgradeLevels,0,sizeof(teamUpgradeLevels));
    memset(reviveCount,0,sizeof(reviveCount));
    for(int i=0;i<MAX_COMBATANTS;i++) goliathStates[i].count=0;
    for(int i=0;i<firefightHumans;i++) ff_restore_human(i,MATTER_MAX_DEFAULT);
    practiceStep=0;
    if(firefightPractice) {
        players[0].pos=(Vector3){0,BASE_EYE_HEIGHT+2,1.5f};
        players[0].yaw=0;players[0].pitch=0;
    }
    start_firefight_wave(1);
}

static void ff_wave_complete(void) {
    if(enemyLab) return;
    enemyCleanupPending=true;
    if(firefightPractice) {
        if(practiceStep>=4) { firefightWaveState=FIREFIGHT_STATE_VICTORY; gameState=GAME_STATE_GAMEOVER; }
        else { ffRemaining[0]=1; firefightWaveEnemiesTotal++; firefightWaveEnemiesRemaining=1; ffSpawnTimer=2; }
        return;
    }
    for(int i=0;i<firefightHumans;i++)grenadeInventory[i]=grenadeInventory[i]<1?1:grenadeInventory[i];
    if(firefightFoundry && !firefightPractice) {
        enemy_stock_scrap(false);
        for(int p=0;p<MAX_PICKUPS;p++)if(pickups[p].type==PICKUP_GOLD_TETHER || pickups[p].type==PICKUP_DYNAMIC_SHOT || pickups[p].type==PICKUP_GRENADE) {
            pickups[p].active=true;pickups[p].respawnTimer=0;
        }
    }
    if(!firefightEndless && firefightWave>=FF_WAVE_COUNT) {
        firefightWaveState=FIREFIGHT_STATE_VICTORY; gameState=GAME_STATE_GAMEOVER; return;
    }
    firefightWaveState=FIREFIGHT_STATE_INTERMISSION;
    firefightWaveBannerTimer=0; firefightIntermissionTimer=FF_PREP_SECONDS;
    if(firefightWave%2==0) firefightTeamLives++;
    for(int i=0;i<firefightHumans;i++) {
        ff_restore_human(i,MATTER_MAX_DEFAULT); prepVotes[i]=ff_upgrade_choice(2,1); prepReady[i]=false; prepMenuOpen[i]=true;
    }
}

static bool ff_spawn_position(int slot, Vector3 *out) {
    if(firefightPractice) { *out=(Vector3){0,BASE_EYE_HEIGHT+2,-1}; return true; }
    if(!firefightFoundry && !firefightPractice) { *out=pick_player_spawn(slot); return spawn_position_clear_ex(*out,spawn_clear_half(slot),4,slot); }
    for(int entrance=0;entrance<2;entrance++) {
        int side=(ffEntrance+entrance)%2 ? 1 : -1;
        for(int n=0;n<5;n++) {
            Vector3 pos={(float)side*36.0f,BASE_EYE_HEIGHT, (float)(n-2)*3.0f};
            pos=nudge_spawn_clear(pos,spawn_clear_half(slot),4.0f,slot);
            if(spawn_position_clear_ex(pos,spawn_clear_half(slot),4.0f,slot)) { *out=pos; return true; }
        }
    }
    return false;
}

static void update_firefight_logic(float dt) {
    if(enemyLab) return;
    enemy_tick(dt);
    if(firefightWaveState==FIREFIGHT_STATE_VICTORY || firefightWaveState==FIREFIGHT_STATE_DEFEAT) return;
    if(firefightWaveBannerTimer>0) firefightWaveBannerTimer-=dt;
    ff_update_coop(dt);
    if(gameState==GAME_STATE_GAMEOVER) return;
    if(firefightWaveState==FIREFIGHT_STATE_INTERMISSION) {
        if(ff_upgrades_maxed()){start_firefight_wave(firefightWave+1);return;}
        firefightIntermissionTimer-=dt;
        bool ready=true;
        for(int i=0;i<firefightHumans;i++) {
            bool pad=playerInput[i]==INPUT_TYPE_GAMEPAD; int d=seat_device(i);
            if(pad?IsGamepadButtonPressed(d,GAMEPAD_BUTTON_RIGHT_FACE_UP):IsKeyPressed(i==0?KEY_V:KEY_N)) prepMenuOpen[i]=!prepMenuOpen[i];
            if(!prepMenuOpen[i]) {ready &= prepReady[i];continue;}
            bool left=pad?IsGamepadButtonPressed(d,GAMEPAD_BUTTON_LEFT_FACE_LEFT):IsKeyPressed(i==0?KEY_F:KEY_LEFT);
            bool right=pad?IsGamepadButtonPressed(d,GAMEPAD_BUTTON_LEFT_FACE_RIGHT):IsKeyPressed(i==0?KEY_H:KEY_RIGHT);
            if(left || right) { prepVotes[i]=ff_upgrade_choice(prepVotes[i],left?-1:1); prepReady[i]=false; }
            if(pad?IsGamepadButtonPressed(d,GAMEPAD_BUTTON_RIGHT_FACE_DOWN):IsKeyPressed(i==0?KEY_SPACE:KEY_RIGHT_SHIFT)) prepReady[i]=!prepReady[i];
            ready &= prepReady[i];
        }
        if(firefightIntermissionTimer<=0 || (ready && firefightIntermissionTimer<=15)) {
            teamUpgrade=ff_vote(prepVotes,firefightHumans);
            if(teamUpgradeLevels[teamUpgrade]>=teamUpgradeCaps[teamUpgrade])teamUpgrade=ff_upgrade_choice(teamUpgrade,1);
            if((int)teamUpgrade>=0 && teamUpgrade<UPGRADE_NONE)teamUpgradeLevels[teamUpgrade]++;
            start_firefight_wave(firefightWave+1);
        }
        return;
    }
    if(firefightPractice && practiceStep==0) return;
    ffCleanupTimer+=dt;
    ffSpawnTimer-=dt;
    if(ffSpawnTimer>0) return;
    int living=0,goliaths=0,specialists=0,slot=-1;
    for(int i=MAX_PLAYERS;i<MAX_COMBATANTS;i++) {
        if(players[i].respawn_timer<=0) { living++; if(enemy_is_special(i)) specialists++; if(players[i].enemyType==ENEMY_TYPE_GOLIATH) goliaths++; }
        else if(slot<0) slot=i;
    }
    int waveCap=ff_wave_cap(firefightWave,firefightHumans);
    if(slot<0 || living>=waveCap) return;
    int type=-1;
    for(int n=0;n<FF_TYPE_COUNT;n++) { int t=(ffSpawnOrdinal+n)%FF_TYPE_COUNT;
        if(t==ENEMY_TYPE_SWARMER && (firefightWave-1)%FF_WAVE_COUNT==FF_WAVE_COUNT-1) continue;
        int specialistCap=(firefightWave-1)%FF_WAVE_COUNT==10?2:1;
        if(ffRemaining[t]>0 && (t!=2 || goliaths<2) && (t<3 || specialists<specialistCap)) {type=t;break;} }
    if(type<0) return;
    int pack=type==ENEMY_TYPE_SWARMER?3:1;
    for(int n=0;n<pack && ffRemaining[type]>0 && living<waveCap;n++) {
        slot=-1;
        for(int i=MAX_PLAYERS;i<MAX_COMBATANTS;i++) if(players[i].respawn_timer>0) {slot=i;break;}
        if(slot<0) break;
        players[slot].enemyType=(EnemyType)type;
        Vector3 position;
        if(!ff_spawn_position(slot,&position)) break;
        init_firefight_bot_entity(slot,(EnemyType)type,firefightWave==1?INPUT_TYPE_BOT_EASY:firefightWave<=3?INPUT_TYPE_BOT_MEDIUM:INPUT_TYPE_BOT_HARD);
        /* Attached armor follows its owner in the normal gameplay update. */
        players[slot].pos=position;
        players[slot].death_pos=position;
        ffRemaining[type]--; firefightEnemiesSpawned++; ffSpawnOrdinal++; living++;
    }
    ffEntrance=(firefightWave>=4)?1-ffEntrance:0;
    ffSpawnTimer=firefightWave<=2?4:2;
}
