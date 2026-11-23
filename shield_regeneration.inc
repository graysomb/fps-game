static void player_vitals_tick(float dt) {
    if (!isfinite(dt) || dt <= 0.0f) return;
    for (int i = 0; i < MAX_PLAYERS && i < activePlayers; ++i) {
        if (!uses_player_vitals(i)) continue;
        Player *p = &players[i];
        p->shieldHitFlash = fmaxf(0.0f, p->shieldHitFlash - dt);
        p->shieldBreakFlash = fmaxf(0.0f, p->shieldBreakFlash - dt);
        p->healthHitFlash = fmaxf(0.0f, p->healthHitFlash - dt);
        if (netTransport.role == NET_ROLE_CLIENT) continue;

        p->rechargeFlags = 0;
        if (gameState != GAME_STATE_PLAYING || p->respawn_timer > 0.0f || p->health <= 0.0f ||
            (enemyLab && enemyLabPaused && !enemyLabStep) || (gameMode == GAME_MODE_FIREFIGHT &&
            (combatantLife[i] != LIFE_ALIVE || firefightWaveState == FIREFIGHT_STATE_RESTART ||
             firefightWaveState == FIREFIGHT_STATE_DEFEAT || firefightWaveState == FIREFIGHT_STATE_VICTORY))) continue;

        float suppression = fmaxf(0.0f, exSuppression[i]);
        float shieldStart = fmaxf(suppression, p->recoveryDelay - (PLAYER_HEALTH_DELAY - PLAYER_SHIELD_DELAY));
        float healthStart = fmaxf(suppression, p->recoveryDelay);
        p->recoveryDelay = fmaxf(0.0f, p->recoveryDelay - dt);

        if (p->shields < PLAYER_SHIELDS_MAX) {
            float rate = (equipment[i].upgrades & EQ_SHIELD) ? PLAYER_SHIELD_UPGRADED_RATE : PLAYER_SHIELD_RATE;
            // Health can use only the part of this frame after shields reach full.
            healthStart = fmaxf(healthStart, shieldStart + (PLAYER_SHIELDS_MAX - p->shields) / rate);
            if (dt > shieldStart) {
                p->shields = fminf(PLAYER_SHIELDS_MAX, p->shields + rate * (dt - shieldStart));
                if (p->shields < PLAYER_SHIELDS_MAX) p->rechargeFlags |= RECHARGE_SHIELDS;
            }
        }
        if (p->health < PLAYER_HEALTH_MAX && dt > healthStart) {
            p->health = fminf(PLAYER_HEALTH_MAX, p->health + PLAYER_HEALTH_RATE * (dt - healthStart));
            if (p->health < PLAYER_HEALTH_MAX) p->rechargeFlags |= RECHARGE_HEALTH;
        }
    }
}
static void eq_shield_icon(int x,int y,Color color){
 Vector2 v[5]={{x+12,y},{x+24,y+5},{x+21,y+17},{x+12,y+26},{x+3,y+17}};
 Vector2 left={x,y+5};
 DrawLineEx(left,v[0],2,color);for(int n=0;n<4;n++)DrawLineEx(v[n],v[n+1],2,color);DrawLineEx(v[4],left,2,color);
 DrawLineEx((Vector2){x+12,y+6},(Vector2){x+12,y+19},2,color);
 DrawLineEx((Vector2){x+6,y+12},(Vector2){x+18,y+12},2,color);
}
static void eq_shield_pickup(Vector3 pos,float time,WeaponRenderPass pass){
 // Six fixed vertices: no mesh allocations, with identical geometry in bloom.
 Vector3 local[6]={{0,.46f,0},{.37f,.28f,0},{.3f,-.16f,0},{0,-.46f,0},{-.3f,-.16f,0},{-.37f,.28f,0}},v[6];
 Matrix m=MatrixMultiply(MatrixRotateY(time*.6f),MatrixTranslate(pos.x,pos.y,pos.z));
 for(int n=0;n<6;n++)v[n]=Vector3Transform(local[n],m);
 Color rim={90,225,255,255};
 for(int n=0;n<6;n++){
  int next=(n+1)%6;
  if(pass==WEAPON_GLASS){Color c=Fade(SKYBLUE,.24f);DrawTriangle3D(pos,v[n],v[next],c);DrawTriangle3D(pos,v[next],v[n],c);}
  else DrawCylinderEx(v[n],v[next],.025f,.025f,4,rim);
 }
 if(pass!=WEAPON_GLASS){
  Vector3 a=Vector3Transform((Vector3){0,-.18f,0},m),b=Vector3Transform((Vector3){0,.22f,0},m);
  DrawCylinderEx(a,b,.022f,.022f,4,rim);
  a=Vector3Transform((Vector3){-.16f,.02f,0},m);b=Vector3Transform((Vector3){.16f,.02f,0},m);DrawCylinderEx(a,b,.022f,.022f,4,rim);
 }
}
