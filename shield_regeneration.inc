/* A player-owned passive upgrade. Simulation is authoritative, never per viewport. */
static void eq_shield_damage(int i){
 if(i<0 || i>=MAX_PLAYERS)return;
 equipment[i].shieldDelay=eqShieldDelay;equipment[i].shieldRegenerating=false;
}
static void eq_shield_tick(float dt){
 if(netTransport.role==NET_ROLE_CLIENT || !isfinite(dt) || dt<=0)return;
 for(int i=0;i<MAX_PLAYERS;i++){
  Equipment *e=&equipment[i];Player *p=&players[i];e->shieldRegenerating=false;
  if(!(e->upgrades&EQ_SHIELD) || is_player_bot(i) || i>=activePlayers || p->respawn_timer>0 || gameState!=GAME_STATE_PLAYING ||
     (enemyLab && enemyLabPaused && !enemyLabStep) || (gameMode==GAME_MODE_FIREFIGHT &&
     (combatantLife[i]!=LIFE_ALIVE || firefightWaveState==FIREFIGHT_STATE_RESTART || firefightWaveState==FIREFIGHT_STATE_DEFEAT || firefightWaveState==FIREFIGHT_STATE_VICTORY)))continue;
  // Only the fraction of this tick after both delays expire may restore Matter.
  float blocked=fmaxf(e->shieldDelay,exSuppression[i]);
  e->shieldDelay=fmaxf(0,e->shieldDelay-dt);
  float recovery=fmaxf(0,dt-blocked),room=player_max_matter(p)*eqShieldCeiling-p->matter;
  if(recovery<=0 || room<=0 || ex_suppressed(i))continue;
  float before=p->matter;add_player_matter(p,fminf(room,eqShieldRate*recovery));
  e->shieldRegenerating=p->matter>before;
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
