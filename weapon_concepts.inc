/* Cached, flat-shaded concept meshes. Local forward is -Z, matching the
 * original weapons. No gameplay state or frame-time mesh allocation here. */
static void wc_box(WrBuilder *b,float x,float y,float z,float w,float h,float d,int material) {
    float bevel=fminf(.045f,fminf(w,fminf(h,d))*.22f);
    wr_prism(b,(Vector3){x,y,z},w,h,d,bevel,material,0);
}
static void wc_grip(WrBuilder *b) {
    wr_prism(b,(Vector3){0,-.65f,.28f},.23f,.68f,.28f,.035f,WR_DARK,-.25f);
    wc_box(b,0,-.98f,.36f,.29f,.12f,.34f,WR_SILVER);
}
static void wc_build(WrBuilder *b,int kind) {
    if(kind!=WEAPON_ADHESIVE_GRENADE && kind!=WEAPON_SWORD)wc_grip(b);
    switch(kind) {
    case WEAPON_IMPULSE_SHOTGUN:
        wc_box(b,0,0,-.15f,.85f,.47f,1.2f,WR_SILVER);
        wc_box(b,0,0,-.82f,1.12f,.48f,.25f,WR_SILVER);
        for(int n=-1;n<=1;n++){
            wc_box(b,n*.34f,0,-.96f,.29f,.32f,.12f,WR_DARK);
            wc_box(b,n*.34f,0,-1.03f,.20f,.22f,.04f,WR_ENERGY);
        }
        for(int n=-1;n<=1;n+=2)wc_box(b,n*.43f,0,-.1f,.04f,.18f,.65f,WR_ENERGY);
        wc_box(b,0,.33f,.16f,.25f,.18f,.52f,WR_DARK);
        break;
    case WEAPON_CUTTER:
        wc_box(b,0,0,-.03f,.42f,.40f,1.0f,WR_SILVER);
        wc_box(b,0,.24f,.02f,.28f,.12f,.7f,WR_DARK);
        for(int n=-1;n<=1;n+=2){
            wc_box(b,0,n*.19f,-.91f,.20f,.14f,.8f,WR_SILVER);
            wc_box(b,0,n*.19f,-1.32f,.16f,.10f,.035f,WR_ENERGY);
            wc_box(b,n*.215f,0,.1f,.03f,.12f,.4f,WR_ENERGY);
        }
        wc_box(b,0,0,-.9f,.035f,.035f,.8f,WR_ENERGY);
        break;
    case WEAPON_ADHESIVE_GRENADE:
        wc_box(b,0,0,0,.62f,.62f,.62f,WR_DARK);
        wc_box(b,0,0,0,.65f,.065f,.65f,WR_ENERGY);
        for(int x=-1;x<=1;x+=2)for(int y=-1;y<=1;y+=2){
            wc_box(b,x*.19f,y*.17f,-.34f,.23f,.20f,.08f,WR_CHARCOAL);
            /* Pads get fixed purple vertex color, independent of the fuse. */
            int first=b->count;
            wc_box(b,x*.19f,y*.17f,-.39f,.20f,.17f,.04f,WR_SILVER);
            for(int v=first;v<b->count;v++){b->colors[v*4]=132;b->colors[v*4+1]=65;b->colors[v*4+2]=205;}
        }
        for(int x=-1;x<=1;x+=2)for(int y=-1;y<=1;y+=2)
            wc_box(b,x*.28f,y*.28f,-.30f,.12f,.12f,.14f,WR_SILVER);
        break;
    case WEAPON_RAILGUN:
        wc_box(b,0,0,.05f,.46f,.40f,.85f,WR_SILVER);
        for(int side=-1;side<=1;side+=2){
            wc_box(b,0,side*.22f,-1.05f,.23f,.17f,1.65f,WR_SILVER);
            for(int n=0;n<5;n++)wc_box(b,0,side*.125f,-.45f-n*.28f,.17f,.035f,.18f,WR_ENERGY);
        }
        wc_box(b,0,.36f,.08f,.23f,.18f,.38f,WR_DARK);
        wc_box(b,0,.36f,-.12f,.16f,.1f,.035f,WR_ENERGY);
        wc_box(b,0,0,.65f,.28f,.48f,.14f,WR_DARK);
        break;
    case WEAPON_ROCKET_LAUNCHER:
        /* Four walls form a real open square bore, not a painted end cap. */
        for(int side=-1;side<=1;side+=2){
            wc_box(b,side*.33f,0,-.28f,.16f,.78f,1.85f,WR_SILVER);
            wc_box(b,0,side*.33f,-.28f,.52f,.16f,1.85f,WR_SILVER);
        }
        wc_box(b,0,0,-.6f,.44f,.44f,.75f,WR_DARK);
        wc_box(b,0,0,-1.0f,.32f,.32f,.16f,WR_ENERGY);
        wc_box(b,0,.52f,0,.22f,.20f,.65f,WR_DARK);
        wc_box(b,0,.65f,0,.15f,.05f,.25f,WR_ENERGY);
        wc_box(b,0,0,.69f,.76f,.76f,.13f,WR_DARK);
        break;
    case WEAPON_SWORD:
        wc_box(b,0,0,.52f,.19f,.18f,.65f,WR_DARK);
        wc_box(b,0,0,.88f,.26f,.24f,.12f,WR_SILVER);
        wc_box(b,0,0,.96f,.15f,.14f,.035f,WR_ENERGY);
        wc_box(b,0,0,.15f,.65f,.23f,.14f,WR_SILVER);
        wc_box(b,0,0,-.72f,.38f,.13f,1.62f,WR_SILVER);
        wc_box(b,.035f,.075f,-.66f,.16f,.025f,1.42f,WR_DARK);
        wc_box(b,-.20f,0,-.72f,.075f,.13f,1.62f,WR_ENERGY);
        /* Clip the broad blade diagonally to a point without extra detail. */
        for(int v=0;v<b->count;v++){
            float *p=b->vertices+v*3;
            if(p[2]<-1.15f)p[0]=-.20f+(p[0]+.20f)*fmaxf(.08f,1-( -p[2]-1.15f)/.4f);
        }
        for(int v=0;v<b->count;v+=3){
            Vector3 a,c,d;memcpy(&a,b->vertices+v*3,sizeof(a));memcpy(&c,b->vertices+(v+1)*3,sizeof(c));memcpy(&d,b->vertices+(v+2)*3,sizeof(d));
            Vector3 n=Vector3Normalize(Vector3CrossProduct(Vector3Subtract(c,a),Vector3Subtract(d,a)));
            for(int p=0;p<3;p++)memcpy(b->normals+(v+p)*3,&n,sizeof(n));
        }
        break;
    case WEAPON_REINFORCED_BUILDER:
        wc_box(b,0,0,.04f,.58f,.48f,.85f,WR_SILVER);
        wc_box(b,0,0,-.42f,.38f,.32f,.08f,WR_ENERGY);
        for(int side=-1;side<=1;side+=2)for(int top=-1;top<=1;top+=2){
            wc_box(b,side*.48f,top*.26f,-.45f,.72f,.14f,.16f,WR_SILVER);
            wc_box(b,side*.78f,top*.26f,-.62f,.14f,.14f,.35f,WR_DARK);
            wc_box(b,side*.78f,top*.26f,-.81f,.12f,.12f,.04f,WR_ENERGY);
        }
        for(int n=0;n<3;n++)wc_box(b,.30f,(n-1)*.13f,.1f,.025f,.085f,.13f,WR_ENERGY);
        break;
    default: break;
    }
}
static Vector3 wc_energy(int kind) {
    switch(kind){
    case WEAPON_IMPULSE_SHOTGUN:return (Vector3){1,.32f,.015f};
    case WEAPON_CUTTER:return (Vector3){.35f,1,.025f};
    case WEAPON_ADHESIVE_GRENADE:return (Vector3){1,.55f,.015f};
    case WEAPON_RAILGUN:return (Vector3){.65f,.08f,1};
    case WEAPON_ROCKET_LAUNCHER:return (Vector3){1,.12f,.015f};
    case WEAPON_REINFORCED_BUILDER:return (Vector3){1,.75f,.035f};
    default:return (Vector3){0,.95f,1};
    }
}
