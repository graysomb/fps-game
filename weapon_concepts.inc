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
/* Small inset details use six faces; broad silhouette pieces keep their bevels. */
static void wc_detail(WrBuilder *b,float x,float y,float z,float w,float h,float d,int m){
    Vector3 c={x,y,z};float a=w*.5f,e=h*.5f,f=d*.5f;
    Vector3 p[8]={{-a,-e,-f},{a,-e,-f},{a,e,-f},{-a,e,-f},{-a,-e,f},{a,-e,f},{a,e,f},{-a,e,f}};
    for(int i=0;i<8;i++)p[i]=Vector3Add(p[i],c);
    const int faces[6][4]={{0,3,2,1},{4,5,6,7},{0,4,7,3},{1,2,6,5},{0,1,5,4},{3,7,6,2}};
    for(int i=0;i<6;i++)wr_quad(b,p[faces[i][0]],p[faces[i][1]],p[faces[i][2]],p[faces[i][3]],m,1);
}
static void wc_handle(WrBuilder *b,float y,float z,float length){
    wc_detail(b,0,y,z,.18f,.10f,length,WR_DARK);
    for(int n=-1;n<=1;n+=2)wc_detail(b,0,y-.13f,z+n*(length*.5f-.05f),.18f,.26f,.10f,WR_DARK);
}
static void wc_side_frame(WrBuilder *b,float x,float y,float z,float w,float h){
    int start=b->count;wr_tube(b,(Vector3){0},w,h,.07f,.035f,.015f,WR_SILVER);
    for(int v=start;v<b->count;v++){
        float *p=b->vertices+v*3,*n=b->normals+v*3;float px=p[0],nx=n[0];
        p[0]=x+p[2];p[1]+=y;p[2]=z-px;n[0]=n[2];n[2]=-nx;
    }
}
static void wc_rocket(WrBuilder *b){
    wc_box(b,0,0,0,.26f,.26f,.65f,WR_CHARCOAL);
    int start=b->count;wc_box(b,0,0,-.43f,.26f,.26f,.25f,WR_SILVER);
    for(int v=start;v<b->count;v++){
        float *p=b->vertices+v*3;float t=Clamp((p[2]+.56f)/.25f,.04f,1);p[0]*=t;p[1]*=t;
    }
    /* Recompute the tapered nose's flat face normals. */
    for(int v=start;v<b->count;v+=3){Vector3 a,c,d;memcpy(&a,b->vertices+v*3,sizeof(a));memcpy(&c,b->vertices+(v+1)*3,sizeof(c));memcpy(&d,b->vertices+(v+2)*3,sizeof(d));Vector3 n=Vector3Normalize(Vector3CrossProduct(Vector3Subtract(c,a),Vector3Subtract(d,a)));for(int k=0;k<3;k++)memcpy(b->normals+(v+k)*3,&n,sizeof(n));}
    wc_detail(b,0,0,-.26f,.28f,.28f,.06f,WR_ENERGY);
    for(int side=-1;side<=1;side+=2){wc_detail(b,side*.19f,0,.22f,.18f,.035f,.28f,WR_SILVER);wc_detail(b,0,side*.19f,.22f,.035f,.18f,.28f,WR_SILVER);}
    wc_box(b,0,0,.35f,.19f,.19f,.04f,WR_ENERGY);
}
static void wc_flash(WrBuilder *b,bool rear){
    if(wr_fallback)return; // The static vertex-color fallback never draws a frozen flash.
    b->part=rear?12:11;
    float z=rear?.85f:-1.10f,tip=rear?1.45f:-1.60f;
    for(int n=-1;n<=1;n++){
        Vector3 a={n*.25f-.09f,0,z},c={n*.25f+.09f,0,z},d={n*.42f,0,tip};
        wr_triangle(b,a,c,d,WR_ENERGY,1,0);wr_triangle(b,c,a,d,WR_ENERGY,1,0);
    }b->part=0;
}
static void wc_build(WrBuilder *b,int kind) {
    if(kind!=WEAPON_ADHESIVE_GRENADE && kind!=WEAPON_SWORD)wc_grip(b);
    switch(kind) {
    case WEAPON_IMPULSE_SHOTGUN:
        wc_box(b,0,0,-.1f,.88f,.48f,1.10f,WR_CHARCOAL);
        for(int side=-1;side<=1;side+=2)wc_box(b,side*.34f,.20f,-.1f,.24f,.14f,1.12f,WR_SILVER);
        wc_box(b,0,.20f,.19f,.50f,.14f,.54f,WR_SILVER);
        wc_detail(b,0,.19f,-.37f,.40f,.025f,.40f,WR_ENERGY);
        wc_detail(b,0,.21f,-.35f,.12f,.025f,.22f,WR_DARK);
        wc_box(b,0,-.22f,-.1f,.92f,.12f,1.12f,WR_SILVER);
        for(int side=-1;side<=1;side+=2){
            wc_side_frame(b,side*.46f,0,-.12f,.70f,.30f);
            wc_detail(b,side*.465f,0,-.12f,.018f,.20f,.58f,WR_ENERGY);
            wc_detail(b,side*.478f,.025f,-.12f,.022f,.12f,.40f,WR_DARK);
            wc_box(b,side*.47f,0,.34f,.10f,.44f,.16f,WR_SILVER);
        }
        /* Rear cartridge stays readable in the actual first-person view. */
        wr_tube(b,(Vector3){0,0,.49f},.88f,.46f,.12f,.055f,.025f,WR_SILVER);
        wc_detail(b,0,0,.48f,.065f,.25f,.025f,WR_ENERGY);
        for(int n=-1;n<=1;n+=2)wc_detail(b,n*.24f,0,.48f,.12f,.045f,.025f,WR_DARK);
        wc_handle(b,.51f,.05f,.58f);
        b->part=1;
        wc_box(b,0,0,-.70f,1.12f,.53f,.20f,WR_SILVER);
        for(int n=-1;n<=1;n++){
            wr_tube(b,(Vector3){n*.35f,0,-.91f},.30f,.36f,.29f,.055f,.025f,WR_DARK);
            wc_detail(b,n*.35f,0,-1.01f,.17f,.22f,.025f,WR_ENERGY);
        }
        b->part=0;wc_flash(b,false);
        break;
    case WEAPON_CUTTER:
        wc_box(b,0,0,.02f,.48f,.42f,.90f,WR_CHARCOAL);
        for(int side=-1;side<=1;side+=2){
            wc_box(b,side*.23f,0,.05f,.09f,.38f,.75f,WR_SILVER);
            wc_detail(b,side*.28f,.02f,.08f,.015f,.12f,.44f,WR_DARK);
            for(int n=0;n<3;n++)wc_detail(b,side*.29f,.02f,.21f-n*.14f,.016f,.07f,.08f,WR_ENERGY);
        }
        wr_tube(b,(Vector3){0,0,.50f},.48f,.40f,.10f,.05f,.02f,WR_SILVER);
        for(int n=-1;n<=1;n+=2)wc_detail(b,n*.075f,0,.50f,.045f,.17f,.02f,WR_ENERGY);
        wc_box(b,0,.25f,.02f,.28f,.13f,.7f,WR_SILVER);
        for(int n=0;n<5;n++)wc_detail(b,0,.32f,.25f-n*.12f,.20f,.015f,.055f,WR_DARK);
        wc_box(b,0,0,-.49f,.38f,.52f,.16f,WR_DARK);
        wc_box(b,0,-.38f,-.31f,.20f,.26f,.22f,WR_CHARCOAL);
        for(int n=-1;n<=1;n+=2){
            b->part=n>0?2:3;
            wc_box(b,0,n*.23f,-.91f,.27f,.16f,.82f,WR_SILVER);
            wc_detail(b,0,n*.148f,-.96f,.15f,.018f,.61f,WR_DARK);
            for(int k=0;k<4;k++)wc_detail(b,0,n*.134f,-.72f-k*.16f,.10f,.014f,.09f,WR_ENERGY);
            wc_box(b,0,n*.23f,-1.31f,.29f,.19f,.10f,WR_DARK);
            wc_detail(b,0,n*.23f,-1.37f,.17f,.10f,.022f,WR_ENERGY);
        }b->part=0;
        break;
    case WEAPON_ADHESIVE_GRENADE:
        wc_box(b,0,0,0,.62f,.62f,.62f,WR_CHARCOAL);
        wc_box(b,0,0,0,.65f,.07f,.65f,WR_ENERGY);
        for(int face=0;face<4;face++)for(int y=-1;y<=1;y+=2){
            b->part=5+face;int start=b->count;float sign=face%2?1:-1;
            if(face<2)wc_box(b,sign*.335f,y*.18f,0,.075f,.20f,.40f,WR_DARK);
            else wc_box(b,0,y*.18f,sign*.335f,.40f,.20f,.075f,WR_DARK);
            start=b->count;
            if(face<2)wc_box(b,sign*.38f,y*.18f,0,.04f,.15f,.33f,WR_SILVER);
            else wc_box(b,0,y*.18f,sign*.38f,.33f,.15f,.04f,WR_SILVER);
            for(int v=start;v<b->count;v++){b->colors[v*4]=132;b->colors[v*4+1]=65;b->colors[v*4+2]=205;}
        }b->part=0;
        for(int x=-1;x<=1;x+=2)for(int y=-1;y<=1;y+=2)for(int z=-1;z<=1;z+=2)
            wc_box(b,x*.28f,y*.28f,z*.28f,.12f,.12f,.12f,WR_SILVER);
        wc_detail(b,0,.32f,0,.20f,.035f,.20f,WR_DARK);
        wc_detail(b,0,.345f,0,.10f,.025f,.10f,WR_ENERGY);
        break;
    case WEAPON_RAILGUN:
        wc_box(b,0,0,.05f,.46f,.40f,.85f,WR_CHARCOAL);
        for(int side=-1;side<=1;side+=2){
            wc_box(b,side*.23f,0,.10f,.09f,.36f,.65f,WR_SILVER);
            wc_detail(b,side*.28f,0,.12f,.018f,.11f,.35f,WR_DARK);
            wc_detail(b,side*.29f,0,.20f,.02f,.07f,.10f,WR_ENERGY);
        }
        b->part=1;
        for(int side=-1;side<=1;side+=2){
            wc_box(b,0,side*.23f,-1.05f,.27f,.17f,1.65f,WR_SILVER);
            wc_detail(b,0,side*.137f,-1.05f,.20f,.025f,1.53f,WR_DARK);
            for(int n=0;n<6;n++){
                wc_detail(b,0,side*.12f,-.39f-n*.27f,.15f,.02f,.16f,WR_ENERGY);
                for(int x=-1;x<=1;x+=2)wc_detail(b,x*.14f,side*.23f,-.39f-n*.27f,.025f,.09f,.09f,WR_DARK);
            }
        }b->part=0;
        wc_detail(b,0,0,.49f,.07f,.21f,.025f,WR_ENERGY);
        wc_box(b,0,.33f,.08f,.19f,.14f,.38f,WR_DARK);
        wr_tube(b,(Vector3){0,.40f,-.03f},.22f,.22f,.20f,.035f,.015f,WR_SILVER);
        wc_detail(b,0,.40f,.06f,.09f,.09f,.015f,WR_ENERGY);
        for(int n=-1;n<=1;n+=2)wc_detail(b,0,n*.17f,.68f,.19f,.08f,.46f,WR_DARK);
        wc_box(b,0,0,.91f,.27f,.46f,.13f,WR_SILVER);
        break;
    case WEAPON_ROCKET_LAUNCHER:
        wr_tube(b,(Vector3){0,0,-.27f},.76f,.76f,1.83f,.09f,.06f,WR_CHARCOAL);
        for(int n=0;n<3;n++)wr_tube(b,(Vector3){0,0,-.87f+n*.58f},.84f,.84f,.32f,.055f,.045f,WR_SILVER);
        wr_tube(b,(Vector3){0,0,-1.19f},.89f,.89f,.18f,.13f,.045f,WR_SILVER);
        wr_tube(b,(Vector3){0,0,.68f},.82f,.82f,.20f,.10f,.04f,WR_DARK);
        wr_tube(b,(Vector3){0,0,.79f},.70f,.70f,.08f,.06f,.02f,WR_SILVER);
        wc_handle(b,.67f,-.05f,.65f);
        for(int side=-1;side<=1;side+=2){
            wc_detail(b,side*.405f,.06f,.06f,.018f,.16f,.24f,WR_DARK);
            wc_detail(b,side*.42f,.06f,.06f,.018f,.055f,.17f,WR_ENERGY);
        }
        b->part=4;{int start=b->count;wc_rocket(b);for(int v=start;v<b->count;v++)b->vertices[v*3+2]-=.72f;}b->part=0;
        wc_flash(b,true);
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
        wc_box(b,0,0,.04f,.58f,.48f,.85f,WR_CHARCOAL);
        wc_box(b,0,.23f,.06f,.62f,.13f,.82f,WR_SILVER);
        wc_box(b,0,-.23f,.06f,.62f,.13f,.82f,WR_SILVER);
        wr_tube(b,(Vector3){0,0,.50f},.58f,.43f,.10f,.05f,.02f,WR_SILVER);
        for(int n=-1;n<=1;n++)wc_detail(b,n*.11f,0,.50f,.055f,.13f,.025f,WR_ENERGY);
        wr_tube(b,(Vector3){0,0,-.44f},.38f,.34f,.28f,.065f,.025f,WR_DARK);
        wc_box(b,0,0,-.30f,.22f,.20f,.04f,WR_ENERGY);
        for(int side=-1;side<=1;side+=2){
            wc_box(b,side*.32f,0,-.31f,.16f,.48f,.20f,WR_SILVER);
            for(int n=0;n<3;n++)wc_detail(b,side*.30f,(n-1)*.13f,.16f,.018f,.07f,.22f,WR_ENERGY);
            b->part=side<0?9:10;
            for(int top=-1;top<=1;top+=2){
                wc_box(b,side*.53f,top*.21f,-.40f,.43f,.14f,.18f,WR_SILVER);
                wc_box(b,side*.72f,top*.21f,-.65f,.14f,.17f,.55f,WR_SILVER);
                wc_detail(b,side*.72f,top*.30f,-.65f,.09f,.025f,.30f,WR_DARK);
                wc_box(b,side*.72f,top*.21f,-.93f,.18f,.20f,.12f,WR_CHARCOAL);
                wc_detail(b,side*.72f,top*.21f,-1.0f,.12f,.12f,.025f,WR_ENERGY);
            }b->part=0;
        }
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
