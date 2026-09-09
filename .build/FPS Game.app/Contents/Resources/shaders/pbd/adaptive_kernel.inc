/* Portable scalar-vector kernels. Included by C and shader wrappers.
 * No hidden particle ownership and no atomic scatter: each job owns a cell. */
AD_FN AdVec ad_vec(float x, float y, float z) {
    AdVec r; r.x=x; r.y=y; r.z=z; r.w=0.0; return r;
}
AD_FN AdVec ad_add(AdVec a, AdVec b) { return ad_vec(a.x+b.x,a.y+b.y,a.z+b.z); }
AD_FN AdVec ad_sub(AdVec a, AdVec b) { return ad_vec(a.x-b.x,a.y-b.y,a.z-b.z); }
AD_FN AdVec ad_mul(AdVec a, float s) { return ad_vec(a.x*s,a.y*s,a.z*s); }
AD_FN float ad_dot(AdVec a, AdVec b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
AD_FN AdVec ad_cross(AdVec a, AdVec b) { return ad_vec(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x); }
AD_FN float ad_length(AdVec a) { return sqrt(ad_dot(a,a)); }
AD_FN AdVec ad_normalize(AdVec a) { return ad_mul(a,1.0/AD_MAX(ad_length(a),0.000001)); }
AD_FN float ad_limit(float x,float lo,float hi) { return x<lo?lo:(x>hi?hi:x); }
AD_FN float ad_component(AdVec a, int k) { return k==0?a.x:(k==1?a.y:a.z); }
AD_FN AdVec ad_center(AdCell c) {
    AdVec p=ad_vec(0.0,0.0,0.0);
    for(int k=0;k<8;k++) p=ad_add(p,c.x[k]);
    return ad_mul(p,0.125);
}
AD_FN AdVec ad_axis(AdCell c, int axis) {
    AdVec a=ad_vec(0.0,0.0,0.0);
    int bit=1<<axis;
    for(int k=0;k<8;k++) if((k&bit)==0) a=ad_add(a,ad_sub(c.x[k|bit],c.x[k]));
    return ad_mul(a,0.25);
}
AD_FN float ad_weight(AdVec uv, int corner) {
    return ((corner&1)!=0?uv.x:1.0-uv.x)*((corner&2)!=0?uv.y:1.0-uv.y)*((corner&4)!=0?uv.z:1.0-uv.z);
}
AD_FN AdVec ad_sample(AdCell c, AdVec uv) {
    AdVec p=ad_vec(0.0,0.0,0.0);
    for(int k=0;k<8;k++) p=ad_add(p,ad_mul(c.x[k],ad_weight(uv,k)));
    return p;
}
AD_FN AdVec ad_local(AdCell c, AdVec rest) {
    return ad_vec((rest.x-float(c.grid[0]))/float(c.grid[3]),
                  (rest.y-float(c.grid[1]))/float(c.grid[3]),
                  (rest.z-float(c.grid[2]))/float(c.grid[3]));
}
AD_FN AdVec ad_rest_corner(AdCell c,int k) {
    return ad_vec(float(c.grid[0]+((k&1)!=0?c.grid[3]:0)),
                  float(c.grid[1]+((k&2)!=0?c.grid[3]:0)),
                  float(c.grid[2]+((k&4)!=0?c.grid[3]:0)));
}
AD_FN AdCell ad_integrate_cell(AdCell c,float dt,float gravity) {
    for(int k=0;k<8;k++) {
        c.old[k]=c.x[k];
        if(c.v[k].w>0.0 && (c.meta[3]&6u)==0u) {
            c.v[k].y+=gravity*dt;
            c.x[k]=ad_add(c.x[k],ad_mul(c.v[k],dt));
        }
    }
    return c;
}
AD_FN AdCell ad_shape_cell(AdCell c,float unit) {
    if((c.meta[3]&7u)!=0u) return c;
    AdVec center=ad_center(c);
    float edge=unit*float(c.grid[3]);
    AdVec a=ad_axis(c,0),b=ad_axis(c,1),d=ad_axis(c,2);
    AdVec initial_a=a,initial_b=b,initial_d=d;
    /* Symmetric Gram-Schmidt relaxation with a volume correction. */
    for(int it=0;it<3;it++) {
        AdVec na=ad_sub(a,ad_mul(ad_add(ad_mul(b,ad_dot(a,b)/AD_MAX(ad_dot(b,b),0.000001)),ad_mul(d,ad_dot(a,d)/AD_MAX(ad_dot(d,d),0.000001))),0.5));
        AdVec nb=ad_sub(b,ad_mul(ad_add(ad_mul(a,ad_dot(b,a)/AD_MAX(ad_dot(a,a),0.000001)),ad_mul(d,ad_dot(b,d)/AD_MAX(ad_dot(d,d),0.000001))),0.5));
        AdVec nd=ad_sub(d,ad_mul(ad_add(ad_mul(a,ad_dot(d,a)/AD_MAX(ad_dot(a,a),0.000001)),ad_mul(b,ad_dot(d,b)/AD_MAX(ad_dot(b,b),0.000001))),0.5));
        a=na;b=nb;d=nd;
    }
    a=ad_mul(ad_normalize(a),ad_limit(0.1*edge+0.9*ad_length(a),0.8*edge,1.2*edge));
    b=ad_mul(ad_normalize(b),ad_limit(0.1*edge+0.9*ad_length(b),0.8*edge,1.2*edge));
    d=ad_mul(ad_normalize(d),ad_limit(0.1*edge+0.9*ad_length(d),0.8*edge,1.2*edge));
    float volume=ad_dot(ad_cross(a,b),d);
    if(volume<edge*edge*edge*0.1) {
        /* Recover an inverted or collapsed cell instead of leaving an
           unsolvable shape to accumulate arbitrarily large deformation. */
        a=ad_length(initial_a)>0.000001?ad_normalize(initial_a):ad_vec(1.0,0.0,0.0);
        b=ad_sub(initial_b,ad_mul(a,ad_dot(a,initial_b)));
        if(ad_length(b)<0.000001)b=ad_cross(initial_d,a);
        if(ad_length(b)<0.000001)b=ad_cross(a,(a.y<0.0?-a.y:a.y)<0.9?ad_vec(0.0,1.0,0.0):ad_vec(1.0,0.0,0.0));
        b=ad_normalize(b);d=ad_cross(a,b);
        a=ad_mul(a,edge);b=ad_mul(b,edge);d=ad_mul(d,edge);
        volume=edge*edge*edge;
    }
    float scale=pow(edge*edge*edge/volume,1.0/3.0);
    a=ad_mul(a,scale*0.5);b=ad_mul(b,scale*0.5);d=ad_mul(d,scale*0.5);
    for(int k=0;k<8;k++) if(c.v[k].w>0.0)
        c.x[k]=ad_add(center,ad_add(ad_mul(a,(k&1)!=0?1.0:-1.0),ad_add(ad_mul(b,(k&2)!=0?1.0:-1.0),ad_mul(d,(k&4)!=0?1.0:-1.0))));
    return c;
}
AD_FN AdCell ad_solve_links(int index, AD_LINK_PARAMS) {
    AdCell result=AD_CELL(index);
    AdVec delta[8];
    for(int k=0;k<8;k++){delta[k]=ad_vec(0.0,0.0,0.0);}
    AdRange range=AD_RANGE(index);
    for(int r=0;r<range.count;r++) {
        AdLink link=AD_LINK(AD_ID(range.begin+r));
        AdCell a=AD_CELL(link.a),b=AD_CELL(link.b);
        bool mineA=index==link.a;
        bool smallA=a.grid[3]<=b.grid[3];
        AdCell fine=a;
        if(!smallA) fine=b;
        int bit=1<<link.axis;
        int side=smallA?link.side:1-link.side;
        for(int c=0;c<8;c++) if(((c&bit)!=0?1:0)==side) {
            AdVec rest=ad_rest_corner(fine,c),ua=ad_local(a,rest),ub=ad_local(b,rest);
            AdVec error=ad_sub(ad_sample(a,ua),ad_sample(b,ub));
            float denom=0.0;
            for(int k=0;k<8;k++) {float wa=ad_weight(ua,k),wb=ad_weight(ub,k);denom+=wa*wa*a.v[k].w+wb*wb*b.v[k].w;}
            if(denom<=0.0) continue;
            for(int k=0;k<8;k++) {
                float wgt=mineA?ad_weight(ua,k):ad_weight(ub,k);
                if(wgt<=0.0 || result.v[k].w<=0.0) continue;
                delta[k]=ad_add(delta[k],ad_mul(error,(mineA?-1.0:1.0)*result.v[k].w*wgt*link.relaxation/denom));
            }
        }
    }
    for(int k=0;k<8;k++) result.x[k]=ad_add(result.x[k],delta[k]);
    return result;
}
AD_FN AdCell ad_solve_contacts(int index, AD_CONTACT_PARAMS) {
    AdCell result=AD_CELL(index); AdRange range=AD_CRANGE(index);
    AdVec delta[8];
    for(int k=0;k<8;k++){delta[k]=ad_vec(0.0,0.0,0.0);}
    for(int r=0;r<range.count;r++) {
        AdContact ct=AD_CONTACT(AD_CID(range.begin+r));
        AdCell a=AD_CELL(ct.a), b=a;
        if(ct.b>=0) b=AD_CELL(ct.b);
        AdVec pa=ad_vec(0.0,0.0,0.0),pb=ct.pb;
        if(ct.b>=0) pb=ad_vec(0.0,0.0,0.0);
        float denom=0.0;
        for(int k=0;k<8;k++) {
            pa=ad_add(pa,ad_mul(a.x[k],ct.wa[k]));denom+=ct.wa[k]*ct.wa[k]*a.v[k].w;
            if(ct.b>=0){pb=ad_add(pb,ad_mul(b.x[k],ct.wb[k]));denom+=ct.wb[k]*ct.wb[k]*b.v[k].w;}
        }
        float gap=ad_dot(ad_sub(pa,pb),ct.normal)+ct.normal.w;
        if(gap>=0.0 || denom<=0.0) continue;
        float relaxation=1.0;
        if(ct.b>=0)relaxation=1.0/AD_MAX(1.0,AD_MAX(AD_CRANGE(ct.a).count,AD_CRANGE(ct.b).count));
        for(int k=0;k<8;k++) {
            float wt=index==ct.a?ct.wa[k]:-ct.wb[k];
            if(wt==0.0 || result.v[k].w<=0.0) continue;
            delta[k]=ad_add(delta[k],ad_mul(ct.normal,-gap*wt*result.v[k].w*relaxation/denom));
        }
    }
    for(int k=0;k<8;k++) result.x[k]=ad_add(result.x[k],delta[k]);
    return result;
}
