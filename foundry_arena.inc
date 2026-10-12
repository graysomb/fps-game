/* Expand horizontal grid bounds while retaining half-metre vertical steps. */
static void foundry_box(int x0,int x1,int y0,int y1,int z0,int z1,Color c) {
    add_static_box_at_grid((int)floorf(x0*1.875f),(int)ceilf((x1+1)*1.875f)-1,y0,y1,
        (int)floorf(z0*1.875f),(int)ceilf((z1+1)*1.875f)-1,c);
}
/* Deterministic arena in the same grid units/material palette as the citadel. */
static void buildFoundryHoldout(void) {
    Color stone={82,102,117,255}, trim={127,173,178,255}, support={235,155,45,255};
    /* Raised holdout, with stairs on north and south and an open lower ring. */
    foundry_box(-10,10,0,3,-8,8,stone);
    for(int step=0;step<4;step++) {
        foundry_box(-4,4,0,step,9+(3-step)*2,10+(3-step)*2,trim);
        foundry_box(-4,4,0,step,-10-(3-step)*2,-9-(3-step)*2,trim);
    }
    foundry_box(-10,-9,4,5,-5,-1,trim);
    foundry_box(9,10,4,5,1,5,trim);
    /* Heavy slabs at east/west choke points. Two colored supports each. */
    for(int side=-1;side<=1;side+=2) {
        int x=side*23;
        foundry_box(x-2,x+2,7,9,-5,5,stone);
        foundry_box(x-1,x+1,0,6,-5,-4,support);
        foundry_box(x-1,x+1,0,6,4,5,support);
        /* Loose quarry chunks, separated from supports and alternate routes. */
        for(int z=-10;z<=10;z+=20) foundry_box(x-2,x,0,1,z,z+1,trim);
        /* Gate pylons frame entrances without closing the lower ring. */
        foundry_box(side*36-1,side*36+1,0,9,-12,-10,stone);
        foundry_box(side*36-1,side*36+1,0,9,10,12,stone);
    }
    /* Nearby blocks make the first tutorial harvest reachable. */
    foundry_box(-3,3,4,5,-5,-5,support);
    /* Four spacious factory bays: tall walls and roofs, two open exits each. */
    for(int sx=-1;sx<=1;sx+=2) for(int sz=-1;sz<=1;sz+=2) {
        int x=sx*23,z=sz*22;
        foundry_box(x-6,x+6,0,9,z-5,z-4,stone);
        foundry_box(x-6,x+6,0,9,z+4,z+5,stone);
        foundry_box(x-6,x+6,10,11,z-5,z+5,trim);
        foundry_box(sx*10-2,sx*10+2,0,4,sz*15-1,sz*15+1,stone);
    }
    /* Redundant roof piers leave both bay ends open. */
    for(int sx=-1;sx<=1;sx+=2) for(int sz=-1;sz<=1;sz+=2) {
        int x=sx*23,z=sz*22;
        foundry_box(x-5,x-4,0,9,z-4,z-3,stone);
        foundry_box(x+4,x+5,0,9,z+3,z+4,stone);
        foundry_box(sx*14-2,sx*14+2,0,4,sz*24-1,sz*24,stone);
        foundry_box(sx*28-1,sx*28,0,4,sz*10-2,sz*10+2,stone);
    }
    /* Staggered central cover and freestanding machinery along the flank routes. */
    foundry_box(-8,-5,4,7,2,3,trim);
    foundry_box(5,8,4,7,-3,-2,trim);
    for(int side=-1;side<=1;side+=2) {
        foundry_box(side*35-1,side*35+1,0,5,-3,3,stone);
        foundry_box(-4,4,0,5,side*33-1,side*33+1,stone);
    }
    if(!firefightPractice)enemy_stock_scrap(true);
}
