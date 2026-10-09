/* Deterministic arena in the same grid units/material palette as the citadel. */
static void buildFoundryHoldout(void) {
    Color stone={82,102,117,255}, trim={127,173,178,255}, support={235,155,45,255};
    /* Raised holdout, with stairs on north and south and an open lower ring. */
    add_static_box_at_grid(-10,10,0,3,-8,8,stone);
    for(int step=0;step<4;step++) {
        add_static_box_at_grid(-4,4,0,step,9+(3-step)*2,10+(3-step)*2,trim);
        add_static_box_at_grid(-4,4,0,step,-10-(3-step)*2,-9-(3-step)*2,trim);
    }
    add_static_box_at_grid(-10,-9,4,5,-5,-1,trim);
    add_static_box_at_grid(9,10,4,5,1,5,trim);
    /* Heavy slabs at east/west choke points. Two colored supports each. */
    for(int side=-1;side<=1;side+=2) {
        int x=side*23;
        add_static_box_at_grid(x-2,x+2,7,8,-5,5,stone);
        add_static_box_at_grid(x,x,0,6,-5,-5,support);
        add_static_box_at_grid(x,x,0,6,5,5,support);
        /* Loose quarry chunks, separated from supports and alternate routes. */
        for(int z=-10;z<=10;z+=20) add_static_box_at_grid(x-2,x,0,1,z,z+1,trim);
        /* Gate pylons frame entrances without closing the lower ring. */
        add_static_box_at_grid(side*36-1,side*36+1,0,9,-12,-10,stone);
        add_static_box_at_grid(side*36-1,side*36+1,0,9,10,12,stone);
    }
    /* Nearby blocks make the first tutorial harvest reachable. */
    add_static_box_at_grid(-3,3,4,5,-5,-5,support);
}
