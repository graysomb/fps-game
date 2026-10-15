#version 330
uniform sampler2D previousState;
uniform int generation;
out vec4 finalColor;
uint hashCell(uint x) {
    x^=x>>16u;x*=0x7feb352du;x^=x>>15u;x*=0x846ca68bu;return x^(x>>16u);
}
void main() {
    ivec2 p=ivec2(gl_FragCoord.xy),base=(p/128)*128,local=p%128;
    int neighbors=0;
    for(int y=-1;y<=1;y++)for(int x=-1;x<=1;x++) {
        if(x==0 && y==0)continue;
        ivec2 q=base+(local+ivec2(x,y)+128)%128;
        neighbors+=texelFetch(previousState,q,0).r>0.5?1:0;
    }
    bool alive=texelFetch(previousState,p,0).r>0.5;
    bool next=alive?(neighbors>=1 && neighbors<=5):(neighbors==3);
    // One in 4096 cells is perturbed each step, deterministically, to sustain motion.
    if((hashCell(uint(p.y*256+p.x)^uint(generation)*747796405u)&4095u)==0u)next=!next;
    // Store both generations together so surface shading needs only one fetch.
    finalColor=vec4(next?1.0:0.0,alive?1.0:0.0,0,1);
}
