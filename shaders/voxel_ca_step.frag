#version 330
uniform sampler2D previousState;
out vec4 finalColor;
// Default second-order outer-totalistic rule from xcont.com/cell/single_2nd.js.
// Index = 18 * previous center + 9 * current center + live neighbor count.
const int rule[36]=int[36](
    0,0,1,0,0,0,0,1,0, 0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,1,0,1, 0,1,1,0,0,0,0,0,0);
void main() {
    ivec2 p=ivec2(gl_FragCoord.xy),base=(p/128)*128,local=p%128;
    vec3 old=texelFetch(previousState,p,0).rgb;
    int current=int(round(old.r)),previous=int(round(old.g)),neighbors=0;
    for(int y=-1;y<=1;y++)for(int x=-1;x<=1;x++) {
        if(x==0 && y==0)continue;
        ivec2 q=base+(local+ivec2(x,y)+128)%128;
        neighbors+=int(round(texelFetch(previousState,q,0).r));
    }
    int next=rule[previous*18+current*9+neighbors];
    // A short activity trail uses the spare blue channel; alpha stays opaque.
    float trail=next!=current?1.0:max(0.0,old.b-16.0/255.0);
    finalColor=vec4(float(next),float(current),trail,1.0);
}
