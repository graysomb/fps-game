#version 330
in vec3 normal;
in vec2 faceUV;
in vec4 color;
flat in float assembly;
uniform int renderMode;
uniform int luminous;
uniform float clockTime;
uniform float hitFlash;
uniform float effectOpacity;
out vec4 finalColor;
void main(){
    vec2 width=max(fwidth(faceUV)*1.15,vec2(.0001));
    vec2 d=min(faceUV,1.0-faceUV);
    float interior=min(smoothstep(0.0,width.x,d.x),smoothstep(0.0,width.y,d.y));
    if(renderMode>=2){
        if(assembly>=.999 || interior>.65)discard;
        float pulse=.75+.25*sin(clockTime*5.0);
        finalColor=vec4(vec3(.15,.8,1.0)*pulse,(1.0-interior)*.52*effectOpacity);return;
    }
    if(assembly<=.001)discard;
    float light=.48+.52*max(dot(normalize(normal),normalize(vec3(-.45,.85,-.6))),0.0);
    vec3 base=mix(vec3(.025,.035,.045),color.rgb*light,interior);
    vec3 energy=luminous!=0?color.rgb*(.92+.08*sin(clockTime*2.4)):vec3(0.0);
    if(luminous!=0)base*=.18;
    energy+=vec3(.06,.65,.9)*(1.0-assembly)*interior;
    base=mix(base,vec3(1.0,.59,.22),hitFlash*.7);
    vec3 result=renderMode==1?energy:base+energy;
    finalColor=vec4(result,effectOpacity);
}
