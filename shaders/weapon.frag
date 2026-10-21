#version 330
in vec3 position;
in vec3 normal;
in vec4 color;
in float energy;
in vec3 localPosition;
in vec3 edgeDistance;
flat in int edgeMask;
uniform vec3 eyePosition;
uniform vec3 energyColor;
uniform float energyStrength;
uniform int emissionOnly;
uniform float effectTime;
uniform int weaponKind;
uniform float charge;
uniform float activity;
uniform float viewHide;
uniform float cooldown;
uniform float deployment;
uniform float buildPulse;
uniform float shotAge;
flat in int partId;
in float visibility;
out vec4 finalColor;
void main() {
    if(visibility<.005)discard;
    if (energy > 2.5) {
        float t = localPosition.z;
        float pulse = pow(0.5 + 0.5 * sin(t * 24.0 - effectTime * 9.0 + localPosition.x), 3.0);
        float fade = smoothstep(0.0, 0.12, t) * (1.0 - smoothstep(0.86, 1.0, t));
        float core = 1.0 - smoothstep(0.20, 1.0, abs(localPosition.y));
        finalColor = vec4(mix(vec3(1.0, 0.55, 0.025), vec3(1.0, 0.94, 0.50), core * pulse),
                          fade * core * (0.60 + 0.40 * pulse));
        return;
    }
    vec3 base = color.rgb;
    bool polished=weaponKind>=2 && weaponKind!=7;
    if(polished){
        // Soft sky fill and restrained metal highlights preserve flat bevels.
        vec3 n=normalize(normal),light=normalize(vec3(-.45,.8,-.45));
        float diffuse=max(0.0,dot(n,light));
        vec3 view=normalize(eyePosition-position);
        float spec=pow(max(0.0,dot(n,normalize(light+view))),28.0);
        float metal=smoothstep(.23,.62,max(base.r,max(base.g,base.b)));
        base=base*(.62+.43*diffuse)+vec3(.15,.17,.20)*spec*metal;
    }
    // A roughly one-pixel, antialiased line along real polygon edges only.
    // Derivatives keep the width consistent in first-person and split-screen views.
    vec3 coverage = smoothstep(vec3(0.0), max(fwidth(edgeDistance), vec3(0.000001)), edgeDistance);
    float interior = 1.0;
    if ((edgeMask & 1) != 0) interior = min(interior, coverage.x);
    if ((edgeMask & 2) != 0) interior = min(interior, coverage.y);
    if ((edgeMask & 4) != 0) interior = min(interior, coverage.z);
    base = mix(vec3(0.035, 0.045, 0.065), base, interior);
    // Cubes keep subtle face tints and outlined edges while emitting energy.
    float energyEdge = energy > 1.5 ? mix(0.45, 1.0, interior) : 1.0;
    vec3 emissive = energyColor * color.r * energyStrength * energyEdge;
    // Distinct readable signals; opaque color pass keeps these visible without bloom.
    float signal=1.0;
    if(weaponKind==2)signal=.30+.65*(1.0-cooldown)+.6*max(0.0,1.0-shotAge/.12);
    if(weaponKind==3)signal=.65+activity*(.35+.65*pow(.5+.5*sin(-localPosition.z*24.0-effectTime*20.0),4.0));
    if(weaponKind==4){ // Adhesive fuse accelerates with supplied fuse progress.
        float blink=step(.45,fract(effectTime*(1.0+charge*7.0)));
        signal=mix(.7,.35+blink*1.1,activity);
        emissive=mix(emissive,vec3(1,.04,.01),charge*activity*.85);
    }
    if(weaponKind==5){ // Rear-to-tip rail bank filling, not random flicker.
        float bank=clamp((-localPosition.z-.35)/1.4,0.0,1.0);
        signal=.22+1.1*(1.0-smoothstep(charge,charge+.06,bank));
        signal+=2.0*max(0.0,1.0-shotAge/.18)*exp(-90.0*pow(bank-shotAge*8.0,2.0));
    }
    if(weaponKind==6)signal=.4+.6*(1.0-cooldown);
    if(weaponKind==7)signal=.8+.2*pow(.5+.5*sin(localPosition.z*12.0+effectTime*6.0),3.0)+activity*.4;
    if(weaponKind==8){
        signal=.45+deployment*.35;
        float travel=(1.0-buildPulse)*1.9;
        signal+=2.0*buildPulse*exp(-35.0*pow(-localPosition.z+.30-travel,2.0));
    }
    if(weaponKind==9){
        signal=.75+.18*sin(effectTime*4.0-localPosition.z*8.0)+activity*.20;
        if(partId==14){
            // Analytic bubbles appear beneath the glossy gel; no particles allocated.
            float bubbles=0.0;
            for(int i=0;i<9;i++){
                float f=float(i);
                vec2 center=vec2(-.31+fract(f*.618)*.52,-.14+fract(effectTime*.15+f*.371)*.29);
                float d=length((localPosition.zy-center)*vec2(1,1.5));
                bubbles+=1.0-smoothstep(.018,.030,d);
            }
            float spec=pow(max(0.0,dot(normalize(normal),normalize(normalize(eyePosition-position)+normalize(vec3(-.5,1,-.6))))),48.0);
            emissive=energyColor*(.45+.3*bubbles)*energyStrength;
            if(emissionOnly==0)emissive+=vec3(.5,.7,.3)*spec;
        }
    }
    if(weaponKind==10){signal=.7+.25*sin(effectTime*(2.0+deployment*6.0));if(partId==32)base*=.12;}
    if(weaponKind==12)signal=.6+deployment*.6;
    if(weaponKind==13)signal=.7+.2*activity;
    if(partId==15 || partId==16){
        vec3 n=normalize(normal),view=normalize(eyePosition-position);
        float rim=pow(1.0-abs(dot(n,view)),3.0);
        float highlight=pow(max(0.0,abs(dot(n,normalize(view+normalize(vec3(-.4,1,-.5)))))),48.0);
        finalColor=vec4(mix(vec3(.08,.26,.42),vec3(.7,.9,1),highlight),.20+.20*rim+.25*highlight);return;
    }
    if(partId==11 || partId==12){
        // Fired energy completes even if the last Matter paid for the shot.
        finalColor=vec4(energyColor*(1.2+visibility),visibility);return;
    }
    emissive*=signal;
    if (emissionOnly != 0) {
        finalColor = vec4(energy > 0.5 ? emissive : vec3(0.0), 1.0);
    } else {
        finalColor = vec4(energy > 0.5 ? emissive + vec3(0.16) : base, color.a);
    }
    finalColor.rgb*=1.0-viewHide;
    finalColor.a*=1.0-viewHide;
}
