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
out vec4 finalColor;
void main() {
    if (energy > 2.5) {
        float t = localPosition.z;
        float pulse = pow(0.5 + 0.5 * sin(t * 24.0 - effectTime * 9.0 + localPosition.x), 3.0);
        float fade = smoothstep(0.0, 0.12, t) * (1.0 - smoothstep(0.86, 1.0, t));
        float core = 1.0 - smoothstep(0.20, 1.0, abs(localPosition.y));
        finalColor = vec4(mix(vec3(1.0, 0.55, 0.025), vec3(1.0, 0.94, 0.50), core * pulse),
                          fade * core * (0.60 + 0.40 * pulse));
        return;
    }
    // Match the world's unlit colors: no directional shading, highlights, or spill.
    vec3 base = color.rgb;
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
    if(weaponKind==2)signal=.8+.2*sin(effectTime*5.0)+activity*.6; // impulse chambers
    if(weaponKind==3)signal=.75+.25*sin(localPosition.z*30.0+effectTime*16.0)+activity*.5;
    if(weaponKind==4){ // Adhesive fuse accelerates with supplied fuse progress.
        float blink=step(.45,fract(effectTime*(1.0+charge*7.0)));
        signal=mix(.7,.35+blink*1.1,activity);
        emissive=mix(emissive,vec3(1,.04,.01),charge*activity*.85);
    }
    if(weaponKind==5){ // Rear-to-tip rail bank filling, not random flicker.
        float bank=clamp((-localPosition.z-.35)/1.4,0.0,1.0);
        signal=.22+1.1*(1.0-smoothstep(charge,charge+.06,bank));
    }
    if(weaponKind==6)signal=.8+.15*sin(effectTime*4.0)+activity*.7;
    if(weaponKind==7)signal=.8+.2*pow(.5+.5*sin(localPosition.z*12.0+effectTime*6.0),3.0)+activity*.4;
    if(weaponKind==8)signal=.6+.4*pow(.5+.5*sin(localPosition.y*18.0-effectTime*5.0),2.0);
    emissive*=signal;
    if (emissionOnly != 0) {
        finalColor = vec4(energy > 0.5 ? emissive : vec3(0.0), 1.0);
    } else {
        finalColor = vec4(energy > 0.5 ? emissive + vec3(0.16) : base, color.a);
    }
}
