#version 330
in vec3 worldPosition;
in vec3 normal;
in vec4 color;
uniform float fillHeight;
out vec4 finalColor;
void main() {
    float filled=1.0-smoothstep(fillHeight-.06,fillHeight+.06,worldPosition.y);
    float scan=1.0-smoothstep(0.0,.10,abs(worldPosition.y-fillHeight));
    float light=.55+.45*max(dot(normalize(normal),normalize(vec3(.4,.8,.3))),0.0);
    vec3 tint=mix(vec3(.16,.65,.80),color.rgb,.5)*light;
    finalColor=vec4(tint+vec3(.3,.8,1.0)*scan,.09+.30*filled+.4*scan);
}
