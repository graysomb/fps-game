#version 330
in vec4 fragColor;
uniform vec4 colDiffuse;
uniform float effectTime;
out vec4 finalColor;
void main() {
    float bands=0.5+0.5*sin((gl_FragCoord.x+gl_FragCoord.y)*0.16-effectTime*35.0);
    vec4 tint=fragColor*colDiffuse;
    finalColor=vec4(mix(tint.rgb,vec3(1.0),0.25+0.25*bands),tint.a);
}
