#version 330
in vec3 normal;
in vec4 color;
uniform vec4 colDiffuse;
uniform int emissionOnly;
uniform int luminous;
uniform vec4 statusTint;
uniform float viewHide;
out vec4 finalColor;
void main() {
    vec3 base = color.rgb * colDiffuse.rgb;
    float light = 0.48 + 0.52 * max(dot(normalize(normal), normalize(vec3(-0.4,0.8,-0.5))), 0.0);
    vec3 lit = luminous != 0 ? base * 1.35 : base * light;
    lit = mix(lit, statusTint.rgb, statusTint.a * 0.65);
    finalColor = vec4(emissionOnly != 0 ? (luminous != 0 ? base : vec3(0.0)) : lit, color.a * colDiffuse.a);
    finalColor.rgb*=1.0-viewHide;
    finalColor.a*=1.0-viewHide;
}
