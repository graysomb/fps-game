#version 330
in vec3 position;
in vec3 normal;
in vec4 color;
in float energy;
in vec3 localPosition;
uniform vec3 eyePosition;
uniform vec3 energyColor;
uniform float energyStrength;
uniform int emissionOnly;
out vec4 finalColor;
void main() {
    vec3 N = normalize(normal);
    vec3 V = normalize(eyePosition - position);
    vec3 L = normalize(vec3(-0.4, 0.85, 0.55));
    float diffuse = max(dot(N, L), 0.0);
    float hemisphere = N.y * 0.5 + 0.5;
    float specular = pow(max(dot(N, normalize(L + V)), 0.0), 24.0);
    vec3 base = color.rgb * (mix(vec3(0.24, 0.28, 0.38), vec3(0.55, 0.59, 0.68), hemisphere)
                             + vec3(0.63) * diffuse);
    base += vec3(0.16) * specular;
    // Local spill is restricted to the gun, independent of world lighting.
    float spill = exp(-5.0 * length(localPosition - vec3(0.0, 0.0, -0.38)));
    base += energyColor * spill * 0.13 * energyStrength;
    vec3 emissive = energyColor * color.r * energyStrength;
    if (emissionOnly != 0) {
        finalColor = vec4(energy > 0.5 ? emissive : vec3(0.0), 1.0);
    } else {
        finalColor = vec4(energy > 0.5 ? emissive + vec3(0.16) : base, 1.0);
    }
}
