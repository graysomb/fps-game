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
out vec4 finalColor;
void main() {
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
    if (emissionOnly != 0) {
        finalColor = vec4(energy > 0.5 ? emissive : vec3(0.0), 1.0);
    } else {
        finalColor = vec4(energy > 0.5 ? emissive + vec3(0.16) : base, color.a);
    }
}
