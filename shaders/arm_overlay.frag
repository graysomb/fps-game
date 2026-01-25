#version 330

in vec2 fragTexCoord;
in vec4 fragColor;

out vec4 finalColor;

uniform float progress;
uniform float reach;
uniform vec4 armColor;

void main() {
    vec2 uv = fragTexCoord;
    float edge = min(min(uv.x, 1.0 - uv.x), min(uv.y, 1.0 - uv.y));
    float edgeMask = smoothstep(0.0, 0.12, edge);

    float pulse = 0.82 + 0.18 * sin(progress * 18.0);
    float reachBoost = 0.65 + 0.35 * clamp(reach, 0.0, 1.0);

    vec3 base = armColor.rgb;
    vec3 shaded = mix(base * 0.55, base * 1.05, edgeMask);

    float alpha = armColor.a * edgeMask * pulse * reachBoost;
    finalColor = vec4(shaded, alpha);
}
