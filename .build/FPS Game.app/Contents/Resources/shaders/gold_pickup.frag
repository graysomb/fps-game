#version 330

in vec3 localPosition;
in vec3 localNormal;
uniform float uTime;
out vec4 finalColor;

void main() {
    float bands = 0.5 + 0.5 * sin(localPosition.y * 22.0 - uTime * 5.0);
    float pulse = 0.5 + 0.5 * sin(uTime * 4.0);
    float rim = pow(1.0 - abs(normalize(localNormal).z), 2.0);
    vec3 gold = mix(vec3(0.48, 0.20, 0.015), vec3(1.0, 0.78, 0.18),
                    0.28 + 0.38 * bands + 0.22 * pulse);
    finalColor = vec4(gold + vec3(0.28, 0.16, 0.02) * rim, 1.0);
}
