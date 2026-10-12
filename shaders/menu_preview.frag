#version 330
in vec3 normal;
uniform vec4 colDiffuse;
out vec4 finalColor;
void main() {
    // Key and fill face the showcase camera so both visible sides retain detail.
    vec3 n = normalize(normal);
    float light = 0.30
        + 0.60 * max(dot(n, normalize(vec3(0.8, 1.0, 0.3))), 0.0)
        + 0.18 * max(dot(n, normalize(vec3(-0.3, 0.2, 1.0))), 0.0);
    finalColor = vec4(colDiffuse.rgb * light, colDiffuse.a);
}
