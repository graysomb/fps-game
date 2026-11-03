#version 330
in vec3 normal;
in vec3 localPosition;
uniform vec4 colDiffuse;
out vec4 finalColor;
void main() {
    // Key and fill face the showcase camera so both visible sides retain detail.
    vec3 n = normalize(normal);
    float light = 0.30
        + 0.60 * max(dot(n, normalize(vec3(0.8, 1.0, 0.3))), 0.0)
        + 0.18 * max(dot(n, normalize(vec3(-0.3, 0.2, 1.0))), 0.0);
    // Second nearest box face identifies true cube edges without triangulation lines.
    vec3 a=abs(localPosition);
    float edge=max(min(a.x,a.y),max(min(a.y,a.z),min(a.z,a.x)));
    float width=max(fwidth(edge),.001);
    float rim=smoothstep(.475-width,.475+width,edge);
    finalColor = vec4(colDiffuse.rgb * light * mix(1.0,.35,rim), colDiffuse.a);
}
