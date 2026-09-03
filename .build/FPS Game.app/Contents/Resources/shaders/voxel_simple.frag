#version 330
in vec2 fragTexCoord;
in vec4 fragColor;
out vec4 finalColor;
void main() {
    float edgeWidth = 0.02;
    float u = fragTexCoord.x;
    float v = fragTexCoord.y;
    // For greedy mesh, UVs might be larger than 1.0 if it's tiled, 
    // but gen_greedy_mesh usually gives 0-1 per quad.
    // We use fract() to handle potential tiling.
    float fu = fract(u);
    float fv = fract(v);
    bool isEdge = (fu < 0.02 || fu > 0.98 || fv < 0.02 || fv > 0.98);
    
    vec3 result = fragColor.rgb;
    if (isEdge) result *= 0.1;
    finalColor = vec4(result, fragColor.a);
}
