#version 330
in vec2 fragTexCoord;
in vec4 fragColor;
out vec4 finalColor;
void main() {
    float u = fragTexCoord.x;
    float v = fragTexCoord.y;
    // Greedy mesh UVs are measured in voxel cells, matching the per-cube
    // edge pattern in instanced_voxel_hack.frag without splitting geometry.
    float fu = fract(u);
    float fv = fract(v);
    bool isEdge = (fu < 0.02 || fu > 0.98 || fv < 0.02 || fv > 0.98);
    
    vec3 result = fragColor.rgb;
    if (isEdge) result *= 0.1;
    finalColor = vec4(result, fragColor.a);
}
