#version 330
in vec2 fragTexCoord;
in vec4 fragColor;
out vec4 finalColor;

in vec3 caCoord;
in vec3 caFace;
uint voxelPatternKey() {
    // Bias toward the interior so both faces of a cell use its grid identity.
    ivec3 cell=ivec3(floor(caCoord-caFace*0.001));
    uvec3 c=uvec3(cell);
    return (c.x*73856093u ^ c.y*19349663u ^ c.z*83492791u)&0xffffffu;
}
uniform sampler2D texture0;
uniform int caEnabled;
uniform float caBlend;
// Keep this sampling function identical in both voxel fragment shaders.
float cellularShade() {
    if (caEnabled == 0) return 1.0;
    vec3 n = abs(caFace);
    int axis = n.x > 0.5 ? 0 : (n.y > 0.5 ? 1 : 2);
    vec2 face = axis == 0 ? caCoord.yz : (axis == 1 ? caCoord.xz : caCoord.xy);
    float footprint = max(fwidth(face.x * 32.0), fwidth(face.y * 32.0));
    float detail = 1.0 - smoothstep(0.5, 2.0, footprint);
    if (detail <= 0.0) return 1.0;
    uvec3 rgb = uvec3(round(fragColor.rgb * 255.0));
    uint h = rgb.x * 73856093u ^ rgb.y * 19349663u ^ rgb.z * 83492791u;
    h ^= uint(axis * 2 + (caFace[axis] > 0.0 ? 1 : 0)) * 2654435761u;
    uint key=voxelPatternKey();
    key^=key>>16u;key*=0x7feb352du;key^=key>>15u;key*=0x846ca68bu;key^=key>>16u;
    h^=key;
    vec2 uv = fract(face);
    for (uint r=0u;r<(h & 3u);r++) uv=vec2(1.0-uv.y,uv.x);
    ivec2 tile=ivec2(int((h>>2u)&1u),int((h>>3u)&1u))*128;
    ivec2 offset=ivec2(int((h>>4u)&127u),int((h>>11u)&127u));
    ivec2 cell=(ivec2(floor(uv*32.0))+offset)%128;
    vec3 states=texelFetch(texture0,tile+cell,0).rgb;
    float live=mix(states.g,states.r,caBlend);
    float trail=max(0.0,states.b-caBlend*(16.0/255.0));
    return 1.0+detail*(-0.10+0.20*live+0.08*trail);
}

void main() {
    float u = fragTexCoord.x;
    float v = fragTexCoord.y;
    // Greedy mesh UVs are measured in voxel cells, matching the per-cube
    // edge pattern in instanced_voxel_hack.frag without splitting geometry.
    float fu = fract(u);
    float fv = fract(v);
    bool isEdge = (fu < 0.02 || fu > 0.98 || fv < 0.02 || fv > 0.98);
    
    vec3 result = fragColor.rgb * cellularShade();
    if (isEdge) result *= 0.1;
    finalColor = vec4(result, fragColor.a);
}
