#version 330

in vec2 fragTexCoord;
in vec3 fragColor;
out vec4 finalColor;

void main()
{
    // Match the original voxel's flat color and dark lines along each face edge.
    vec2 uv = fract(fragTexCoord);
    bool edge = uv.x < 0.05 || uv.x > 0.95 || uv.y < 0.05 || uv.y > 0.95;
    finalColor = vec4(edge ? fragColor * 0.1 : fragColor, 1.0);
}
