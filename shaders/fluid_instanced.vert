#version 330

in vec3 vertexPosition;
in vec3 vertexNormal;
// xyz = world position, w = uniform scale applied to the sphere mesh
in vec4 instancePosRadius;

uniform mat4 matView;
uniform mat4 matProjection;

out vec3 fragPosition;
out vec3 fragNormal;

void main()
{
    vec4 worldPosition = vec4(vertexPosition * instancePosRadius.w + instancePosRadius.xyz, 1.0);
    fragPosition = worldPosition.xyz;
    fragNormal = normalize(vertexNormal);
    gl_Position = matProjection * matView * worldPosition;
}
