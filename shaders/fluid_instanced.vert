#version 330

in vec3 vertexPosition;
in vec3 vertexNormal;
in mat4 instanceTransform;

uniform mat4 matView;
uniform mat4 matProjection;

out vec3 fragPosition;
out vec3 fragNormal;

void main()
{
    vec4 worldPosition = instanceTransform * vec4(vertexPosition, 1.0);
    fragPosition = worldPosition.xyz;
    fragNormal = normalize((instanceTransform * vec4(vertexNormal, 0.0)).xyz);
    gl_Position = matProjection * matView * worldPosition;
}
