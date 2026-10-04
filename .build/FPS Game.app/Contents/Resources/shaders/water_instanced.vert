#version 330

in vec3 vertexPosition;
in vec3 vertexNormal;
in mat4 instanceTransform;

uniform mat4 matView;
uniform mat4 matProjection;

out vec3 fragNormal;
out vec3 fragPosition;

void main()
{
    vec4 worldPos = instanceTransform * vec4(vertexPosition, 1.0);
    fragPosition = worldPos.xyz;
    fragNormal = normalize((instanceTransform * vec4(vertexNormal, 0.0)).xyz);
    gl_Position = matProjection * matView * worldPos;
}
