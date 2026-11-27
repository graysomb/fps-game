#version 330

in vec3 vertexPosition;
in vec2 vertexTexCoord;
in mat4 instanceTransform;

uniform mat4 matView;
uniform mat4 matProjection;

out vec2 fragTexCoord;
out vec3 fragColor;

void main()
{
    fragColor = vec3(instanceTransform[0][3], instanceTransform[1][3], instanceTransform[2][3]);

    mat4 transform = instanceTransform;
    transform[0][3] = 0.0;
    transform[1][3] = 0.0;
    transform[2][3] = 0.0;
    transform[3][3] = 1.0;

    fragTexCoord = vertexTexCoord;
    gl_Position = matProjection * matView * transform * vec4(vertexPosition, 1.0);
}
