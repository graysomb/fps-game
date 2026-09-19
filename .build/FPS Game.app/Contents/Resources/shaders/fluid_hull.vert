#version 330

in vec3 vertexPosition;
in vec3 vertexNormal;
in vec4 vertexColor;

uniform mat4 mvp;

out vec3 fragNormal;
out vec4 fragColor;

void main()
{
    fragNormal = vertexNormal;
    fragColor = vertexColor;
    gl_Position = mvp * vec4(vertexPosition, 1.0);
}
