#version 330

in vec3 vertexPosition;
in vec3 vertexNormal;
uniform mat4 mvp;
out vec3 localPosition;
out vec3 localNormal;

void main() {
    localPosition = vertexPosition;
    localNormal = vertexNormal;
    gl_Position = mvp * vec4(vertexPosition, 1.0);
}
