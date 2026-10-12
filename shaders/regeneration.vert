#version 330
in vec3 vertexPosition;
in vec3 vertexNormal;
in vec4 vertexColor;
uniform mat4 mvp;
out vec3 worldPosition;
out vec3 normal;
out vec4 color;
void main() {
    worldPosition=vertexPosition;
    normal=vertexNormal;
    color=vertexColor;
    gl_Position=mvp*vec4(vertexPosition,1.0);
}
