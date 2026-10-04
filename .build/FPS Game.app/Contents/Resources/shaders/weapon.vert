#version 330
in vec3 vertexPosition;
in vec3 vertexNormal;
in vec2 vertexTexCoord;
in vec4 vertexColor;
uniform mat4 mvp;
uniform mat4 matModel;
uniform mat4 matNormal;
out vec3 position;
out vec3 normal;
out vec4 color;
out float energy;
out vec3 localPosition;
void main() {
    localPosition = vertexPosition;
    position = vec3(matModel * vec4(vertexPosition, 1.0));
    normal = normalize(vec3(matNormal * vec4(vertexNormal, 0.0)));
    color = vertexColor;
    energy = vertexTexCoord.x;
    gl_Position = mvp * vec4(vertexPosition, 1.0);
}
