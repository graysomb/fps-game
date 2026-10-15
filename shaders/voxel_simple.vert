#version 330
in vec3 vertexPosition;
in vec2 vertexTexCoord;
in vec4 vertexColor;
in vec3 vertexNormal;
uniform mat4 mvp;
out vec2 fragTexCoord;
out vec4 fragColor;
out vec3 caCoord;
out vec3 caFace;
void main() {
    fragTexCoord = vertexTexCoord;
    fragColor = vertexColor;
    caCoord = vertexPosition / 0.5;
    caFace = vertexNormal;
    gl_Position = mvp * vec4(vertexPosition, 1.0);
}
