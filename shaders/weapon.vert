#version 330
in vec3 vertexPosition;
in vec3 vertexNormal;
in vec2 vertexTexCoord;
in vec2 vertexTexCoord2;
in vec4 vertexColor;
uniform mat4 mvp;
uniform mat4 matModel;
uniform mat4 matNormal;
uniform float effectTime;
uniform int launcher;
uniform int weaponKind;
uniform float activity;
out vec3 position;
out vec3 normal;
out vec4 color;
out float energy;
out vec3 localPosition;
out vec3 edgeDistance;
flat out int edgeMask;
void main() {
    localPosition = vertexPosition;
    vec3 p = vertexPosition;
    // Construction arms deploy laterally without moving the grip or receiver.
    if(weaponKind==8 && abs(p.x)>.35) p.x=sign(p.x)*(.35+(abs(p.x)-.35)*(.35+.65*activity));
    if (vertexTexCoord.x > 2.5) {
        float t = vertexPosition.z;
        float radius = (launcher != 0 ? 0.68 : 0.54) + 0.025 * sin(t * 31.0 - effectTime * 8.0);
        float angle = vertexPosition.x + effectTime * 1.8 +
                      0.06 * sin(t * 43.0 - effectTime * 11.0) + vertexPosition.y * 0.040 / radius;
        p = vec3(cos(angle) * radius, sin(angle) * radius * 0.85 + (launcher != 0 ? 0.14 : 0.0),
                 mix(launcher != 0 ? 0.94 : 0.70, launcher != 0 ? -1.60 : -1.02, t));
    }
    position = vec3(matModel * vec4(p, 1.0));
    normal = normalize(vec3(matNormal * vec4(vertexNormal, 0.0)));
    color = vertexColor;
    energy = vertexTexCoord.x;
    edgeDistance = vec3(vertexTexCoord2, 1.0 - vertexTexCoord2.x - vertexTexCoord2.y);
    edgeMask = int(vertexTexCoord.y + 0.5);
    gl_Position = mvp * vec4(p, 1.0);
}
