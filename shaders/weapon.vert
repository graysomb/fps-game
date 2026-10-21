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
uniform float charge;
uniform float cooldown;
uniform float deployment;
uniform float attached;
uniform float shotAge;
flat out int partId;
out float visibility;
mat3 turnX(float a){float c=cos(a),s=sin(a);return mat3(1,0,0,0,c,s,0,-s,c);}
mat3 turnY(float a){float c=cos(a),s=sin(a);return mat3(c,0,-s,0,1,0,s,0,c);}
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
    vec3 n=vertexNormal;
    partId=int(floor(vertexTexCoord.x/8.0));
    energy=mod(vertexTexCoord.x,8.0);
    visibility=1.0;
    float kick=max(0.0,1.0-shotAge/.18);
    if(partId==1)p.z+=kick*kick*(weaponKind==2?.11:.055);
    if(partId==2 || partId==3){
        float side=partId==2?1.0:-1.0;
        vec3 pivot=vec3(0,side*.23,-.49);
        mat3 r=turnX(-side*.075*activity);
        p=pivot+r*(p-pivot);n=r*n;
    }
    if(partId==4){ // Rocket becomes visible as its cooldown completes.
        visibility=1.0-smoothstep(.15,.85,cooldown);
        p=vec3(0,0,-.72)+(p-vec3(0,0,-.72))*visibility;
    }
    if(partId>=5 && partId<=8){
        float side=(partId==5 || partId==7)?-1.0:1.0;
        if(partId<=6)p.x+=side*.06*attached;else p.z+=side*.06*attached;
    }
    if(partId==9 || partId==10){
        float side=partId==9?-1.0:1.0;
        vec3 pivot=vec3(side*.32,0,-.31);
        mat3 r=turnY(side*.85*(1.0-deployment));
        p=pivot+r*(p-pivot);n=r*n;
    }
    if(partId==11 || partId==12){
        visibility=max(0.0,1.0-shotAge/.10);
        p.y+=sin(effectTime*91.0+p.z*13.0)*.025*visibility;
    }
    if (energy > 2.5) {
        float t = vertexPosition.z;
        float radius = (launcher != 0 ? 0.68 : 0.54) + 0.025 * sin(t * 31.0 - effectTime * 8.0);
        float angle = vertexPosition.x + effectTime * 1.8 +
                      0.06 * sin(t * 43.0 - effectTime * 11.0) + vertexPosition.y * 0.040 / radius;
        p = vec3(cos(angle) * radius, sin(angle) * radius * 0.85 + (launcher != 0 ? 0.14 : 0.0),
                 mix(launcher != 0 ? 0.94 : 0.70, launcher != 0 ? -1.60 : -1.02, t));
    }
    position = vec3(matModel * vec4(p, 1.0));
    normal = normalize(vec3(matNormal * vec4(n, 0.0)));
    color = vertexColor;
    edgeDistance = vec3(vertexTexCoord2, 1.0 - vertexTexCoord2.x - vertexTexCoord2.y);
    edgeMask = int(vertexTexCoord.y + 0.5);
    gl_Position = mvp * vec4(p, 1.0);
}
