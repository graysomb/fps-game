#version 330
in vec3 vertexPosition;
in vec3 vertexNormal;
in vec2 vertexTexCoord;
in vec4 vertexTangent; // cube center.xyz, deterministic assembly order.w
in vec4 vertexColor;
uniform mat4 mvp;
uniform mat4 matNormal;
uniform float rebuild;
uniform int renderMode; // 0 solid, 1 emission, 2 ghost, 3 ghost emission
out vec3 normal;
out vec2 faceUV;
out vec4 color;
flat out float assembly;
void main(){
    assembly=clamp((rebuild-vertexTangent.w*.7)/.3,0.0,1.0);
    vec3 p=vertexPosition;
    if(renderMode<2 && assembly<1.0){
        vec3 offset=vec3(sin(vertexTangent.w*39.0)*.3,.12,cos(vertexTangent.w*27.0)*.3);
        p=vertexTangent.xyz+(p-vertexTangent.xyz)*(.3+.7*assembly)+offset*(1.0-assembly);
    }
    normal=normalize((matNormal*vec4(vertexNormal,0.0)).xyz);
    faceUV=vertexTexCoord;color=vertexColor;
    gl_Position=mvp*vec4(p,1.0);
}
