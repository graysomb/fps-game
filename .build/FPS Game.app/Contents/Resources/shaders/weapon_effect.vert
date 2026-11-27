#version 330
in vec3 vertexPosition;
in vec3 vertexNormal;
in vec4 vertexColor;
in vec2 vertexTexCoord;
in vec2 vertexTexCoord2;
uniform mat4 mvp;
uniform mat4 matModel;
uniform mat4 matNormal;
uniform int effectKind;
uniform float age;
uniform float effectTime;
uniform float seed;
out vec3 localPosition;
out vec3 normal;
out vec3 position;
out vec4 color;
out vec3 edge;
flat out int course;
flat out int edgeMask;
void main(){
 vec3 p=vertexPosition,n=vertexNormal;
 localPosition=p;
 if(effectKind==2){
  float phase=seed+p.z*1.0472;
  float t=p.x,angle=t*8.0+phase+effectTime*1.1;
  float r=.18+2.2*(1.0-t);
  p=vec3(cos(angle)*r,sin(t*9.0+phase)*r*.25,sin(angle)*r);
  p.y+=vertexPosition.y*.018;
  n=vec3(0,1,0);
 }
 position=vec3(matModel*vec4(p,1));normal=normalize(vec3(matNormal*vec4(n,0)));
 course=int(floor(vertexTexCoord.x/8.0));edgeMask=int(vertexTexCoord.y+.5);
 edge=vec3(vertexTexCoord2,1-vertexTexCoord2.x-vertexTexCoord2.y);
 color=vertexColor;gl_Position=mvp*vec4(p,1);
}
