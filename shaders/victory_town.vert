#version 330
in vec3 vertexPosition;
in vec3 vertexNormal;
in vec4 vertexColor;
uniform mat4 mvp;
uniform float townTime;
uniform int townMode;
out vec3 p;
out vec3 n;
out vec4 tint;
void main(){
 p=vertexPosition;n=vertexNormal;tint=vertexColor;
 if(townMode==1 || townMode==6){float w=sin(townTime*1.570796+p.x*1.3+p.z)*.045;p.x+=w;p.z+=w*.4;}
 if(townMode==3)p.y+=.012*sin(p.x*3.+townTime*3.14159)*cos(p.z*2.-townTime*1.570796);
 gl_Position=mvp*vec4(p,1.);
}
