#version 330
in vec3 vertexPosition;
in vec3 vertexNormal;
in vec2 vertexTexCoord;
in vec4 vertexColor;
uniform mat4 mvp;
uniform float campTime;
uniform int campMode;
out vec3 point;
out vec3 normal;
out vec2 uv;
out vec4 tint;
void main(){
 vec3 p=vertexPosition;point=p;normal=vertexNormal;uv=vertexTexCoord;tint=vertexColor;
 if(campMode==1){float h=max(0.0,p.y-.15);p.x+=sin(campTime*6.0+p.z*19.0)*h*.15;p.z+=cos(campTime*4.7+p.x*17.0)*h*.13;p.y+=h*.22*sin(campTime*8.0+p.x*23.0+p.z*15.0);}
 if(campMode==2){p.x+=p.y*(.65*sin(p.y*5.0-campTime*.65)+.3*sin(p.y*10.0-campTime*1.2));}
 gl_Position=mvp*vec4(p,1.0);
}
