#version 330
in vec3 point;
in vec3 normal;
in vec2 uv;
in vec4 tint;
uniform vec4 colDiffuse;
uniform float campTime;
uniform int campMode;
uniform int emissionPass;
out vec4 finalColor;
float hash(vec2 p){return fract(sin(dot(p,vec2(127.1,311.7)))*43758.5453);}
float noise(vec2 p){vec2 i=floor(p),f=fract(p);f=f*f*(3.0-2.0*f);return mix(mix(hash(i),hash(i+vec2(1,0)),f.x),mix(hash(i+vec2(0,1)),hash(i+vec2(1,1)),f.x),f.y);}
void main(){
 vec4 c=tint*colDiffuse;
 if(campMode==2){
  float h=uv.y,flow=noise(vec2(uv.x*5.0,h*7.0-campTime*.75));
  float strand=exp(-pow((uv.x-.5+.12*sin(h*13.0-campTime))*(7.0+flow*6.0),2.0));
  float alpha=strand*smoothstep(0.0,.13,h)*(1.0-smoothstep(.3,1.0,h))*(.08+.14*flow);
  finalColor=vec4(c.rgb,alpha*c.a);return;
 }
 if(campMode==1){float h=clamp((point.y-.15)/.5,0.0,1.0);vec3 fire=mix(vec3(1.0,.84,.22),vec3(1.0,.24,.015),h);finalColor=vec4(fire*(1.0+.18*sin(campTime*9.0+point.x*17.0)),1);return;}
 float light=.48+.52*max(dot(normalize(normal),normalize(vec3(-.4,.8,.3))),0.0);
 float warmth=exp(-abs(point.y-.4)*3.0)*(.18+.035*sin(campTime*8.0)+.025*sin(campTime*13.0));
 vec3 color=c.rgb*light+vec3(1.0,.27,.025)*warmth;
 if(point.y>1.18 && point.y<1.21 && length(point.xz)<.55)color*=.97+.03*sin(length(point.xz)*35.0-campTime*2.0);
 finalColor=vec4(emissionPass!=0?vec3(0):color,c.a);
}
