#version 330
in vec3 localPosition;
in vec3 position;
in vec3 normal;
in vec4 color;
in vec3 edge;
flat in int course;
flat in int edgeMask;
uniform vec3 eyePosition;
uniform vec3 tint;
uniform int effectKind;
uniform int emissionOnly;
uniform float age;
uniform float lifetime;
uniform float progress;
uniform float seed;
out vec4 finalColor;
float hash(vec3 p){return fract(sin(dot(floor(p*55.0)+seed,vec3(127.1,311.7,74.7)))*43758.5453);}
void main(){
 vec3 n=normalize(normal),v=normalize(eyePosition-position),l=normalize(vec3(-.45,.8,-.45));
 float diffuse=.55+.45*max(0.0,dot(n,l));
 vec3 light=vec3(0),base=color.rgb*tint*diffuse;
 if(effectKind==0){
  float dissolve=smoothstep(lifetime-.65,lifetime,age);
  if(hash(localPosition)<dissolve)discard;
  light=tint*.36;
  base=tint*(.4+.6*diffuse)+pow(max(0.0,dot(n,normalize(v+l))),48.0)*vec3(.7);
 }
 if(effectKind==1||effectKind==2){
  light=tint*(.7+.3*sin(localPosition.x*8.0-age*9.0));base=light;
  if(effectKind==2){
   // A moving packet travels down a continuous inward spiral; no wrap seam.
   float packet=fract(localPosition.x*2.0-age*.9+seed*.13+localPosition.z*.16);
   if(packet>.62 || localPosition.x<.02||localPosition.x>.98)discard;
   light*=smoothstep(0.0,.12,packet)*(1.0-smoothstep(.42,.62,packet));base=light;
  }
 }
 if(effectKind==3){
  if(float(course)+.01>progress*7.0)discard;
  vec3 cov=smoothstep(vec3(0),max(fwidth(edge),vec3(.00001)),edge);
  float interior=1.0;
  if((edgeMask&1)!=0)interior=min(interior,cov.x);
  if((edgeMask&2)!=0)interior=min(interior,cov.y);
  if((edgeMask&4)!=0)interior=min(interior,cov.z);
  float fresh=clamp(1.0-(progress*7.0-float(course))*1.5,0,1)*(1.0-step(.999,progress));
  light=vec3(.1,.55,1)*(1.0-interior)*fresh*1.6;
  base=mix(vec3(.09,.12,.15),base,interior)+light;
 }
 finalColor=vec4(emissionOnly!=0?light:base,1);
}
