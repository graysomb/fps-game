#version 330
in vec3 position;
in vec3 normal;
in vec3 localPosition;
in vec4 color;
in float visibility;
flat in int partId;
uniform vec3 eyePosition;
uniform int emissionOnly;
uniform float viewHide;
out vec4 finalColor;
float hash(vec3 p){return fract(sin(dot(p,vec3(127.1,311.7,74.7)))*43758.5453);}
void main(){
 if(visibility<.005)discard;
 if(partId==11){finalColor=vec4(vec3(1,.65,.24)*(1.0+visibility),visibility);return;}
 if(emissionOnly!=0){finalColor=vec4(0,0,0,1);return;}
 bool polymer=partId==31 || color.r<.17;
 vec3 n=normalize(normal),l=normalize(vec3(-.45,.8,-.45)),v=normalize(eyePosition-position);
 float grain=hash(floor(localPosition*650.0));
 // Derivative attenuation avoids sparkling when the textured grip is distant.
 float detail=1.0-smoothstep(.002,.012,length(fwidth(localPosition)));
 float roughness=polymer?.78:.27;
 float spec=pow(max(0.0,dot(n,normalize(v+l))),mix(90.0,8.0,roughness));
 vec3 base=color.rgb*(.78+detail*(grain-.5)*.12);
 float fresnel=.08+.32*pow(1.0-max(0.0,dot(n,v)),5.0);
 base=base*(.65+.55*max(0.0,dot(n,l)))+vec3(.56,.62,.70)*spec*(polymer?.045:.32)+fresnel*.06;
 finalColor=vec4(base*(1.0-viewHide),1.0-viewHide);
}
