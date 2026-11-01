#version 330
out vec4 finalColor;
void main(){vec3 encoded=fract(gl_FragCoord.z*vec3(1.,255.,65025.));encoded-=encoded.yzz*vec3(1./255.,1./255.,0.);finalColor=vec4(encoded,1.);}

