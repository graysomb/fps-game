#version 330
in vec3 p;
in vec3 n;
in vec4 tint;
uniform vec3 eye;
uniform int townMode;
uniform int emissionPass;
uniform float townTime;
uniform float townPower;
uniform mat4 lightMatrix;
uniform sampler2D shadowMap;
out vec4 finalColor;
float shadow(){vec4 clip=lightMatrix*vec4(p,1.);vec3 uv=clip.xyz/clip.w*.5+.5;if(any(lessThan(uv,vec3(0.)))||any(greaterThan(uv,vec3(1.))))return 1.;float sum=0.;for(int x=-1;x<=1;x++)for(int y=-1;y<=1;y++){float d=dot(texture(shadowMap,uv.xy+vec2(x,y)/2048.).rgb,vec3(1.,1./255.,1./65025.));sum+=uv.z-.0012<d?1.:.22;}return sum/9.;}
float hash(vec3 q){return fract(sin(dot(q,vec3(127.1,311.7,74.7)))*43758.5453);}
void main(){
 float power=(p.x < -17. && p.x > -24. && p.z < -27. && p.z > -30.)?townPower:1.;
 vec3 normal=normalize(n),base=tint.rgb;float sun=max(dot(normal,normalize(vec3(-.45,.8,-.35))),0.);
 if(townMode==0){float wear=hash(floor(p*13.));base*=.94+wear*.12;if(tint.a<.999){vec3 grid=abs(fract(p*2.+.5)-.5)*(1.-abs(normal));base*=mix(1.,.91,smoothstep(.47,.495,max(grid.x,max(grid.y,grid.z))));}}
 if(townMode==4){float grain=sin(p.x*22.+sin(p.z*9.)*2.+sin(p.y*11.));base*=.96+.04*grain;}
 if(townMode==6){float weave=sin(p.x*160.)*sin(p.y*160.);base*=.96+.04*weave;}
 if(townMode==8)base*=.75+.25*hash(floor(p*35.));
 vec3 c=base*(vec3(.38,.43,.50)+vec3(.56,.57,.53)*sun*shadow());
 vec3 lamps[7]=vec3[7](vec3(0.,2.6,0.),vec3(21.5,2.6,14.),vec3(-21.5,2.6,14.),vec3(12.,1.5,12.),vec3(21.,3.,-27.),vec3(-21.,3.,-27.),vec3(15.,2.,8.));
 for(int i=0;i<7;i++){vec3 delta=lamps[i]-p;float incidence=.3+.7*max(dot(normal,normalize(delta)),0.);c+=base*vec3(1.1,.65,.30)*incidence*2.1*(i==5?townPower:1.)/(1.+dot(delta,delta)*.55);}
 if(townMode==5){vec3 h=normalize(normalize(eye-p)+normalize(vec3(-.4,.8,-.3)));c+=vec3(.22,.26,.28)*pow(max(dot(normal,h),0.),32.);}
 if(townMode==2)c=base*mix(.12,1.,power)*(1.08+.05*sin(townTime*1.570796+p.x));
 if(townMode==3){float ripple=sin(length(p.xz-vec2(20.2,9.3))*18.-townTime*6.283185);c+=vec3(.12,.22,.23)*pow(max(0.,ripple),10.);c+=vec3(.08,.13,.17)*pow(1.-abs(dot(normal,normalize(eye-p))),3.);}
 if(townMode==7)c=mix(vec3(.13,.23,.29),vec3(.7,.85,.9),pow(1.-abs(dot(normal,normalize(eye-p))),4.));
 float fog=1.-exp(-length(p-eye)*.006);c=mix(c,vec3(.075,.12,.17),fog);
 if(emissionPass==1)c=townMode==2?base*.28*power:vec3(0.);
 finalColor=vec4(c,townMode==7?tint.a:1.);
}
