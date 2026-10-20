#version 330
uniform vec2 resolution;
uniform float blend;
uniform int style;
uniform float railCharge;
out vec4 finalColor;
float stroke(float d,float width){return 1.0-smoothstep(width,width+1.3/resolution.y,abs(d));}
float segment(vec2 p,vec2 a,vec2 b){vec2 q=p-a,v=b-a;return length(q-v*clamp(dot(q,v)/dot(v,v),0.0,1.0));}
void main(){
 vec2 p=(gl_FragCoord.xy-resolution*.5)/resolution.y;
 float r=length(p),px=1.0/resolution.y,ink=0.0;
 vec3 color=vec3(.15,.88,1.0);

 if(style==4){
  color=vec3(.78,.35,1.0);

  ink=max(ink,stroke(p.x,.55*px)*step(.018,abs(p.y))*step(abs(p.y),.12));
  ink=max(ink,stroke(p.y,.55*px)*step(.018,abs(p.x))*step(abs(p.x),.12));
  for(int n=1;n<=5;n++){float v=float(n)*.02;ink=max(ink,stroke(abs(p.y)-v,px)*step(abs(p.x),.007));ink=max(ink,stroke(abs(p.x)-v,px)*step(abs(p.y),.007));}
  ink=max(ink,stroke(p.y+.085,2.0*px)*step(abs(p.x),.08)*step(p.x,-.08+.16*railCharge));
 }else if(style==1||style==5){
  color=style==1?vec3(1,.78,.16):vec3(1,.48,.15);
  vec2 q=abs(p);ink=max(stroke(q.x-.045,px)*step(q.y,.03),stroke(q.y-.03,px)*step(.025,q.x)*step(q.x,.045));
  if(style==5)for(int n=1;n<=3;n++)ink=max(ink,stroke(p.y+float(n)*.025,px)*step(abs(p.x),.018-float(n)*.003));
 }else if(style==2){
  color=vec3(1,.62,.1);ink=stroke(r-.065,px)*step(.026,min(abs(p.x),abs(p.y)));
 }else if(style==3){
  color=vec3(.3,1,.48);ink=max(stroke(p.x,.5*px)*step(.009,abs(p.y))*step(abs(p.y),.035),stroke(p.y,.5*px)*step(.009,abs(p.x))*step(abs(p.x),.035));
 }else{
  ink=max(stroke(segment(p,vec2(-.016,-.012),vec2(0,.004)),px),stroke(segment(p,vec2(.016,-.012),vec2(0,.004)),px));
 }
 ink=max(ink,1.0-smoothstep(1.0*px,2.0*px,r));
 finalColor=vec4(color,ink*.95*blend);
}
