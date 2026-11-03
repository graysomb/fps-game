#version 330
in vec2 fragTexCoord;
in vec4 fragColor;
uniform float effectTime;
out vec4 finalColor;
void main() {
    float tail=smoothstep(0.0,0.32,fragTexCoord.x);
    float edge=sin(3.14159265*fragTexCoord.y);
    float core=pow(max(0.0,edge),8.0);
    float ripple=0.85+0.15*sin(fragTexCoord.x*36.0-effectTime*28.0);
    vec3 blue=mix(vec3(0.025,0.20,1.0),vec3(0.4,0.85,1.0),core);
    finalColor=vec4(blue,fragColor.a*tail*edge*ripple*(1.0-smoothstep(0.94,1.0,fragTexCoord.x)));
}

