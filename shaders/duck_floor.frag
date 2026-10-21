#version 330
in vec2 uv;
uniform float duckTime;
uniform float swimSpeed;
uniform int showWake;
out vec4 finalColor;
void main() {
    vec2 p=(uv-.5)*vec2(2.4,3.5);
    // Mesh center sits slightly behind the duck. Shadow stays under its belly.
    vec2 contact=(p+vec2(0,.38))/vec2(.35,.59);
    float shadow=exp(-dot(contact,contact)*2.4)*.32;
    float wake=0.0;
    if(showWake!=0) {
        vec2 q=p+vec2(0,.32);
        float radius=length(q/vec2(1.0,1.35));
        float phase=radius*13.0-duckTime*3.14159265;
        float ring=1.0-smoothstep(.05,.14,abs(sin(phase)));
        float rear=smoothstep(-.35,.35,q.y);
        wake=ring*smoothstep(.35,.65,radius)*(1.0-smoothstep(.65,1.45,radius))*rear*.16*swimSpeed;
    }
    float alpha=shadow+wake;
    if(alpha<.004)discard;
    finalColor=vec4(mix(vec3(.07,.11,.12),vec3(.55,.71,.70),wake/max(alpha,.001)),alpha);
}
