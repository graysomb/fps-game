#version 330
uniform vec3 viewForward, viewRight, viewUp;
uniform vec4 viewport; // width, height, tan(fov/2), aspect
uniform vec3 hazeColor;
uniform sampler2D spaceTexture;
uniform vec4 meteorTrail; // head.xy, tail.zw in the celestial projection
uniform float meteorGlow;
out vec4 finalColor;
void main() {
    vec2 screen=(gl_FragCoord.xy/viewport.xy)*2.-1.;
    vec3 ray=normalize(viewForward+screen.x*viewport.z*viewport.w*viewRight+screen.y*viewport.z*viewUp);
    float altitude=max(ray.y,0.);
    // Stereographic upper hemisphere: no pole pinching or longitude seam.
    // Preserve the original ceiling's richly patterned texture at a distant scale.
    // Texture repetition belongs to the celestial sphere, never to world tiles.
    vec2 uv=.5+2.*ray.xz/(1.+altitude);
    vec3 space=texture(spaceTexture,uv).rgb;
    if(meteorGlow>0.) {
        vec2 sky=ray.xz/(1.+altitude),line=meteorTrail.xy-meteorTrail.zw;
        float along=clamp(dot(sky-meteorTrail.zw,line)/max(dot(line,line),.000001),0.,1.);
        float distance=length(sky-mix(meteorTrail.zw,meteorTrail.xy,along));
        float width=mix(.0006,.0045,along),aa=max(length(fwidth(sky))*.5,.0001);
        float streak=(1.-smoothstep(width,width+aa,distance))*along;
        float halo=exp(-distance*distance/.00015)*along*.6;
        float head=exp(-dot(sky-meteorTrail.xy,sky-meteorTrail.xy)/.000045)*2.2;
        space+=vec3(.75,.87,1.)*(streak*1.5+halo+head)*meteorGlow;
    }
    // Keep the original fog bank and silhouettes, revealing space high overhead.
    vec3 color=mix(hazeColor,space,smoothstep(.35,.94,altitude));
    finalColor=vec4(color,1.);
}
