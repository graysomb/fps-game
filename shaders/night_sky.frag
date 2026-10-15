#version 330
uniform vec3 viewForward, viewRight, viewUp;
uniform vec4 viewport; // width, height, tan(fov/2), aspect
uniform vec3 hazeColor;
uniform sampler2D spaceTexture;
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
    // Keep the original fog bank and silhouettes, revealing space high overhead.
    vec3 color=mix(hazeColor,space,smoothstep(.35,.94,altitude));
    finalColor=vec4(color,1.);
}
