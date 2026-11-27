#version 330
in vec2 fragTexCoord;
uniform sampler2D texture0;
uniform sampler2D worldDepth;
uniform vec4 settings; // aspect, tan(fov/2), near, far
uniform vec3 fogColor;
out vec4 finalColor;
void main() {
    vec4 base=texture(texture0,fragTexCoord);
    float depth=texture(worldDepth,fragTexCoord).r;
    if(depth>=0.999999){finalColor=base;return;}
    float z=depth*2.0-1.0;
    float linearDepth=2.0*settings.z*settings.w/(settings.w+settings.z-z*(settings.w-settings.z));
    vec2 xy=(fragTexCoord*2.0-1.0)*vec2(settings.x,1.0)*settings.y;
    float distance=linearDepth*sqrt(1.0+dot(xy,xy));
    float fog=smoothstep(24.0,56.0,distance);
    finalColor=vec4(mix(base.rgb,fogColor,fog),base.a);
}
