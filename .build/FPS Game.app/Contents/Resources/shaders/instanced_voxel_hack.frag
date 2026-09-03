#version 330

in vec3 fragPosition;
in vec2 fragTexCoord;
in vec4 fragColor;
in vec3 fragNormal;

out vec4 finalColor;

void main()
{
    float fu = fract(fragTexCoord.x);
    float fv = fract(fragTexCoord.y);
    bool isEdge = (fu < 0.02 || fu > 0.98 || fv < 0.02 || fv > 0.98);
    
    vec3 result = fragColor.rgb;
    if (isEdge) result *= 0.1;
    finalColor = vec4(result, fragColor.a);
}
