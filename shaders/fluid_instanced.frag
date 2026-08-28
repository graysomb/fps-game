#version 330

in vec3 fragPosition;
in vec3 fragNormal;

uniform vec3 viewPos;
uniform int uSurfaceMode;

out vec4 finalColor;

void main()
{
    vec3 normal = normalize(fragNormal);
    vec3 viewDirection = normalize(viewPos - fragPosition);
    float fresnel = pow(clamp(1.0 - dot(normal, viewDirection), 0.0, 1.0), 2.2);
    float light = 0.35 + 0.65 * max(dot(normal, normalize(vec3(-0.35, 0.85, 0.4))), 0.0);
    vec3 deepColor = vec3(0.035, 0.24, 0.52);
    vec3 rimColor = vec3(0.42, 0.82, 1.0);
    vec3 color = mix(deepColor * light, rimColor, fresnel * 0.75);
    float alpha = (uSurfaceMode != 0) ? (0.82 + 0.16 * fresnel) : 1.0;
    finalColor = vec4(color, alpha);
}
