#version 330

in vec3 fragNormal;
in vec3 fragPosition;

out vec4 finalColor;

void main()
{
    vec3 lightDir = normalize(vec3(0.35, 0.85, 0.25));
    float diffuse = 0.62 + 0.38 * max(dot(normalize(fragNormal), lightDir), 0.0);
    float heightTint = clamp(fragPosition.y / 40.0, 0.0, 1.0);
    vec3 water = mix(vec3(0.02, 0.30, 0.72), vec3(0.08, 0.55, 0.95), heightTint);
    finalColor = vec4(water * diffuse, 0.65);
}
