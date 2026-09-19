#version 330

in vec3 fragNormal;

out vec4 finalColor;

void main()
{
    vec3 normal = normalize(fragNormal);
    vec3 lightDirection = normalize(vec3(-0.35, 0.85, 0.4));
    float light = 0.38 + 0.62 * max(dot(normal, lightDirection), 0.0);
    vec3 water = vec3(0.035, 0.30, 0.68);
    finalColor = vec4(water * light, 0.96);
}
