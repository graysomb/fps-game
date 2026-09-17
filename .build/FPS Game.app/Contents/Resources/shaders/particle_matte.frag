#version 330

in vec3 fragPosition;
in vec3 fragNormal;
in vec4 fragColor;

out vec4 finalColor;

void main()
{
    vec3 N = normalize(fragNormal);
    vec3 L = normalize(vec3(0.35, 0.85, 0.40));

    // Soft wrapped diffuse for smooth matte shading (no harsh specular)
    float NdotL = dot(N, L);
    float wrappedDiff = clamp((NdotL + 0.25) / 1.25, 0.0, 1.0);

    // Two-tone hemisphere ambient (sky fill from above, soft ground bounce from below)
    float hemi = N.y * 0.5 + 0.5;
    vec3 ambient = mix(vec3(0.20, 0.20, 0.24), vec3(0.42, 0.40, 0.38), hemi);

    // Pure matte diffuse lighting
    vec3 lighting = ambient + wrappedDiff * vec3(0.70, 0.68, 0.65);
    vec3 result = fragColor.rgb * lighting;

    finalColor = vec4(result, fragColor.a);
}
