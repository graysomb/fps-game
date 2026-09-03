#version 330

// Input vertex attributes (from vertex shader)
in vec3 fragPosition;
in vec2 fragTexCoord;
in vec3 fragNormal;
in vec3 viewPos;

// Input uniform values
uniform vec4 colDiffuse;

// Output fragment color
out vec4 finalColor;

void main()
{
    // Base orb color (glowing blue)
    vec3 baseColor = vec3(0.2, 0.6, 1.0);
    
    // View direction
    vec3 viewDir = normalize(viewPos - fragPosition);
    
    // Fresnel effect (rim lighting)
    float fresnel = dot(fragNormal, viewDir);
    fresnel = clamp(1.0 - fresnel, 0.0, 1.0);
    fresnel = pow(fresnel, 2.0); // Sharpen the rim
    
    // Center glow (inverse fresnel)
    float centerGlow = 1.0 - fresnel;
    
    // Combine
    vec3 emission = baseColor * 2.0; // Bright emission
    vec3 color = mix(baseColor, vec3(0.8, 0.9, 1.0), fresnel); // White-ish rim
    
    // Alpha for transparency (ghostly look)
    float alpha = 0.6 + 0.4 * fresnel;
    
    finalColor = vec4(color * emission, alpha);
}
