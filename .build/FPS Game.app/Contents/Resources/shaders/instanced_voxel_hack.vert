#version 330

// Input vertex attributes
in vec3 vertexPosition;
in vec2 vertexTexCoord;
in vec3 vertexNormal;

// Input instance attributes (matrix columns)
in mat4 instanceTransform;

// Input uniform values
uniform mat4 mvp;
uniform mat4 matModel;
uniform mat4 matView;
uniform mat4 matProjection;

// Output vertex attributes (to fragment shader)
out vec3 fragPosition;
out vec2 fragTexCoord;
out vec4 fragColor;
out vec3 fragNormal;

void main()
{
    // Unpack color from the bottom row components (which are usually 0 in affine transforms)
    // Col 0, Row 3 -> m3
    // Col 1, Row 3 -> m7
    // Col 2, Row 3 -> m11
    float r = instanceTransform[0][3];
    float g = instanceTransform[1][3];
    float b = instanceTransform[2][3];
    
    vec4 color = vec4(r, g, b, 1.0);

    // Reconstruct pure transform matrix for geometry
    // We take the top-left 3x3 for rotation/scale/shear
    // And the translation column
    mat4 cleanTransform = instanceTransform;
    
    // Clear the color slots to 0.0 for the geometry transform
    cleanTransform[0][3] = 0.0;
    cleanTransform[1][3] = 0.0;
    cleanTransform[2][3] = 0.0;
    cleanTransform[3][3] = 1.0; // Ensure w is 1

    // Calculate vertex position
    vec4 worldPos = cleanTransform * vec4(vertexPosition, 1.0);
    fragPosition = worldPos.xyz;
    fragTexCoord = vertexTexCoord;
    fragColor = color;
    
    // Normal transform (using the 3x3 submatrix)
    fragNormal = normalize((cleanTransform * vec4(vertexNormal, 0.0)).xyz);

    // Calculate final vertex position
    gl_Position = matProjection * matView * worldPos;
}
