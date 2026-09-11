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
out vec3 fragNormal;
out vec4 fragColor;

void main()
{
    // Unpack color from 4th row (indices [0][3], [1][3], [2][3])
    float r = instanceTransform[0][3];
    float g = instanceTransform[1][3];
    float b = instanceTransform[2][3];
    fragColor = vec4(r, g, b, 1.0);

    // Clean transform matrix for geometry
    mat4 cleanTransform = instanceTransform;
    cleanTransform[0][3] = 0.0;
    cleanTransform[1][3] = 0.0;
    cleanTransform[2][3] = 0.0;
    cleanTransform[3][3] = 1.0;

    // World position and normal
    vec4 worldPos = cleanTransform * vec4(vertexPosition, 1.0);
    fragPosition = worldPos.xyz;
    fragNormal = normalize(mat3(cleanTransform) * vertexNormal);

    // Clip space position
    gl_Position = matProjection * matView * worldPos;
}
