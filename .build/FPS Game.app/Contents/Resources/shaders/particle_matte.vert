#version 330

// Input vertex attributes
in vec3 vertexPosition;
in vec2 vertexTexCoord;
in vec3 vertexNormal;

// Input instance attribute: xyz = position, w = radius
in vec4 instancePosRadius;

// Input uniform values
uniform mat4 mvp;
uniform mat4 matModel;
uniform mat4 matView;
uniform mat4 matProjection;
uniform vec4 particleColor;

// Output vertex attributes (to fragment shader)
out vec3 fragPosition;
out vec3 fragNormal;
out vec4 fragColor;

void main()
{
    fragColor = (particleColor.a > 0.0) ? particleColor : vec4(1.0, 0.627, 0.157, 1.0);

    // World position and normal (uniform sphere scale preserves normals)
    vec3 worldPos = vertexPosition * instancePosRadius.w + instancePosRadius.xyz;
    fragPosition = worldPos;
    fragNormal = vertexNormal;

    // Clip space position
    gl_Position = matProjection * matView * vec4(worldPos, 1.0);
}
