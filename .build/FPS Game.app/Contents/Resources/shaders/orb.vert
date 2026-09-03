#version 330

// Input vertex attributes
in vec3 vertexPosition;
in vec2 vertexTexCoord;
in vec3 vertexNormal;

// Input uniform values
uniform mat4 mvp;
uniform mat4 matModel;
uniform mat4 matView;

// Output vertex attributes (to fragment shader)
out vec3 fragPosition;
out vec2 fragTexCoord;
out vec3 fragNormal;
out vec3 viewPos;

void main()
{
    // Calculate fragment position in world space
    fragPosition = vec3(matModel * vec4(vertexPosition, 1.0));
    
    // Calculate fragment normal in world space (simplified)
    fragNormal = normalize(vec3(matModel * vec4(vertexNormal, 0.0)));
    
    // Pass texture coordinates
    fragTexCoord = vertexTexCoord;
    
    // Get view position from view matrix (inverse translation)
    // This assumes standard view matrix construction
    viewPos = vec3(inverse(matView)[3]);

    // Calculate final vertex position
    gl_Position = mvp * vec4(vertexPosition, 1.0);
}
