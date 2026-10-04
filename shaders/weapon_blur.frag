#version 330
in vec2 fragTexCoord;
uniform sampler2D texture0;
uniform vec2 direction;
out vec4 finalColor;
void main() {
    vec3 result = texture(texture0, fragTexCoord).rgb * 0.227027;
    result += texture(texture0, fragTexCoord + direction * 1.384615).rgb * 0.316216;
    result += texture(texture0, fragTexCoord - direction * 1.384615).rgb * 0.316216;
    result += texture(texture0, fragTexCoord + direction * 3.230769).rgb * 0.070270;
    result += texture(texture0, fragTexCoord - direction * 3.230769).rgb * 0.070270;
    finalColor = vec4(result, 1.0);
}
