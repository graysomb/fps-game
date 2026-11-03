#version 330
in vec3 normal;
in vec3 localPosition;
in vec4 color;
uniform int duckMaterial; // 0 matte plumage, 1 iridescent feathers, 2 bill/eyes
out vec4 finalColor;
void main() {
    vec3 n=normalize(normal),light=normalize(vec3(-0.45,0.85,-0.55));
    float diffuse=max(dot(n,light),0.0);
    vec3 base=color.rgb;
    if(duckMaterial==1) {
        // Restrained green/teal sheen that follows the facets as the duck turns.
        float sheen=pow(max(dot(n,normalize(vec3(.65,.45,.2))),0.0),3.0);
        base=mix(base,base*vec3(.7,1.15,1.6),sheen*.35);
    }
    vec3 lit=base*(.55+.45*diffuse);
    float shine=pow(max(dot(n,normalize(vec3(-.3,.9,-.6))),0.0),duckMaterial==2?40.0:18.0);
    lit+=vec3(.10,.11,.10)*shine*(duckMaterial==0?.12:.55);
    finalColor=vec4(lit,color.a);
}
