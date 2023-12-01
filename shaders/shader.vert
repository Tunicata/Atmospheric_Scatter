#version 450 core

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec2 texCoord;

layout(location = 0) out vec3 fsPosition;

layout(binding = 0) uniform UniformBufferObject {
    mat4 AtmosModel;        // Atmosphere Model matrix
    mat4 MVP;               // Model View Projection matrix

    vec3 viewPos;           // Position of the viewer
    vec3 sunPos;            // Position of the sun, light direction

    int viewSamples;        // view Ray (view -> sample point) Sample
    int lightSamples;       // light Ray (sample point) Sample

    vec3 sunIntensity;     // Intensity of the sun
    float planetRadius;     // Radius of the planet [m]
    float atmosphereRadius; // Radius of the atmosphere [m]
    vec3  scatterRayleigh;  // Rayleigh scattering coefficient
    float scatterMie;       // Mie scattering coefficient

    float hDensityRayleigh; // Rayleigh scale height
    float hDensityMie;      // Mie scale height
    float g;                // Mie scattering direction - Mie Asymmetry Coefficient

    float toneMappingFactor;
} ubo;

void main()
{
    vec4 posVec4 = vec4(position, 1.0);
    fsPosition = vec3(ubo.AtmosModel * posVec4);
	gl_Position = ubo.MVP * posVec4;
}
