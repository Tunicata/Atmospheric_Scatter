#pragma once
#include <glm/vec3.hpp>

typedef struct Sun
{
    glm::vec3 sunDir;
    glm::vec3 sunIntensity;
    float sunAngle;
    
}Sun;

typedef struct Atmosphere
{
    alignas(4) float toneMappingFactor;    ///< Whether tone mapping is applied
    
    glm::vec3 scatterRayleigh = glm::vec3(5.8e-3f, 13.5e-3f, 33.1e-3f);
    float scatterMie = 21e-3f;
    
    float  hDensityRayleigh = 8.f;
    float hDensityMie  = 1.2f;

    float asymmetryMie = 0.8f;
    
    float  planetRadius = 6360;
    float  atmosphereRadius = 6460;

    glm::vec3 absorbOzone = glm::vec3(0.65f, 1.881f, 0.085f);
    float  ozoneCenterHeight = 25;

    float  ozoneThickness = 30;
}Atmosphere;
