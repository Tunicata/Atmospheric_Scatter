#version 450 core

#define M_PI 3.1415926535897932384626433832795

layout(location = 0) in vec3 fsPosition;

layout(location = 0) out vec4 outColor;

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

vec2 RaySphereIntersection(vec3 rayOrigin, vec3 rayDir, float sphereRadius)
{
    // Assume sphere is centered at the origin

    float a = dot(rayDir, rayDir);
    float b = 2.0f * dot(rayOrigin, rayDir);
    float c = dot(rayOrigin, rayOrigin) - (sphereRadius * sphereRadius);
    float d = b * b - 4.0f * a * c;

    if (d < 0.0f) {
        return vec2(1e5, -1e5);
    }

    float sqrtD = sqrt(d);
    return (vec2(-b - sqrtD, -b + sqrtD) / (2.0f * a));
}

vec3 IntegrateScattering(vec3 rayDir, vec3 rayStart)
{
    vec3 sunDir = normalize(ubo.sunPos);
    vec2 intersection = RaySphereIntersection(rayStart, rayDir, ubo.atmosphereRadius);
    
    // If Intersects behind
    if (intersection.x > intersection.y)
        return vec3(0.0, 0.0, 0.0);

    // Rayleigh and Mie Phase functions
    float mu = dot(rayDir, sunDir);
    float mu_2 = mu * mu;
    
    float phaseR = 3.0 / (16.0 * M_PI) * (1.0 + mu_2);

    float g_2 = ubo.g * ubo.g;
    float phaseM = 3.0 / (8.0 * M_PI) *((1.0 - g_2) * (1.0 + mu_2)) /((2.0 + g_2) * pow(1.0 + g_2 - 2.0 * ubo.g * mu, 1.5));

    // compute density along the ray
    float viewRayStart = 0.0f;

    // Segment size
    intersection.y = min(intersection.y, RaySphereIntersection(rayStart, rayDir, ubo.planetRadius).x);
    float viewStepSize = (intersection.y - intersection.x) / float(ubo.viewSamples);
    
    vec3 scatterR = vec3(0);
    vec3 scatterM = vec3(0);

    // Optical depth 
    float viewOpticalDepthR = 0.0;
    float viewOpticalDepthM = 0.0;

    // Ray Marching Part I: View ray
    for (int i = 0; i < ubo.viewSamples; ++i)
    {
        vec3 viewSample = rayStart + rayDir * (viewRayStart + viewStepSize * 0.5);
        float heightView = length(viewSample) - ubo.planetRadius;

        float viewDensityR = exp(-heightView / ubo.hDensityRayleigh) * viewStepSize;
        float viewDensityM = exp(-heightView / ubo.hDensityMie) * viewStepSize;
        
        viewOpticalDepthR += viewDensityR;
        viewOpticalDepthM += viewDensityM;

        // Ray Marching Part II: Sun ray
        float sunRayStep = RaySphereIntersection(viewSample, sunDir, ubo.atmosphereRadius).y / float(ubo.lightSamples);
        float sunRayStart = 0.0;

        // Optical depth 
        float sunOpticalDepthLightR = 0.0;
        float sunOpticalDepthLightM = 0.0;
        
        for (int j = 0; j < ubo.lightSamples; ++j)
        {
            vec3 sunSample = viewSample + sunDir * (sunRayStart + sunRayStep * 0.5);
            float heightSun = length(sunSample) - ubo.planetRadius;

            sunOpticalDepthLightR += exp(-heightSun / ubo.hDensityRayleigh) * sunRayStep;
            sunOpticalDepthLightM += exp(-heightSun / ubo.hDensityMie) * sunRayStep;

            // Next sun ray sample
            sunRayStart += sunRayStep;
        }

        //        vec3 sigma_R = ubo.scatterRayleigh * viewDensityR;
        //        float mieDensity = viewDensityM;
        //        float sigma_mieS = ubo.scatterMie  * mieDensity;
        //        float sigma_mieT = (ubo.scatterMie  + ubo.Absorb_M) * mieDensity;
        //        vec3 ozone = vec3(0.65e-3f, 1.881e-3f, 0.085e-3f) * max(0.0f, 1 - 0.5 * abs(height - 25) / 15);
        //
        //        vec3 sigmaS = sigma_R + sigma_mieS;
        //        vec3 sigmaT = sigma_R + sigma_mieT + ozone;
        //        
        //        vec3 eyeTrans = exp(-sumSigmaT - 0.5 * sigmaT);
        
        vec3 Tr = ubo.scatterRayleigh * (viewOpticalDepthR + sunOpticalDepthLightR);
        //  Local Mie extinction coeff. = 1.1
        float Tm = ubo.scatterMie * 1.1f * (viewOpticalDepthM + sunOpticalDepthLightM);
        
        vec3 extinction = exp(-(Tr + Tm));
        
        // Accumulate the scattering 
        scatterR += viewDensityR * extinction;
        scatterM += viewDensityM * extinction;

        // Next view ray sample
        viewRayStart += viewStepSize;
    }

    return ubo.sunIntensity * (scatterR  * ubo.scatterRayleigh * phaseR + scatterM * ubo.scatterMie * phaseM);
}

void main()
{
    vec3 aColor = IntegrateScattering(normalize(fsPosition - ubo.viewPos), ubo.viewPos);
    
    aColor = mix(aColor, (1.0 - exp(-1.0 * aColor)), ubo.toneMappingFactor);

    outColor = vec4(aColor, 1.0);
}
