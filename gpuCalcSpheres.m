function [A, B] = gpuCalcSpheres(params)
%calculate all of the constants necessary to evaluate the scattered field
%estimate the order required to represent the scattered field for each sphere

%calculate the required order
numOrd = computeN_l(params.sphereRadius, params.material(1));

%set the refractive index for the sphere
n = params.material(2) + 1i*params.material(3);

%calculate the scattering coefficients
[A, B] = gpuScatteringCoefficients(params.material,params.sphereRadius);
        
