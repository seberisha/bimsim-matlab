function [A, B] = cpuCalcSpheres(params)
%calculate all of the constants necessary to evaluate the scattered field
%estimate the order required to represent the scattered field for each sphere

%calculate the scattering coefficients
[A, B] = cpuScatteringCoefficients(params.material,params.sphereRadius);
        
