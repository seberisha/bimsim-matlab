function [E_t, E_f] = gpuSimulateScattering(params)

%compute a set of plane waves for Monte-Carlo simulation
params.k_j = gpuArray(monteCarlo(params.s,params.samples, params.kVec, params.condenserNA_in, params.condenserNA_out));

% simulate the focused/incident near field
E_f = gpuCalcEf(params);

[params.A, params.B] = gpuCalcSpheres(params);

[E_s, E_i] = gpuCalcEsEi(params);

E_t = gpuSumEf(params, E_s, E_i, E_f);
