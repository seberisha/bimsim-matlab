function [E_t, E_f] = cpuSimulateScattering(params)

%compute a set of plane waves for Monte-Carlo simulation
params.k_j = (monteCarlo(params.s,params.samples, params.kVec, params.condenserNA_in, params.condenserNA_out));

% simulate the focused/incident near field
E_f = cpuCalcEf(params);

[params.A, params.B] = cpuCalcSpheres(params);

[E_s, E_i] = cpuCalcEsEi(params);

E_t = cpuSumEf(params, E_s, E_i, E_f);
