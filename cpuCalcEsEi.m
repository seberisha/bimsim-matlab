function [E_s, E_i] = cpuCalcEsEi(params)

%scattered and internal fields
[E_s, E_i] = cpuComputeScatteredFields(params, params.material);
