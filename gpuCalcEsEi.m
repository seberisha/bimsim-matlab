function [E_s, E_i] = gpuCalcEsEi(params)

%scattered and internal fields
[E_s, E_i] = gpuComputeScatteredFields(params, params.material);
