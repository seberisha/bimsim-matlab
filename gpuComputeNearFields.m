function [E_f, E_t] = gpuComputeNearFields(params, materialIdx)

%focused field
E_f = gpuComputeEf(params, params.kVecs(materialIdx,:));

E_f(isnan(E_f))=0; E_f(isinf(E_f))=0;   %do this if computing fft later

%scattered and internal fields
[E_s, E_i] = gpuComputeScatteredFields(params, params.material(materialIdx,:));

E_t = zeros(params.simRes,params.simRes,'gpuArray');    %total field

%total field
E_t(params.psMask<params.a) = E_i(params.psMask<params.a);
E_t(params.psMask>=params.a) = E_f(params.psMask>=params.a) + E_s(params.psMask>=params.a);
E_t(isnan(E_t))=0; E_t(isinf(E_t))=0;
