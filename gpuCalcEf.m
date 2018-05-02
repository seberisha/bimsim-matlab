function E_f = gpuCalcEf(params)

%focused field
E_f = gpuComputeEf(params);

E_f(isnan(E_f))=0; E_f(isinf(E_f))=0;   %do this if computing fft later