function E_f = cpuCalcEf(params)

%focused field
E_f = cpuComputeEf(params);

E_f(isnan(E_f))=0; E_f(isinf(E_f))=0;   %do this if computing fft later