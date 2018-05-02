function E_t = cpuSumEf(params, E_s, E_i, E_f)
%add the incident field to the sum of scattered fields

E_t = zeros(params.simRes,params.simRes);    %total field

%total field
E_t(params.psMask<params.sphereRadius) = E_i(params.psMask<params.sphereRadius);
E_t(params.psMask>=params.sphereRadius) = E_f(params.psMask>=params.sphereRadius) + E_s(params.psMask>=params.sphereRadius);
E_t(isnan(E_t))=0; E_t(isinf(E_t))=0;

