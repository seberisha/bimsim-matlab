function [E_s , E_i] = gpuExtendedSourceScatteredFields(params, pEs, pEi)


%the vector from the focal point to the center of the sphere
c = params.gpu_ps - params.gpu_pf;

for i=1:params.samples
    
    %multiply by the phase
    phase = exp(1i.*params.wavNum.*params.k_j(:,i)'*c');
    params.gpu_Es = params.gpu_Es + phase.*pEs(:,:,i);
    params.gpu_E_i = params.gpu_E_i + phase.*pEi(:,:,i);

end

params.gpu_Es(params.psMask<params.a) = 0;%E_i(r<a);
E_s = params.gpu_Es;

params.gpu_E_i(params.psMask>params.a) = 0;
E_i = params.gpu_E_i;
