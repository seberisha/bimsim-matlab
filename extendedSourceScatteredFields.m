function [E_s , E_i] = extendedSourceScatteredFields(params,pEs, pEi)

E_s = zeros(params.simRes,params.simRes);
E_i = E_s;

%the vector from the focal point to the center of the sphere
c = params.ps - params.pf;
for i=1:params.samples
    %multiply by the scattering coefficients B
    %sum all orders
    phase = exp(1i.*params.wavNum.*params.k_j(:,i)'*c');
    E_s = E_s + phase.*pEs(:,:,i);
    E_i = E_i + phase.*pEi(:,:,i);
end

%close(h)
%save 2d_computeEsEi.mat
E_s(params.psMask<params.a) = 0;%E_i(r<a);

E_i(params.psMask>=params.a) = 0;