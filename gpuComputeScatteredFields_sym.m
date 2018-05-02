function [E_s , E_i] = gpuComputeScatteredFields_sym(params)

kr = params.wavNum.*params.r_ps;
gpu_hl_kr = gpuArray(shank1(params.numOrd, kr, 'multiple'));

%The internal field specifying the field inside of a sphere
%for an incident plane wave

%compute the argument for the spherical bessel function k*n*r
knr = params.wavNum * params.n .* params.r_ps;

%compute the spherical bessel function
gpu_jl_knr = gpuArray(sphbesselj(params.numOrd,knr,'multiple'));


%h = waitbar(0, 'Monte Carlo integration...');


% r = sqrt(normKvec(1)^2 + normKvec(2)^2 + normKvec(3)^2);
% theta =  atan2(normKvec(2), normKvec(1));
% phi = acos(normKvec(2)/r);
% kSpherical = [r theta phi]

%the vector from the focal point to the center of the sphere
c = params.gpu_ps - params.gpu_pf;
% gpu_hl_kr = gpuArray(hl_kr);
% gpu_jl_knr = gpuArray(jl_knr);


scale = (1/params.samples).*params.subA;
for i=1:params.samples
    cos_theta = ((params.k_j(:,i)'*params.normPMinPs')');
%    cos_theta = reshape(cos_theta, params.simRes, params.simRes);
    %compute the legendre polynomials needed for computing E_s and E_i
    gpu_Pl_costheta = gpuLegendre(params.numOrd,cos_theta,params.P);
    
%     test_1 = rand(params.simRes,params.simRes,params.numOrd,'gpuArray');
%     test_2 = rand(params.simRes,params.simRes,params.numOrd,'gpuArray');
%     test  = test_1.*test_2;
    gpu_hlkr_Plcostheta = gpu_hl_kr.*gpu_Pl_costheta;
    
    %hlkr_Plcostheta = gather(gpu_hlkr_Plcostheta);
    
    
    %multiply by the legendre polynomial
    gpu_jlknr_Plcostheta = gpu_jl_knr.*gpu_Pl_costheta;
    
    %multiply by the scattering coefficients B
    %sum all orders
    phase = exp(1i.*params.wavNum.*params.k_j(:,i)'*c');
    params.gpu_Es = params.gpu_Es + scale.*phase.*sum(bsxfun(@times, gpu_hlkr_Plcostheta, params.B'),2);
    %E_s(gpu_r<params.a) = 0;%E_i(r<a);
    %subplot(1,2,1), imagesc((abs((E_s)))),title('E_s'), colorbar, axis image, colormap(brewer)
    params.gpu_E_i = params.gpu_E_i + scale.*phase.*sum(bsxfun(@times, gpu_jlknr_Plcostheta, params.A'),2);
    %E_i(gpu_r>params.a) = 0;
    %subplot(1,2,2), imagesc((abs((E_i)))),title('E_i'), colorbar, axis image, colormap(brewer)
    %pause(0.1)
    %waitbar(i/params.samples);
end

%close(h)

params.gpu_Es(params.psMask<params.a) = 0;%E_i(r<a);
E_s = params.gpu_Es;

params.gpu_E_i(params.psMask>params.a) = 0;
E_i = params.gpu_E_i;
