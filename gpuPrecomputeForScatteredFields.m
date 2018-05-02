function [E_s , E_i] = gpuPrecomputeForScatteredFields(params)

kr = params.wavNum.*params.r_ps;
gpu_hl_kr = gpuArray(gpuShank1(params.numOrd, kr, 'multiple'));

%The internal field specifying the field inside of a sphere
%for an incident plane wave

%compute the argument for the spherical bessel function k*n*r
knr = params.wavNum * params.n .* params.r_ps;

%compute the spherical bessel function
gpu_jl_knr = gpuArray(sphbesselj(params.numOrd,knr,'multiple'));

scale = (1/params.samples).*params.subA;
E_s = gpuArray(zeros(params.simRes,params.simRes, params.samples));
E_i = E_s;

for i=1:params.samples
    cos_theta = ((params.k_j(:,i)'*params.normPMinPs')');
    cos_theta = reshape(cos_theta, params.simRes, params.simRes);
    %compute the legendre polynomials needed for computing E_s and E_i
    gpu_Pl_costheta = gpuLegendre(params.numOrd,cos_theta,params.P);
    
    gpu_hlkr_Plcostheta = gpu_hl_kr.*gpu_Pl_costheta;
    
    %multiply by the legendre polynomial
    gpu_jlknr_Plcostheta = gpu_jl_knr.*gpu_Pl_costheta;
    
    %multiply by the scattering coefficients B
    %sum all orders
    E_s(:,:,i) = scale.*sum(bsxfun(@times, gpu_hlkr_Plcostheta, params.B),3);
    E_i(:,:,i) = scale.*sum(bsxfun(@times, gpu_jlknr_Plcostheta, params.A),3);
end