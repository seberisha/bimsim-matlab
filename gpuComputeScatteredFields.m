function [gpu_Es , gpu_Ei] = gpuComputeScatteredFields(params,material)



n = material(2) + 1i*material(3);

numOrd = computeN_l(params.sphereRadius, material(1));

kr = params.magKVector.*params.r_ps;
gpu_hl_kr = gpuArray(gpuShank1(numOrd, kr, 'multiple'));

%The internal field specifying the field inside of a sphere
%for an incident plane wave

%compute the argument for the spherical bessel function k*n*r
knr = params.magKVector * n .* params.r_ps;

%compute the spherical bessel function
gpu_jl_knr = gpuArray(sphbesselj(numOrd,knr,'multiple'));

%the vector from the focal point to the center of the sphere
c = gpuArray(params.ps) - gpuArray(params.pf);
gpu_Es = zeros(params.simRes,params.simRes,'gpuArray');
gpu_Ei = gpu_Es;

scale = (1/params.samples).*params.subA;
for i=1:params.samples
    cos_theta = ((params.k_j(:,i)'*params.normPMinPs')');
    cos_theta = reshape(cos_theta, params.simRes, params.simRes);
    %compute the legendre polynomials needed for computing E_s and E_i
    gpu_Pl_costheta = gpuLegendre(numOrd,cos_theta);
    
    gpu_hlkr_Plcostheta = gpu_hl_kr.*gpu_Pl_costheta;
    
    %multiply by the legendre polynomial
    gpu_jlknr_Plcostheta = gpu_jl_knr.*gpu_Pl_costheta;
    
    %multiply by the scattering coefficients B
    %sum all orders
    phase = exp(1i.*params.magKVector.*params.k_j(:,i)'*c');
    gpu_Es = gpu_Es + scale.*phase.*sum(bsxfun(@times, gpu_hlkr_Plcostheta, params.B),3);
    gpu_Ei = gpu_Ei + scale.*phase.*sum(bsxfun(@times, gpu_jlknr_Plcostheta, params.sphereRadius),3);
end

gpu_Es(params.psMask<params.sphereRadius) = 0;%E_i(r<a);

gpu_Ei(params.psMask>=params.sphereRadius) = 0;
