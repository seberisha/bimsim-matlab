function [cpu_Es , cpu_Ei] = cpuComputeScatteredFields(params,material)



n = material(2) + 1i*material(3);

numOrd = computeN_l(params.sphereRadius, material(1));

kr = params.magKVector.*params.r_ps;
cpu_hl_kr = (cpuShank1(numOrd, kr, 'multiple'));

%The internal field specifying the field inside of a sphere
%for an incident plane wave

%compute the argument for the spherical bessel function k*n*r
knr = params.magKVector * n .* params.r_ps;

%compute the spherical bessel function
cpu_jl_knr = (sphbesselj(numOrd,knr,'multiple'));

%the vector from the focal point to the center of the sphere
c = (params.ps) - (params.pf);
cpu_Es = zeros(params.simRes,params.simRes);
cpu_Ei = cpu_Es;

scale = (1/params.samples).*params.subA;
for i=1:params.samples
    cos_theta = ((params.k_j(:,i)'*params.normPMinPs')');
    cos_theta = reshape(cos_theta, params.simRes, params.simRes);
    %compute the legendre polynomials needed for computing E_s and E_i
    cpu_Pl_costheta = cpuLegendre(numOrd,cos_theta);
    
    cpu_hlkr_Plcostheta = cpu_hl_kr.*cpu_Pl_costheta;
    
    %multiply by the legendre polynomial
    cpu_jlknr_Plcostheta = cpu_jl_knr.*cpu_Pl_costheta;
    
    %multiply by the scattering coefficients B
    %sum all orders
    phase = exp(1i.*params.magKVector.*params.k_j(:,i)'*c');
    cpu_Es = cpu_Es + phase.*sum(bsxfun(@times, cpu_hlkr_Plcostheta, params.B),3);
    cpu_Ei = cpu_Ei + phase.*sum(bsxfun(@times, cpu_jlknr_Plcostheta, params.sphereRadius),3);
end

cpu_Es = scale.*cpu_Es;
cpu_Ei = scale.*cpu_Ei;
cpu_Es(params.psMask<params.sphereRadius) = 0;%E_i(r<a);

cpu_Ei(params.psMask>=params.sphereRadius) = 0;
