function [E_s, E_i] = precomputeForExtendedSource(params)

kr = params.wavNum.*params.r_ps;

hl_kr = shank1(params.numOrd, kr, 'multiple');

%The internal field specifying the field inside of a sphere
%for an incident plane wave

%compute the argument for the spherical bessel function k*n*r
%knr = params.wavNum * params.n .* params.r;
knr = params.wavNum * params.n .* params.r_ps;

%compute the spherical bessel function
jl_knr = sphbesselj(params.numOrd,knr,'multiple');

E_s = zeros(params.simRes,params.simRes,params.samples);
E_i = E_s;

for i=1:params.samples
    cos_theta = (params.k_j(:,i)'*params.normPMinPs')';
    cos_theta = reshape(cos_theta, params.simRes, params.simRes);
    %compute the legendre polynomials needed for computing E_s and E_i
    Pl_costheta = myLegendre(params.numOrd,cos_theta);
    hlkr_Plcostheta = hl_kr.*Pl_costheta;
    
    %multiply by the legendre polynomial
    jlknr_Plcostheta = jl_knr.*Pl_costheta;
    
    %multiply by the scattering coefficients B
    %sum all orders
    E_s(:,:,i) = (1/params.samples).*params.subA.*sum(bsxfun(@times, hlkr_Plcostheta, reshape(params.B,[1 1 params.numOrd+1])),3);
    E_i(:,:,i) = (1/params.samples).*params.subA.*sum(bsxfun(@times, jlknr_Plcostheta, reshape(params.A,[1 1 params.numOrd+1])),3);
end