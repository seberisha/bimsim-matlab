function [E_s , E_i] = computeScatteredFields_sym(params)
%normalize the position vectors
%normRvecs = bsxfun(@rdivide, params.rVecs,  sqrt(sum(params.rVecs.^2,2)));


kr = params.wavNum.*params.normPMinPs;
hl_kr = shank1(params.numOrd, kr, 'multiple');

%The internal field specifying the field inside of a sphere
%for an incident plane wave

%compute the argument for the spherical bessel function k*n*r
knr = params.wavNum * params.n .* params.normPMinPs;

%compute the spherical bessel function
jl_knr = sphbesselj(params.numOrd,knr,'multiple');

E_s = zeros(params.simRes,1);
E_i = E_s;
%h = waitbar(0, 'Monte Carlo integration...');


% r = sqrt(normKvec(1)^2 + normKvec(2)^2 + normKvec(3)^2);
% theta =  atan2(normKvec(2), normKvec(1));
% phi = acos(normKvec(2)/r);
% kSpherical = [r theta phi]

%the vector from the focal point to the center of the sphere
c = params.ps - params.pf;
for i=1:params.samples
    cos_theta = (params.k_j(:,i)'*params.normPsVecs')';
    %compute the legendre polynomials needed for computing E_s and E_i
    Pl_costheta = myLegendre(params.numOrd,cos_theta);
    hlkr_Plcostheta = hl_kr.*Pl_costheta;

    %multiply by the legendre polynomial
    jlknr_Plcostheta = jl_knr.*Pl_costheta;
    
    %multiply by the scattering coefficients B
    %sum all orders
    phase = exp(1i.*params.wavNum.*params.k_j(:,i)'*c');
    E_s = E_s + (1/params.samples).*params.subA.*phase.*sum(bsxfun(@times, hlkr_Plcostheta, params.B'),2);
    %E_s(params.r_ps<params.a) = 0;%E_i(r<a);
    %subplot(1,2,1), imagesc((abs((E_s)))),title('E_s'), colorbar, axis image, colormap(brewer)
    E_i = E_i + (1/params.samples).*params.subA.*phase.*sum(bsxfun(@times, jlknr_Plcostheta, params.A'),2);
    %E_i(params.r_ps>params.a) = 0;
    %subplot(1,2,2), imagesc((abs((E_i)))),title('E_i'), colorbar, axis image, colormap(brewer)
    %pause(0.1)
    %waitbar(i/params.samples);
end

%close(h)

E_s(params.psMask<params.a) = 0;%E_i(r<a);

E_i(params.psMask>=params.a) = 0;