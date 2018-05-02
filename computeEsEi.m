function [E_s , E_i] = computeEsEi(params)

%normalize the position vectors
normRvecs = bsxfun(@rdivide, params.rVecs,  sqrt(sum(params.rVecs.^2,2)));

% compute alpha1 and alpha2 from NA_in and NA_out, respectively
alpha1 = asin(params.NA_in); alpha2 = asin(params.NA_out);


%create a vector of orders [0 1 ... Nl]
ordVec = (0:params.numOrd)';

%The scattered field for a single incident plane-wave k produced by
%a sphere with radius r positioned at point pf

%calculate the prefix term (2l + 1)*i^l
twolp1 = 2.*ordVec+1;
il = 1i.^ordVec;
twolp1_il = twolp1.*il;

%compute the arguments needed to evaluate spherical bessel functions,
%hankel functions, and their derivatives
ka=params.wavNum*params.a;
kr = params.wavNum.*params.r;
kna = params.wavNum*params.n*params.a;

%evaluate the spherical bessel functions of the first kind at ka
jl_ka = sphbesselj(params.numOrd,ka,'multiple');
%evaluate the derivate of the spherical bessel functions of the first kind at kna
jl_kna_p = derivSphBes(params.numOrd, kna);
%evaluate the spherical bessel functions of the first kind at kna
jl_kna = sphbesselj(params.numOrd,kna,'multiple');
%evaluate the derivative of the spherical bessel functions of the first kind at ka
jl_ka_p = derivSphBes(params.numOrd, ka);

%compute the numerator for B coefficients
numB = jl_ka.*jl_kna_p.*params.n - jl_kna.*jl_ka_p;
%evaluate the derivative of the hankel functions of the first kind at ka
hl_ka_p = derivSphHan(params.numOrd, ka);
%evaluate the hankel functions of the first kind at ka
hl_ka = shank1(params.numOrd, ka, 'multiple');
%compute the denominator for coefficients A and B
denAB = jl_kna.*hl_ka_p - hl_ka.*jl_kna_p*params.n;

%compute B
B = twolp1_il.*(numB./denAB);
hl_kr = shank1(params.numOrd, kr, 'multiple');

%The internal field specifying the field inside of a sphere
%for an incident plane wave

%compute the argument for the spherical bessel function k*n*r
knr = params.wavNum * params.n .* params.r;

%compute the spherical bessel function
jl_knr = sphbesselj(params.numOrd,knr,'multiple');



%calculate the numerator for the scattering coefficients A
numA = jl_ka.*hl_ka_p - jl_ka_p.*hl_ka;
%calculate the scattering coefficients A
A = twolp1_il.*(numA./denAB);


E_s = zeros(params.simRes,params.simRes);
E_i = E_s;
h = waitbar(0, 'Monte Carlo integration...');
%compute the amplitude that makes it through the condenser
subA = 2 * pi * params.E0 * ( (1 - cos(alpha2)) - (1 - cos(alpha1)) );

params.samples=200

% r = sqrt(normKvec(1)^2 + normKvec(2)^2 + normKvec(3)^2);
% theta =  atan2(normKvec(2), normKvec(1));
% phi = acos(normKvec(2)/r);
% kSpherical = [r theta phi]
k_j = monteCarlo(params.samples, params.kVec./params.wavNum, params.NA_in, params.NA_out);

%the vector from the focal point to the center of the sphere
c = params.ps - params.pf;
for i=1:params.samples
    cos_theta = (k_j(:,i)'*normRvecs')';
    cos_theta = reshape(cos_theta, params.simRes, params.simRes);
    %compute the legendre polynomials needed for computing E_s and E_i
    Pl_costheta = myLegendre(params.numOrd,cos_theta);
    hlkr_Plcostheta = hl_kr.*Pl_costheta;
    
    %multiply by the legendre polynomial
    jlknr_Plcostheta = jl_knr.*Pl_costheta;
    
    %multiply by the scattering coefficients B
    %sum all orders
    phase = exp(1i.*params.kVec*c');
    E_s = E_s + (1/params.samples).*subA.*phase.*sum(bsxfun(@times, hlkr_Plcostheta, reshape(B,[1 1 params.numOrd+1])),3);
    %E_s(params.r<params.a) = 0;%E_i(r<a);
    %subplot(1,2,1), imagesc((abs((E_s)))),title('E_s'), colorbar, axis image, colormap(brewer)
    E_i = E_i + (1/params.samples).*subA.*phase.*sum(bsxfun(@times, jlknr_Plcostheta, reshape(A,[1 1 params.numOrd+1])),3);
    %E_i(params.r>params.a) = 0;
    %subplot(1,2,2), imagesc((abs((E_i)))),title('E_i'), colorbar, axis image, colormap(brewer)
    %pause(0.1)
    waitbar(i/params.samples);
end
close(h)

E_s(params.r<params.a) = 0;%E_i(r<a);

E_i(params.r>params.a) = 0;