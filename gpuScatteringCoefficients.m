function [A, B] = gpuScatteringCoefficients(material,a)

%refractive index
n = material(2) + 1i*material(3);

%wavelength in micrometers
lambda = material(1);

wavNum = 2*pi/lambda;        %wavenumber

numOrdForMC = computeN_l(a, lambda);

%create a vector of orders [0 1 ... Nl]
ordVec = gpuArray((0:numOrdForMC)');

%calculate the prefix term (2l + 1)*i^l
twolp1 = 2.*ordVec+1;
il = 1i.^ordVec;
twolp1_il = twolp1.*il;

%compute the arguments needed to evaluate spherical bessel functions,
%hankel functions, and their derivatives
ka = wavNum*a;
kna = wavNum*n*a;

%evaluate the spherical bessel functions of the first kind at ka
jl_ka = gpuArray(sphbesselj(numOrdForMC,ka,'multiple'));

%evaluate the derivate of the spherical bessel functions of the first kind at kna
jl_kna_p = gpuArray(derivSphBes(numOrdForMC, kna));

%evaluate the spherical bessel functions of the first kind at kna
jl_kna = gpuArray(sphbesselj(numOrdForMC,kna,'multiple'));

%evaluate the derivative of the spherical bessel functions of the first kind at ka
jl_ka_p = gpuArray(derivSphBes(numOrdForMC, ka));

%compute the numerator for B coefficients
numB = jl_ka.*jl_kna_p.*n - jl_kna.*jl_ka_p;

%evaluate the derivative of the hankel functions of the first kind at ka
hl_ka_p = gpuArray(derivSphHan(numOrdForMC, ka));

%evaluate the hankel functions of the first kind at ka
hl_ka = gpuArray(shank1(numOrdForMC, ka, 'multiple'));

%compute the denominator for coefficients A and B
denAB = jl_kna.*hl_ka_p - hl_ka.*jl_kna_p*n;

%compute B
B = twolp1_il.*(numB./denAB);
B = reshape(B,[1 1 numOrdForMC+1]);
%calculate the numerator for the scattering coefficients A
numA = jl_ka.*hl_ka_p - jl_ka_p.*hl_ka;
%calculate the scattering coefficients A
A = twolp1_il.*(numA./denAB);
A = reshape(A,[1 1 numOrdForMC+1]);
