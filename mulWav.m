%% setup
% single plane wave case
clear
addpath(genpath('~/source/stim-matlab/'))

theta=0; phi=0;
kdir = [phi theta 1]
params.E0 = 1;
params.NA_in = 0;
params.NA_out = .3;

params.samples=200;
params.a = 1;                     %radius of the sphere
%wavelength
params.lam = 1;
params.numOrd = computeN_l(params.a, params.lam);
%params.numOrd = 50;
params.n=1.4;
%specify the size of the field plane in wavelength units (microns)
params.gridSize = 5;
%specify the spatial resolution of the field plane
params.spatRes = 256;

params.wavNum = 2*pi/params.lam;        %wavenumber

%direction of the incident light
[x,y,z] = sph2cart(theta,phi,1);

% r=1;
% x = r * cos(theta) * sin(phi); y = (r* sin(theta) * sin(phi));
% z = (r * cos(phi));


params.kVec=[x y z]*params.wavNum

%% get r and rVecs
%create a grid of points representing pixel positions in the field plane
gridPoints = linspace(-params.gridSize,params.gridSize,params.spatRes);
[x,y] = meshgrid(gridPoints, gridPoints);
z = zeros(params.spatRes,params.spatRes);   %field plane z = 0


%convert the field plane pixel positions to cartesian coordinates
[theta, phi, r] = cart2sph(x,y,z);



%create a parameter structure for the simulation
params.r=r;                         %r value at each pixel position

%store the point positions in rVecs
rVecs = zeros(params.spatRes*params.spatRes, 3);
rVecs(:,1) = x(:); rVecs(:,2) = y(:); rVecs(:,3) = z(:);

%% world space coordinates
% pMin=[-5 0 -5]; pMax = [-5 0 5]; normal = [5 0 5];
% size=256;
% [a,b] = meshgrid(linspace(0,1,size),linspace(0,1,size));
% idx=1;
% for i=1:size
%     for j=1:size
%         rVecs(idx,:) = worldSpace(a(i,j),b(i,j), pMin, pMax, normal);
%         params.r(i,j)=norm(rVecs(idx,:));
%         idx=idx+1;
%     end
% end


%% focused field
%apply the partial wave expansion to get the focused field

%params.kVec=k_j(:,1)'.*params.wavNum;

% this is needed for the plane wave expansion
orderEf=200;

normRvecs = bsxfun(@rdivide, rVecs,  sqrt(sum(rVecs.^2,2)));
normKvec = params.kVec./params.wavNum;
cos_theta = (normKvec*normRvecs')';
cos_theta = reshape(cos_theta, params.spatRes, params.spatRes);

%calculate the prefix term (2l + 1)*i^l
ordVecEf=(0:orderEf)';
il = 1i.^ordVecEf;

jl_kr = sphbesselj(orderEf,params.wavNum.*params.r,'multiple');
Pl_costheta = myLegendre(orderEf,cos_theta);
jlkr_Pcostheta = jl_kr.*Pl_costheta;

il_jlkr_Pcostheta = (bsxfun(@times, jlkr_Pcostheta, reshape(il,[1 1 orderEf+1])));
alpha1 = asin(params.NA_in); alpha2 = asin(params.NA_out);
Pl_cosalpha1 = myLegendre(orderEf+1,cos(alpha1));
Pl_cosalpha2 = myLegendre(orderEf+1,cos(alpha2));

ord=0;
il_jlkr_Pcostheta(:,:,ord+1) = il_jlkr_Pcostheta(:,:,ord+1).*(Pl_cosalpha1(ord+2)-Pl_cosalpha2(ord+2)-Pl_cosalpha1(1)+Pl_cosalpha2(1));

ord=1;
il_jlkr_Pcostheta(:,:,ord+1) = il_jlkr_Pcostheta(:,:,ord+1).*(Pl_cosalpha1(ord+2)-Pl_cosalpha2(ord+2)-Pl_cosalpha1(1)+Pl_cosalpha2(1));

for ord=2:orderEf
    il_jlkr_Pcostheta(:,:,ord+1) = il_jlkr_Pcostheta(:,:,ord+1).*(Pl_cosalpha1(ord+2)-Pl_cosalpha2(ord+2)-Pl_cosalpha1(ord)+Pl_cosalpha2(ord));
end

E_f = 2*pi*params.E0.*sum(il_jlkr_Pcostheta,3);
brewer = brewermap(1000);
figure, imagesc((abs((E_f)))),title('abs(E_f)'), colorbar, axis image
colormap(brewer)


%%
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


E_s = zeros(params.spatRes,params.spatRes);
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

for i=1:params.samples
    cos_theta = (k_j(:,i)'*normRvecs')';
    cos_theta = reshape(cos_theta, params.spatRes, params.spatRes);
    %compute the legendre polynomials needed for computing E_s and E_i
    Pl_costheta = myLegendre(params.numOrd,cos_theta);
    hlkr_Plcostheta = hl_kr.*Pl_costheta;
    
    %multiply by the legendre polynomial
    jlknr_Plcostheta = jl_knr.*Pl_costheta;
    
    %multiply by the scattering coefficients B
    %sum all orders
    
    E_s = E_s + (1/params.samples).*subA.*sum(bsxfun(@times, hlkr_Plcostheta, reshape(B,[1 1 params.numOrd+1])),3);
    %E_s(params.r<params.a) = 0;%E_i(r<a);
    %subplot(1,2,1), imagesc((abs((E_s)))),title('E_s'), colorbar, axis image, colormap(brewer)
    E_i = E_i + (1/params.samples).*subA.*sum(bsxfun(@times, jlknr_Plcostheta, reshape(A,[1 1 params.numOrd+1])),3);
    %E_i(params.r>params.a) = 0;
    %subplot(1,2,2), imagesc((abs((E_i)))),title('E_i'), colorbar, axis image, colormap(brewer)
    %pause(0.1)
    waitbar(i/params.samples);
end
close(h)

E_s(params.r<params.a) = 0;%E_i(r<a);
figure, imagesc((abs((E_s)))),title('E_s'), colorbar, axis image
colormap(brewer)

E_i(params.r>params.a) = 0;
figure, imagesc((abs((E_i)))),title('E_i'), colorbar, axis image
colormap(brewer)

E_t = zeros(params.spatRes, params.spatRes);
E_t(params.r<params.a) = E_i(params.r<params.a);
E_t(params.r>=params.a) = E_f(params.r>=params.a) + E_s(params.r>=params.a);
figure, imagesc((abs((E_t)))),title('E_t'), colorbar, axis image
colormap(brewer)