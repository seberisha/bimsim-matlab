%% setup
% single plane wave case
clear

theta=0; phi=0;

addpath(genpath('~/source/stim-matlab/'))
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

[x,y,z] = sph2cart(phi, theta, 1);

params.kVec=[x y z].*params.wavNum

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

%%



%params.kVec = [0.6761 -0.4546 0.5798].*params.wavNum

%params.kVec=[1 0 0] .* params.wavNum;%direction of the incident light

%calculate the dot product between k and r
k_dot_r = (params.kVec*rVecs')';
k_dot_r = reshape(k_dot_r, params.spatRes, params.spatRes);

%create a vector of orders [0 1 ... Nl]
params.ordVec = 0:params.numOrd;
ordVec = (0:params.numOrd)';

params.E0 = 1;                      %amplitude of the plane wave

%% focused field

%single plane wave:
E_f = params.E0.*exp(1i.*k_dot_r);

brewer = brewermap(1000);
figure, imagesc((abs((E_f)))),title('abs(E_f)'), colorbar, axis image
colormap(brewer)

%%
%The scattered field for a single incident plane-wave k produced by
%a sphere with radius r positioned at point pf

%calculate the prefix term (2l + 1)*i^l
twolp1 = 2.*ordVec+1;
il = 1i.^ordVec;
twolp1_il = twolp1.*il;

%normalize rVecs to have norm 1
normRvecs = bsxfun(@rdivide, rVecs,  sqrt(sum(rVecs.^2,2)));
normKvec = params.kVec./params.wavNum;
cos_theta = (normKvec*normRvecs')';
cos_theta = reshape(cos_theta, params.spatRes, params.spatRes);

% temp = myLegendre(params.numOrd+1, cos_theta);
% 
% sqrt_one_m_xsq = sqrt(1-cos_theta.^2);
% xSq_m_one = cos_theta.^2 -1;
% 
% for i=0:params.numOrd
%     Pl_costheta(:,:,i+1) = -sqrt_one_m_xsq.*(-((i+1).*temp(:,:,i+1) - temp(:,:,i+2))./(xSq_m_one));
% end


% temp = legendre(params.numOrd, cos_theta);
% for i=1:params.numOrd+1
%  Pl_costheta(:,:,i) = squeeze(temp(i,:,:));
% end
%compute the legendre polynomials needed for computing E_s and E_i
Pl_costheta = myLegendre(params.numOrd,cos_theta);

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
% for i=0:params.numOrd
%    Pl_costheta(:,:,i+1)=legendreP(i,cos_theta); 
% end
hlkr_Plcostheta = hl_kr.*Pl_costheta;

%multiply by the scattering coefficients B
E_s = bsxfun(@times, hlkr_Plcostheta, reshape(B,[1 1 params.numOrd+1]));

%sum all orders
E_s = sum(E_s,3);
E_s(params.r<params.a) = 0;%E_i(r<a);
figure, imagesc((abs((E_s)))),title('E_s'), colorbar, axis image
colormap(brewer)

%%
%The internal field specifying the field inside of a sphere
%for an incident plane wave 

%compute the argument for the spherical bessel function k*n*r
knr = params.wavNum * params.n .* params.r;

%compute the spherical bessel function
jl_knr = sphbesselj(params.numOrd,knr,'multiple');

%multiply by the legendre polynomial
jlknr_Plcostheta = jl_knr.*Pl_costheta;

%calculate the numerator for the scattering coefficients A
numA = jl_ka.*hl_ka_p - jl_ka_p.*hl_ka;
%calculate the scattering coefficients A
A = twolp1_il.*(numA./denAB);

E_i = bsxfun(@times, jlknr_Plcostheta, reshape(A,[1 1 params.numOrd+1]));

E_i = sum(E_i,3);
E_i(params.r>params.a) = 0;
figure, imagesc((abs((E_i)))),title('E_i'), colorbar, axis image
colormap(brewer)

E_t = zeros(params.spatRes, params.spatRes);
E_t(params.r<params.a) = E_i(params.r<params.a);
E_t(params.r>=params.a) = E_f(params.r>=params.a) + E_s(params.r>=params.a);
figure, imagesc((abs((E_t)))),title('E_t'), colorbar, axis image
colormap(brewer)