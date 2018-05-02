%% setup
%clear
addpath(genpath('~/source/stim-matlab'))

M=256; N=256; res=256;
gridIdx=15;
[x,y] = meshgrid(linspace(-gridIdx,gridIdx,res), linspace(-gridIdx,gridIdx,res));
z=zeros(size(x));
[theta, phi, r] = cart2sph(x, y, z);

lambda=1;
k = 2*pi/lambda;
NA_o=0.7;
f_c = (2*pi/lambda)*NA_o; S = 100; deltaU = 1/S;
res=256; lambda=1; 
alpha_1=0; NA=NA_o; E_0=1; a=1; n=1.4; c=0;  
%N_l = computeN_l(a, lambda);
N_l=10; 

alpha_2=asin(NA);
k=2*pi/lambda; kr = k*r; ka=k*a; kna=k*n*a; knr=k*n*r;

%measured data
A = rand(M,N);
%imaginary part of the complex refractive index
k_j=rand;

%z=ones(size(x))*1; 
%[x,y,z] = meshgrid(linspace(-gridIdx,gridIdx,res), linspace(-gridIdx,gridIdx,res),linspace(-gridIdx,gridIdx,res));

%% pre-compute stuff
%spherical Bessel functions for all orders

% phase shift
% c = p_s - p_f;
e_ikc = exp(1i*k.*c);

j_kr = sphbesselj(N_l,kr,'multiple');

% Legendre polynomials 
P_ct = myLegendre(N_l,cos(theta));

%spherical Bessel functions for all orders
j_ka = squeeze(sphbesselj(N_l,ka,'multiple'));

% first derivative of the spherical Bessel function of the first kind
j = sym('sqrt(1/2*pi/x)*besselj(n+1/2,x)');
dj = simplify(diff(j));
dj = vectorize(inline(char(dj),'n','x'));

deriv_j_kna = dj((0:N_l)', kna);

j_kna = squeeze(sphbesselj(N_l,kna,'multiple'));
deriv_j_ka = dj((0:N_l)', ka);
h_ka = shank1((0:N_l)',ka,'one');

% first derivative of the spherical Hankel function of the first kind
h = sym('((pi/(2*x))^.5)*besselj(n+.5,x)+i*((pi/(2*x))^.5)*bessely(n+.5,x)');
dh = simplify(diff(h));
dh = vectorize(inline(char(dh),'n','x'));

deriv_h_ka = dh((0:N_l)', ka);

B = ((2*(0:N_l)' + 1).*1i.^(0:N_l)' ).*(j_ka.*deriv_j_kna*n - j_kna.*deriv_j_ka)./(j_kna.*deriv_h_ka - h_ka.*deriv_j_kna*n);

h_kr = shank1(N_l,kr,'multiple');

deriv_j_kna = dj((0:N_l)', kna);

A = ((2*(0:N_l)' + 1).*1i.^(0:N_l)' ).*(j_ka.*deriv_h_ka - deriv_j_ka.*h_ka)./(j_kna.*deriv_h_ka - h_ka.*deriv_j_kna*n);

j_knr = squeeze(sphbesselj(N_l,knr,'multiple'));
%%
    
[E_t, E_f, E_s, E_i] = forwardModel_singlePlaneWave(N_l, E_0, a, r, j_kr, P_ct, B, h_kr, A, j_knr, e_ikc);

%create a colormap and display the image
brewer = brewermap(1000);
figure(1), imagesc(abs((E_f))),title('E_f'), colorbar, axis image
colormap(brewer);

%create a colormap and display the image
figure(2), imagesc(abs((E_s))), title('E_s'), colorbar, axis image
colormap(brewer);

figure(3), imagesc(abs((E_i))), title('E_i'), colorbar, axis image
colormap(brewer);

figure(4), imagesc(abs((E_t))), title('E_t'), colorbar, axis image
colormap(brewer);