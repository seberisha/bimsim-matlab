function [E_t, E_f, E_s, E_i] = forwardModel(res, lambda, N_l, alpha_1, NA, E_0, a, n, c, gridIdx )
%% setup

addpath(genpath('~/source/stim-matlab'))
%params
%res = 256; lambda=3; N_l = 100; alpha_1=0; NA=.9; E_0=1; a=3; n=3;
%c=0; gridIdx=20;
% [x,y, z] = meshgrid(linspace(-gridIdx,gridIdx,res), linspace(-gridIdx,gridIdx,res), linspace(-gridIdx,gridIdx,res));
% [theta, r] = cart2pol(x,y,z);

[x,y] = meshgrid(linspace(-gridIdx,gridIdx,res), linspace(-gridIdx,gridIdx,res));
[theta, r] = cart2pol(x,y);


alpha_2=asin(NA);
k=2*pi/lambda; kr = k*r; ka=k*a; kna=k*n*a; knr=k*n*r;
% phase shift
% c = p_s - p_f;
e_ikc = exp(1i*k.*c);

%% pre-compute j, p, h, coefficients
%spherical Bessel functions for all orders

j_kr = sphbesselj(N_l,kr,'multiple');

% Legendre polynomials 
P_ct = myLegendre(N_l,cos(theta));
P_ca1 = myLegendre(N_l+1, cos(alpha_1));
P_ca2= myLegendre(N_l+1, cos(alpha_2));

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

B = ((2*(0:N_l)' + 1).*1i.^(0:N_l)' ).*(j_ka.*deriv_j_kna - j_kna.*deriv_j_ka)./(j_kna.*deriv_h_ka - h_ka.*deriv_j_kna*n);

h_kr = shank1(N_l,kr,'multiple');

deriv_j_kna = dj((0:N_l)', kna);

A = ((2*(0:N_l)' + 1).*1i.^(0:N_l)' ).*(j_ka.*deriv_h_ka - deriv_j_ka.*h_ka)./(j_kna.*deriv_h_ka - h_ka.*deriv_j_kna*n);

j_knr = squeeze(sphbesselj(N_l,knr,'multiple'));



%% Compute E_f
E_f = computeEf(j_kr, P_ct, P_ca1, P_ca2, E_0, N_l);

%create a colormap and display the image
brewer = brewermap(1000);
figure(1), imagesc(abs((E_f))),title('E_f'), colorbar, axis image
colormap(brewer);

%% 
%Compute E_s
E_s = computeEs(B, h_kr, P_ct, N_l, e_ikc);

%Compute E_i
E_i = computeEi(A, j_knr, P_ct, e_ikc, N_l);

E_i(r>a) = 0;
E_s(r<a) = E_i(r<a);
E_i(abs(r-a)<=1e-8) = E_f(r==a) + E_s(r==a);
%create a colormap and display the image
brewer = brewermap(1000);
figure(2), imagesc(abs((E_s))), title('E_s'), colorbar, axis image
colormap(brewer);

brewer = brewermap(1000);
figure(3), imagesc(abs((E_i))), title('E_i'), colorbar, axis image
colormap(brewer);

E_t = E_s + E_f;

brewer = brewermap(1000);
figure(4), imagesc(abs((E_t))), title('E_t'), colorbar, axis image
colormap(brewer);



