%% setup
clear

addpath(genpath('~/source/stim-matlab/'))
res=256;
gridIdx=15;

% r=ones(res,res);
% theta = meshgrid(linspace(0,2*pi,res),linspace(0,2*pi,res));
% [x,y,z] = sph2cart(theta,0,r);

[x,y] = meshgrid(linspace(-gridIdx,gridIdx,res), linspace(-gridIdx,gridIdx,res));
z=zeros(size(x));
[theta, phi, r] = cart2sph(x, y, z);
%wavelength of incident light
lambda=1;
NA_in = 0;
NA_out= 1;
S = 100; deltaU = 1/S;
alpha_1=asin(NA_in);
E_0=1;
a=1;
n=1.4;% + 1i*.05; 
c=zeros(3,1);
N_l = computeN_l(a, lambda);
alpha_2=asin(NA_out);
%magnitude of the k-vector
kmag=2*pi/lambda; kr = kmag*r; 
%the arguments k*a and k*n*a
ka=kmag*a; kna=kmag*n*a; 
knr=kmag*n*r;
k_vec=[1 0 0];

%% pre-compute stuff

%spherical Bessel functions for all orders
j_ka = squeeze(sphbesselj(N_l,ka,'multiple'));

% first derivative of the spherical Bessel function of the first kind
j = sym('sqrt(1/2*pi/x)*besselj(n+1/2,x)');
dj = simplify(diff(j));
dj = vectorize(inline(char(dj),'n','x'));

%deriv_j_kna = dricbesj((0:N_l)', kna);

deriv_j_kna = dj((0:N_l)', kna);

j_kna = squeeze(sphbesselj(N_l,kna,'multiple'));

%deriv_j_ka=dricbesj((0:N_l)', ka);

deriv_j_ka = dj((0:N_l)', ka);

%Hankel functions for ka
h_ka = shank1((0:N_l)',ka,'one');

% first derivative of the spherical Hankel function of the first kind
h = sym('((pi/(2*x))^.5)*besselj(n+.5,x)+i*((pi/(2*x))^.5)*bessely(n+.5,x)');
dh = simplify(diff(h));
dh = vectorize(inline(char(dh),'n','x'));

deriv_h_ka = dh((0:N_l)', ka);
h_kr = shank1(N_l,kr,'multiple');
deriv_j_kna = dj((0:N_l)', kna);
j_knr = squeeze(sphbesselj(N_l,knr,'multiple'));

%A (internal scattering coefficient)
%numerator
num_a = j_ka.*deriv_h_ka - deriv_j_ka.*h_ka;
den = j_kna.*deriv_h_ka - h_ka.*deriv_j_kna*n;
A = ((2*(0:N_l)' + 1).*1i.^((0:N_l)') ).*(num_a./den);

%B (external scattering coefficient)
num_b=j_ka.*deriv_j_kna*n - j_kna.*deriv_j_ka;

B = ((2.*(0:N_l)' + 1).*1i.^((0:N_l)') ).*( num_b./den);

%%                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
rVecs = zeros(3,res*res);
idx=1;

for j=1:res
    for i=1:res
        temp = [x(i,j); y(i,j); z(i,j)];
        rVecs(:,idx)=temp./norm(temp);
        idx=idx+1;
    end
end

%%
Nl_Ef = 200;

%spherical Bessel functions
j_kr = squeeze(sphbesselj(Nl_Ef,kr,'multiple'));

% Legendre polynomials
P_ca1 = myLegendre(Nl_Ef+1, cos(alpha_1));
P_ca2= myLegendre(Nl_Ef+1, cos(alpha_2));


theta_Ef = k_vec*rVecs;
cos_theta = reshape(theta_Ef, res, res);
% Legendre polynomials
P_ct_ef = myLegendre(Nl_Ef,cos_theta);

brewer = brewermap(1000);

E_f = computeEf(j_kr, P_ct_ef, P_ca1, P_ca2, E_0, Nl_Ef);
figure, imagesc((abs((E_f)))),title('E_f'), colorbar, axis image
colormap(brewer)

%%
NW=200;
NA_in=0; NA_out=1;
k_j = monteCarlo(NW,k_vec, NA_in, NA_out);
E_s = zeros(res,res);
E_i = E_s;
h = waitbar(0, 'Monte Carlo integration...');
 %compute the amplitude that makes it through the condenser
amplitude=1;
subA = 2 * pi * amplitude * ( (1 - cos(asin(alpha_2))) - (1 - cos(asin(alpha_1))) );
for i=1:NW
    kRvecs=(k_j(:,i)'*rVecs);
    cos_theta=reshape(kRvecs, res,res);
    % phase shift
    % c = p_s - p_f;
    % e_ikc = exp(1i*k_j(:,i)'*c);
    
    % Legendre polynomials
    P_ct = myLegendre(N_l,cos_theta);
    
    E_s = E_s + (1/NW).*subA.*computeEs(B, h_kr, P_ct, N_l, 1);
    
    %Compute E_i
    E_i = E_i + (1/NW).*subA.*computeEi(A, j_knr, P_ct, 1, N_l);
    
    E_i(r>a) = 0;
    E_s(r<a) = E_i(r<a);
    %E_i(abs(r-a)<=1e-2) = E_f(abs(r-a)<=1e-2) + E_s(abs(r-a)<=1e-2);
    waitbar(i/(NW+1), h)
end
% E_s = 2*pi*(1-cos(alpha_2))*E_s;
% E_i = 2*pi*(1-cos(alpha_2))*E_i;
close(h)

figure
subplot(1,2,1), imagesc((abs(E_s))),title('E_s'), colorbar, axis image
subplot(1,2,2), imagesc((abs(E_i))),title('E_i'), colorbar, axis image
colormap(brewer)
E_t = zeros(res,res);
E_t(r>=a) = E_s(r>=a) + E_f(r>=a);
E_t(r<a)=E_i(r<a);
figure, imagesc((abs(E_t))),title('E_t'), colorbar, axis image
colormap(brewer)