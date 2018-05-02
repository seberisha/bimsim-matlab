%% setup
clear
addpath(genpath('~/source/stim-matlab/'))
res=256;
k=[1 0 0];
gridIdx=15;

% [x,y] = meshgrid(linspace(-gridIdx,gridIdx,res), linspace(-gridIdx,gridIdx,res));
% z=zeros(size(x));
% [theta, phi, r] = cart2sph(x, y, z);

r=ones(res,1);
theta = linspace(0,2*pi,res)';

x = r .* cos(zeros(size(r))) .* cos(theta);
y = r .* cos(zeros(size(r))) .* sin(theta);
z = r .* sin(zeros(size(r)));

%wavelength of the incident field
lambda=1;
NA_in = 0;
NA_out= 1;
S = 100; deltaU = 1/S;
alpha_1=asin(NA_in);
E_0=1;
a=1;
%complex refractive index of the sphere
n=1.4; 
c=zeros(3,1);
N_l = computeN_l(a, lambda);
%N_l=200;
alpha_2=asin(NA_out);
%magnitude of the k-vector
kmag=2*pi/lambda;
kr = kmag*r; ka=kmag*a; kna=kmag*n*a; knr=kmag*n*r;
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

h_ka = shank1((0:N_l)',ka,'one');

% first derivative of the spherical Hankel function of the first kind
h = sym('((pi/(2*x))^.5)*besselj(n+.5,x)+i*((pi/(2*x))^.5)*bessely(n+.5,x)');
dh = simplify(diff(h));
dh = vectorize(inline(char(dh),'n','x'));

deriv_h_ka = dh((0:N_l)', ka);

h_kr = squeeze(shank1(N_l,kr,'multiple'));

deriv_j_kna = dj((0:N_l)', kna);

j_knr = squeeze(sphbesselj(N_l,knr,'multiple'));

%A (internal scattering coefficient)
%numerator
num_a = j_ka.*deriv_h_ka - deriv_j_ka.*h_ka;
den = j_kna.*deriv_h_ka - h_ka.*deriv_j_kna*n;
A = ((2*(0:N_l)' + 1).*1i.^(0:N_l)' ).*(num_a./den);

%B (external scattering coefficient)
num_b=j_ka.*deriv_j_kna*n - j_kna.*deriv_j_ka;

B = ((2.*(0:N_l)' + 1).*1i.^(0:N_l)' ).*( num_b./den);


%%
%get the point in world space and then the r vector
% pMin=[-5, 0, -5];
% pMax=[5, 0, 5];
% normal=[0, 1, 0];
% 
% b=pMax;c=normal;
% 
% Y = b - pMin;
% X = c - pMin - Y;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     


% rVecs = zeros(3,res*res);
% idx=1;
% 
% 
% for j=1:res
%     for i=1:res
%         temp = [x(i,j); y(i,j); z(i,j)];
%         rVecs(:,idx)=temp./norm(temp);
%         idx=idx+1;
%     end
% end

for j=1:res
    rVecs(:,j)=[x(j);y(j);z(j)];
    rVecs(:,j)=rVecs(:,j)/norm(rVecs(:,j));
end

% pMin=[-15; 0; -15];
% pMax=[15; 0; 15];
% normal=[0; 1; 0];
% b=pMax;c=normal;
% 
% Y = b - pMin;
% X = c - pMin - Y;
% 
% 
% idx=1;
% for j=1:res
%     for i=1:res
%         su=x(i,j)/res;
%         sv=y(i,j)/res;
%         test(:,idx)=[su;sv;0];
%         rVecs(:,idx) = pMin + X * su + Y * sv;
%         rVecs(:,idx)=rVecs(:,idx)./norm(rVecs(:,idx));
%         idx=idx+1;
%     end
% end



%%
Nl_Ef = 100;
%spherical Bessel functions for all orders

j_kr = squeeze(sphbesselj(Nl_Ef,kr,'multiple'));

% Legendre polynomials
P_ca1 = myLegendre(Nl_Ef+1, cos(alpha_1));
P_ca2= myLegendre(Nl_Ef+1, cos(alpha_2));


theta_Ef = k_vec*rVecs;
%cos_theta = reshape(theta_Ef, res, res);
% Legendre polynomials
%P_ct_ef = myLegendre(Nl_Ef,cos_theta);
P_ct_ef = myLegendre(Nl_Ef,theta_Ef);

E_f = computeEf_vec(j_kr, P_ct_ef, P_ca1, P_ca2, E_0, Nl_Ef);
figure, plot((abs((E_f)))),title('E_f')

%%
NW=1000;
NA_in=0; NA_out=1;
k_j = monteCarlo(NW,k_vec, NA_in, NA_out);
E_s = zeros(res,1);
E_i = E_s;
h = waitbar(0, 'Monte Carlo integration...');
 %compute the amplitude that makes it through the condenser
amplitude=1;
subA = 2 * pi * amplitude * ( (1 - cos(asin(alpha_2))) - (1 - cos(asin(alpha_1))) );
for i=1:NW
    kRvecs=(k_j(:,i)'*rVecs);
    cos_theta=kRvecs;
    % phase shift
    % c = p_s - p_f;
    % e_ikc = exp(1i*k_j(:,i)'*c);
    
    % Legendre polynomials
    P_ct = myLegendre(N_l,cos_theta);
    
    E_s = E_s + (1/NW)*subA.*computeEs_vec(B, h_kr, P_ct, N_l, 1);
    
    %Compute E_i
    E_i = E_i + (1/NW)*subA.*computeEi_vec(A, j_knr, P_ct, 1, N_l);
    
    E_i(r>a) = 0;
    E_s(r<a) = E_i(r<a);
    %E_i(abs(r-a)<=1e-2) = E_f(abs(r-a)<=1e-2) + E_s(abs(r-a)<=1e-2);
    waitbar(i/(NW+1), h)
end
close(h)

figure
subplot(1,2,1), plot((abs(E_s))),title('E_s')
subplot(1,2,2), plot((abs(E_i))),title('E_i')
E_t = zeros(res,1);
E_t(r>=a) = E_s(r>=a) + E_f(r>=a);
E_t(r<a)=E_i(r<a);
figure, plot((abs(E_t))),title('E_t')