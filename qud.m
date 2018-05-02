%% setup
clear
res = 128;
[x,y] = meshgrid(linspace(-20,20,res), linspace(-20,20,res));
[theta, r] = cart2pol(x,y);

lambda=3;
k=2*pi/lambda;
kr = k*r;
sr = 1;
N_l = maxOrder(sr, lambda);

E_f = zeros(res,res);

alpha_1=0;

NA=.7;
alpha_2=asin(NA);
E_0=1;


%% method 6

% Legendre polynomials 
P_ct = myLegendre(N_l,cos(theta));
P_ca1 = myLegendre(N_l+1, cos(alpha_1).*ones(size(theta)));
P_ca2= myLegendre(N_l+1, cos(alpha_2).*ones(size(theta)));

for ord=0:N_l
    j_kr = sphbesselj(ord,kr,'one');
    i_c = 1i^ord;
    legMat = P_ca1(:,:,ord+2) - P_ca2(:,:,ord+2) - P_ca1(:,:,ord+1) + P_ca2(:,:,ord+1);
    E_f = E_f + i_c.*j_kr.*squeeze(P_ct(:,:,ord+1)).*squeeze(legMat); 
end

E_f = 2*pi*E_0*E_f;
figure(6), imagesc(abs((E_f))),colorbar,colorbar
