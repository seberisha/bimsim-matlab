%%
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


% P_ct = legendre(N_l,cos(theta));
% P_ca1 = legendre(N_l+1, cos(alpha_1));%.*ones(size(theta)));
% P_ca2= legendre(N_l+1, cos(alpha_2));%.*ones(size(theta)));
%%

P_ct = myLegendre(N_l,cos(theta));
P_ca1 = myLegendre(N_l+1, cos(alpha_1));%.*ones(size(theta)));
P_ca2= myLegendre(N_l+1, cos(alpha_2));%.*ones(size(theta)));
%%
j_kr = sphbesselj(N_l,kr);
allOrders=j_kr.*P_ct;
ord=0;
allOrders(ord+1,:,:) = 1i^ord.*allOrders(ord+1,:,:).*(P_ca1(ord+2)-P_ca2(ord+2)-P_ca1(1)+P_ca2(1));

for ord=1:N_l
    allOrders(ord+1,:,:) = 1i^ord.*allOrders(ord+1,:,:).*(P_ca1(ord+2)-P_ca2(ord+2)-P_ca1(ord)+P_ca2(ord));
end
E_f = squeeze(sum(allOrders,1));
E_0=1;
E_f = 2*pi*E_0*E_f;
imagesc(abs((E_f)))

%%
% P_ca1_lm1 = zeros(size(P_ct));
% P_ca1_lm1(1,:,:) = P_ca1(1,:,:);
% P_ca1_lm1(2:end,:,:) =  P_ca1(1:end-2,:,:);
% P_ca2_lm1 = zeros(size(P_ct));
% P_ca2_lm1(1,:,:) =P_ca2(1,:,:);
% P_ca2_lm1(2:end,:,:) =  P_ca2(1:end-2,:,:);

% P_ca1_lm1 = zeros(size(P_ct));
% P_ca1_lm1(1,:,:) = P_ca1(1,:,:);
% P_ca1_lm1(2,:,:) = P_ca1(1,:,:);
% P_ca1_lm1(3:end,:,:) =  P_ca1(2:end-2,:,:);
% P_ca2_lm1 = zeros(size(P_ct));
% P_ca2_lm1(1,:,:) =P_ca2(1,:,:);
% P_ca2_lm1(2,:,:) =P_ca2(1,:,:);
% P_ca2_lm1(3:end,:,:) =  P_ca2(2:end-2,:,:);


% legMat = P_ca1(2:end,:,:) - P_ca2(2:end,:,:) - P_ca1_lm1 + P_ca2_lm1;
% j_kr = sphbesselj(N_l,kr);
% allOrders = j_kr.*P_ct.*legMat; 
% for ord=0:N_l
%    E_f = E_f+ 1i^ord.*squeeze(allOrders(ord+1,:,:));
% end
% ord=0;
% j_kr = sphbesselj(ord,kr);
% i_c = 1i^ord;
% legMat = P_ca1(ord+2,:,:) - P_ca2(ord+2,:,:) - P_ca1(ord+1,:,:) + P_ca2(ord+1,:,:);
% E_f = E_f + i_c.*j_kr.*squeeze(P_ct(ord+1,:,:)).*squeeze(legMat);
% imagesc(abs((E_f)))
% 
% ord=1;
% j_kr = sphbesselj(ord,kr);
% i_c = 1i^ord;
% legMat = P_ca1(ord+2,:,:) - P_ca2(ord+2,:,:) - P_ca1(ord,:,:) + P_ca2(ord,:,:);
% E_f = E_f + i_c.*j_kr.*squeeze(P_ct(ord+1,:,:)).*squeeze(legMat);
% imagesc(abs((E_f)))

for ord=0:N_l
    j_kr = sphbesselj(ord,kr);
    i_c = 1i^ord;
    legMat = P_ca1(ord+2,:,:) - P_ca2(ord+2,:,:) - P_ca1(ord+1,:,:) + P_ca2(ord+1,:,:);
    E_f = E_f + i_c.*j_kr.*squeeze(P_ct(ord+1,:,:)).*squeeze(legMat); 
    imagesc(abs((E_f)))
end

E_0=1;
E_f = 2*pi*E_0*E_f;
imagesc(abs((E_f)))
%imagesc((abs(fftshift(E_f))))



