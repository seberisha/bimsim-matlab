function [E_f, coords_r, coords_c] = focusedField(E_0, lambda, NA_c, alpha_1, res)
%Description:
%

%Assumptions:
%   -Incident light is coherent: |k_j|=k, for all j \in J.

%Input: 
%   E_0 - a vector providing the amplitude and polarization direction
%               of the field
%   lambda - wavelength
%   NA_c - condenser aperture
%   k - magnitudes of vectors of the propogating plane wave
%   a - sphere radius
%   alpha_1 - the inner angle

%Authors: S.Berisha
%Last modified: 09/16/15

% %%
% [theta_p,r_p] = meshgrid((linspace(0,360,res))*pi/180,linspace(0,1,res));
% 
% %[theta,r] = pol2cart(theta_p,r_p);
% %% the maximum order required for convergence
% N_l = maxOrder(r_p, lambda);
% 
% %% precompute spherical bessel functions
% 
% lOrders = 0:N_l(1,1); 
% k=2*pi/lambda;
% kr = k*r_p(1,1);
% j_kr = sphbesselj(lOrders,kr);
% 
% %% precompute Legendre polynomials
% 
% P_ct = legendre(N_l(1,1),cos(theta_p(1,1)));
% %for lenses
% alpha_2 = 1/sin(NA_c);
% P_ca1 = legendre(N_l(1,1)+1, cos(alpha_1));
% P_ca2= legendre(N_l(1,1)+1, cos(alpha_2));
% %P_ca1_lp1 = legendre(N_l(1,1)+1, cos(alpha_1));
% % P_ca2_lp1= legendre(nu+1, cos(alpha_2));
% % P_ca1_lm1 = legendre(nu-1, cos(alpha_1));
% % P_ca2_lm1 = legendre(nu-1, cos(alpha_2));
% 
% %% compute E_f
% legVec = zeros(length(lOrders),1);
% legVec(1)=P_ca1(2) - P_ca2(2);
% for l=2:N_l(1,1)+1
%     legVec(l)=P_ca1(l+1) - P_ca2(l+1) - P_ca1(l-1) + P_ca2(l-1);
% end
%   
% E_f = zeros(res,res);
% for row=1:1
%     for col=1:1
%         i_vec = i.^lOrders;
%         E_f(row,col) = 2*pi*E_0 *(i_vec.*j_kr)*legVec;
%     end
% end

%%
[x,y] = meshgrid(linspace(-1,1,res), linspace(-1,1,res));
[theta, r] = cart2pol(x,y);
%%
N_l = maxOrder(abs(r), lambda);

k=2*pi/lambda;
alpha_1=0;
alpha_2=.76;
%alpha_2 = 1/sin(NA_c);
E_f = zeros(res,res);
temp = round(res/2);
coords_r=zeros(res,res);
coords_c=zeros(res,res);
for row=1:res
    for col=1:res
        lOrders = 0:N_l(row,col); 
        kr = k*abs(r(row,col));
        j_kr = sphbesselj(lOrders,kr);

        P_ct = legendre(N_l(row,col),cos(theta(row,col)));
        P_ca1 = legendre(N_l(row,col)+1, cos(alpha_1));
        P_ca2= legendre(N_l(row,col)+1, cos(alpha_2));
        legVec = zeros(length(lOrders),1);
        legVec = P_ca1(2:end) - P_ca2(2:end) - [P_ca1(1);P_ca1(1:end-2)] + [P_ca2(1);P_ca2(1:end-2)];
        
        
%         legVec(1)=P_ca1(2) - P_ca2(2);
%         for l=2:N_l(row,col)+1
%             legVec(l)=P_ca1(l+1) - P_ca2(l+1) - P_ca1(l-1) + P_ca2(l-1);
%         end
        i_vec = i.^lOrders;
        
        [x,y]=pol2cart(theta(row,col),r(row,col));
      
        m_r = ceil(temp - y*temp);
        m_c = ceil(temp + x*temp);
        if(m_r==0),m_r=1;end
        if(m_c==0),m_c=1;end
        coords_r(row,col)=m_r;
        coords_c(row,col)=m_c;
        E_f(m_r,m_c) = 2*pi*E_0 *(i_vec.*j_kr.*P_ct')*legVec;
        
    end
end


