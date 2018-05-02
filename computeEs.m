function E_s = computeEs(B, h_kr, P_ct, N_l, e_ikc)
%h(v*rho)*P(cos(theta)) for all orders
hkr_Pcosth = h_kr.*P_ct;


%scale by B coefficients for all orders


for j=0:N_l
    hkr_Pcosth(:,:,j+1) = hkr_Pcosth(:,:,j+1)*B(j+1);
end

%add up for all orders
E_s = sum(hkr_Pcosth,3);

%phase shift
E_s = e_ikc.*E_s;

