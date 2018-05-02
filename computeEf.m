function E_f = computeEf(j_kr, P_ct, P_ca1, P_ca2, E_0, N_l)
%j(v*rho)*P(cos(theta)) for all orders

jkr_Pcosth=j_kr.*P_ct;

ord=0;
jkr_Pcosth(:,:,ord+1) = 1i^ord.*jkr_Pcosth(:,:,ord+1).*(P_ca1(ord+2)-P_ca2(ord+2)-P_ca1(1)+P_ca2(1));

ord=1;
jkr_Pcosth(:,:,ord+1) = 1i^ord.*jkr_Pcosth(:,:,ord+1).*(P_ca1(ord+2)-P_ca2(ord+2)-P_ca1(1)+P_ca2(1));

for ord=2:N_l
    jkr_Pcosth(:,:,ord+1) = 1i^ord.*jkr_Pcosth(:,:,ord+1).*(P_ca1(ord+2)-P_ca2(ord+2)-P_ca1(ord)+P_ca2(ord));
end

E_f = (sum(jkr_Pcosth,3));
E_f = 2*pi*E_0*E_f;