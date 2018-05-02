function E_f = computeEf_vec(j_kr, P_ct, P_ca1, P_ca2, E_0, N_l)

allOrders=j_kr.*P_ct;

ord=0;
allOrders(:,ord+1) = 1i^ord.*allOrders(:,ord+1).*(P_ca1(ord+2)-P_ca2(ord+2)-P_ca1(1)+P_ca2(1));

ord=1;
allOrders(:,ord+1) = 1i^ord.*allOrders(:,ord+1).*(P_ca1(ord+2)-P_ca2(ord+2)-P_ca1(1)+P_ca2(1));

for ord=2:N_l
    allOrders(:,ord+1) = 1i^ord.*allOrders(:,ord+1).*(P_ca1(ord+2)-P_ca2(ord+2)-P_ca1(ord)+P_ca2(ord));
end

E_f = (sum(allOrders,2));
E_f = 2*pi*E_0*E_f;