function E_s = computeEs_vec(B, h_kr, P_ct, N_l, e_ikc)
hP = h_kr.*P_ct;

for ord=0:N_l
    hP(:,ord+1) = hP(:,ord+1).*B(ord+1);
end

E_s = (sum(hP,2));
E_s = e_ikc*E_s;

