function E_i = computeEi_vec(A, j_knr, P_ct, e_ikc, N_l)
jP = j_knr.*P_ct;

for ord=0:N_l
    jP(:,ord+1) = jP(:,ord+1).*A(ord+1);
end

E_i = (sum(jP,2));
E_i= e_ikc*E_i;

