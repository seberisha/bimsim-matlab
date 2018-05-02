function E_f = Ef_singlePlaneWave(j_kr, P_ct, E_0, N_l)

jP=j_kr.*P_ct;

for ord=0:N_l
    jP(:,:,ord+1) = (2*ord+1)*1i^ord.*jP(:,:,ord+1);
end

E_f = (sum(jP,3));
E_f = E_0*E_f;