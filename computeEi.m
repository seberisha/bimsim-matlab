function E_i = computeEi(A, j_knr, P_ct, N_l, e_ikc)
%j(vnrho)P(costheta) for all orders
jknr_Pcosth = j_knr.*P_ct;

%scale by A coefficients for all orders
for j=0:N_l
    jknr_Pcosth(:,:,j+1)=jknr_Pcosth(:,:,j+1)*A(j+1);
end
%sum across all orders
E_i = sum(jknr_Pcosth,3);%phase shift
E_i= e_ikc.*E_i;

