function E_f = mcEf(j_kr,P_costh,order)

jkr_Pcosth = j_kr.*P_costh;
    
for j=0:order
    jkr_Pcosth(:,:,j+1) = jkr_Pcosth(:,:,j+1)*(2*j+1)*1i^j;
end

E_f = sum(jkr_Pcosth,3);
