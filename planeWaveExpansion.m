% this is needed for the plane wave expansion



normRvecs = bsxfun(@rdivide, rVecs,  sqrt(sum(rVecs.^2,2)));
normKvec = params.kVec./params.wavNum;
cos_theta = (normKvec*normRvecs')';
cos_theta = reshape(cos_theta, params.spatRes, params.spatRes);

orderEf=200;
%calculate the prefix term (2l + 1)*i^l
ordVecEf=(0:orderEf);
twolp1 = 2.*ordVecEf+1;
il = 1i.^ordVecEf;
twolp1_il = twolp1.*il;

jl_kr = sphbesselj(orderEf,params.wavNum.*params.r,'multiple');
Pl_costheta = myLegendre(orderEf,cos_theta);
jlkr_Pcostheta = jl_kr.*Pl_costheta;

il_jlkr_Pcostheta = (bsxfun(@times, jlkr_Pcostheta, reshape(twolp1_il,[1 1 orderEf+1])));
params.NA_in = 0; params.NA_out = 1;
alpha1 = asin(params.NA_in); alpha2 = asin(params.NA_out);
Pl_cosalpha1 = myLegendre(orderEf+1,cos(alpha1));
Pl_cosalpha2 = myLegendre(orderEf+1,cos(alpha2));

ord=0;
il_jlkr_Pcostheta(:,:,ord+1) = il_jlkr_Pcostheta(:,:,ord+1).*(Pl_cosalpha1(ord+2)-Pl_cosalpha2(ord+2)-Pl_cosalpha1(1)+Pl_cosalpha2(1));

ord=1;
il_jlkr_Pcostheta(:,:,ord+1) = il_jlkr_Pcostheta(:,:,ord+1).*(Pl_cosalpha1(ord+2)-Pl_cosalpha2(ord+2)-Pl_cosalpha1(1)+Pl_cosalpha2(1));

for ord=2:orderEf
    il_jlkr_Pcostheta(:,:,ord+1) = il_jlkr_Pcostheta(:,:,ord+1).*(Pl_cosalpha1(ord+2)-Pl_cosalpha2(ord+2)-Pl_cosalpha1(ord)+Pl_cosalpha2(ord));
end

E_f = 2*pi*params.E0.*sum(il_jlkr_Pcostheta,3);
brewer = brewermap(1000);
figure, imagesc((abs((E_f)))),title('abs(E_f)'), colorbar, axis image
colormap(brewer)