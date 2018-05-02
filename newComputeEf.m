function E_f = newComputeEf(params)
%% focused field
%apply the partial wave expansion to get the focused field

normrVecs = bsxfun(@rdivide, params.rVecs,  sqrt(sum(params.rVecs.^2,2)));
normKvec = params.kVec./params.wavNum;
cos_theta = (normKvec*normrVecs')';
cos_theta = reshape(cos_theta, params.simRes, params.simRes);

%calculate the prefix term (2l + 1)*i^l
ordVecEf=(0:params.orderEf)';
il = 1i.^ordVecEf;

jl_kr = sphbesselj(params.orderEf,params.wavNum.*params.r,'multiple');
Pl_costheta = myLegendre(params.orderEf,cos_theta);
jlkr_Pcostheta = jl_kr.*Pl_costheta;

il_jlkr_Pcostheta = (bsxfun(@times, jlkr_Pcostheta, reshape(il,[1 1 params.orderEf+1])));

ord=0;
il_jlkr_Pcostheta(:,:,ord+1) = il_jlkr_Pcostheta(:,:,ord+1).*(params.Pl_cosalpha1(ord+2)-params.Pl_cosalpha2(ord+2)-params.Pl_cosalpha1(1)+params.Pl_cosalpha2(1));

ord=1;
il_jlkr_Pcostheta(:,:,ord+1) = il_jlkr_Pcostheta(:,:,ord+1).*(params.Pl_cosalpha1(ord+2)-params.Pl_cosalpha2(ord+2)-params.Pl_cosalpha1(1)+params.Pl_cosalpha2(1));

for ord=2:params.orderEf
    il_jlkr_Pcostheta(:,:,ord+1) = il_jlkr_Pcostheta(:,:,ord+1).*(params.Pl_cosalpha1(ord+2)-params.Pl_cosalpha2(ord+2)-params.Pl_cosalpha1(ord)+params.Pl_cosalpha2(ord));
end

E_f = 2*pi*params.E0.*sum(il_jlkr_Pcostheta,3);

save 2d_computeEf.mat
