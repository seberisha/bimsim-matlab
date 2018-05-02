function E_f = computeEf_sym(params)
%% focused field
%apply the partial wave expansion to get the focused field

normrVecs = bsxfun(@rdivide, params.rVecs,  params.r);
normKvec = params.kVec./params.wavNum;
cos_theta = (normKvec*normrVecs')';

%calculate the prefix term (2l + 1)*i^l
ordVecEf=(0:params.orderEf)';
il = 1i.^ordVecEf;

jl_kr = sphbesselj(params.orderEf,params.wavNum.*params.r,'multiple');
Pl_costheta = myLegendre(params.orderEf,cos_theta);
jlkr_Pcostheta = jl_kr.*Pl_costheta;

il_jlkr_Pcostheta = (bsxfun(@times, jlkr_Pcostheta, il'));

ord=0;
il_jlkr_Pcostheta(:,ord+1) = il_jlkr_Pcostheta(:,ord+1).*(params.Pl_cosalpha1(ord+2)-params.Pl_cosalpha2(ord+2)-params.Pl_cosalpha1(1)+params.Pl_cosalpha2(1));

ord=1;
il_jlkr_Pcostheta(:,ord+1) = il_jlkr_Pcostheta(:,ord+1).*(params.Pl_cosalpha1(ord+2)-params.Pl_cosalpha2(ord+2)-params.Pl_cosalpha1(1)+params.Pl_cosalpha2(1));

for ord=2:params.orderEf
    il_jlkr_Pcostheta(:,ord+1) = il_jlkr_Pcostheta(:,ord+1).*(params.Pl_cosalpha1(ord+2)-params.Pl_cosalpha2(ord+2)-params.Pl_cosalpha1(ord)+params.Pl_cosalpha2(ord));
end

E_f = 2*pi*params.E0.*sum(il_jlkr_Pcostheta,2);
