function E_f = gpuComputeEf_sym(params)
%% focused field
%apply the partial wave expansion to get the focused field

cos_theta = (params.normKvec*params.normPMinPf')';
%cos_theta = reshape(cos_theta, params.simRes, params.simRes);
jl_kr = (sphbesselj(params.orderEf,params.wavNum.*params.r,'multiple'));
jl_kr=gpuArray(jl_kr);
Pl_costheta = gpuLegendre(params.orderEf,cos_theta,params.P);
jlkr_Pcostheta = jl_kr.*Pl_costheta;

params.il_jlkr_Pcostheta = (bsxfun(@times, jlkr_Pcostheta, params.il'));

ord=0;
params.il_jlkr_Pcostheta(:,ord+1) = params.il_jlkr_Pcostheta(:,ord+1).*(params.Pl_cosalpha1(ord+2)-params.Pl_cosalpha2(ord+2)-params.Pl_cosalpha1(1)+params.Pl_cosalpha2(1));

ord=1;
params.il_jlkr_Pcostheta(:,ord+1) = params.il_jlkr_Pcostheta(:,ord+1).*(params.Pl_cosalpha1(ord+2)-params.Pl_cosalpha2(ord+2)-params.Pl_cosalpha1(1)+params.Pl_cosalpha2(1));

for ord=2:params.orderEf
    params.il_jlkr_Pcostheta(:,ord+1) = params.il_jlkr_Pcostheta(:,ord+1).*(params.Pl_cosalpha1(ord+2)-params.Pl_cosalpha2(ord+2)-params.Pl_cosalpha1(ord)+params.Pl_cosalpha2(ord));
end

E_f = 2*pi*params.E0.*sum(params.il_jlkr_Pcostheta,2);
