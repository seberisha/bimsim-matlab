function E_f = cpuComputeEf(params)

%position vectors with respect to the focal point
rVecs = bsxfun(@minus, params.origRVecs,params.pf);
%rVecs(:,2,:) = params.a;
%normalized position vectors with respect to the focal point
normPMinparams.pf = bsxfun(@rdivide, rVecs,  sqrt(sum(rVecs.^2,2)));

params.r=reshape(sqrt(sum(rVecs.^2,2)), params.simRes, params.simRes); %r value at each pixel position

%normalize light direction vector
normparams.kVec = params.kVec./params.magKVector;

cos_theta = (normparams.kVec*normPMinparams.pf')';
cos_theta = reshape(cos_theta, params.simRes, params.simRes);
jl_kr = (sphbesselj(params.orderEf,params.magKVector.*params.r,'multiple'));
jl_kr=(jl_kr);
Pl_costheta = cpuLegendre(params.orderEf,cos_theta);
jlkr_Pcostheta = jl_kr.*Pl_costheta;

params.il_jlkr_Pcostheta = (bsxfun(@times, jlkr_Pcostheta, params.il));

ord=0;
params.il_jlkr_Pcostheta(:,:,ord+1) = params.il_jlkr_Pcostheta(:,:,ord+1).*(params.Pl_cosalpha1(ord+2)-params.Pl_cosalpha2(ord+2)-params.Pl_cosalpha1(1)+params.Pl_cosalpha2(1));

ord=1;
params.il_jlkr_Pcostheta(:,:,ord+1) = params.il_jlkr_Pcostheta(:,:,ord+1).*(params.Pl_cosalpha1(ord+2)-params.Pl_cosalpha2(ord+2)-params.Pl_cosalpha1(1)+params.Pl_cosalpha2(1));

for ord=2:params.orderEf
    params.il_jlkr_Pcostheta(:,:,ord+1) = params.il_jlkr_Pcostheta(:,:,ord+1).*(params.Pl_cosalpha1(ord+2)-params.Pl_cosalpha2(ord+2)-params.Pl_cosalpha1(ord)+params.Pl_cosalpha2(ord));
end

E_f = 2*pi*params.E0.*sum(params.il_jlkr_Pcostheta,3);



