function E_f = parComputeEf(rVecs,kVec,wavNum,simRes,orderEf,r,Pl_cosalpha1,Pl_cosalpha2,E0 )
%% focused field
%apply the partial wave expansion to get the focused field

normrVecs = bsxfun(@rdivide, rVecs,  sqrt(sum(rVecs.^2,2)));
normKvec = kVec./wavNum;
cos_theta = (normKvec*normrVecs')';
cos_theta = reshape(cos_theta, simRes, simRes);

%calculate the prefix term (2l + 1)*i^l
ordVecEf=(0:orderEf)';
il = 1i.^ordVecEf;

jl_kr = sphbesselj(orderEf,wavNum.*r,'multiple');
Pl_costheta = myLegendre(orderEf,cos_theta);
jlkr_Pcostheta = jl_kr.*Pl_costheta;

il_jlkr_Pcostheta = (bsxfun(@times, jlkr_Pcostheta, reshape(il,[1 1 orderEf+1])));

ord=0;
il_jlkr_Pcostheta(:,:,ord+1) = il_jlkr_Pcostheta(:,:,ord+1).*(Pl_cosalpha1(ord+2)-Pl_cosalpha2(ord+2)-Pl_cosalpha1(1)+Pl_cosalpha2(1));

ord=1;
il_jlkr_Pcostheta(:,:,ord+1) = il_jlkr_Pcostheta(:,:,ord+1).*(Pl_cosalpha1(ord+2)-Pl_cosalpha2(ord+2)-Pl_cosalpha1(1)+Pl_cosalpha2(1));

slicePl_cosalpha1 = zeros(size(Pl_cosalpha1));
slicePl_cosalpha1(:,:,4:end) = Pl_cosalpha1(:,:,4:end);

slicePl_cosalpha2 = zeros(size(Pl_cosalpha2));
slicePl_cosalpha2(:,:,4:end) = Pl_cosalpha2(:,:,4:end);

parfor ord=2:orderEf
    il_jlkr_Pcostheta(:,:,ord+1) = il_jlkr_Pcostheta(:,:,ord+1).*(slicePl_cosalpha1(ord+2)-slicePl_cosalpha2(ord+2)-Pl_cosalpha1(ord)+Pl_cosalpha2(ord));
end

E_f = 2*pi*E0.*sum(il_jlkr_Pcostheta,3);
