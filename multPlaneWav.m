%% setup
% single plane wave case
addpath(genpath('~/source/stim-matlab/'))
params.numOrd = 10;

params.n=1.4;

params.gridSize = 5; params.spatRes = 256; 
params.lam = 1; 
params.wavNum = 2*pi/params.lam; %wavenumber
gridPoints = linspace(-params.gridSize,params.gridSize,params.spatRes);
[x,y] = meshgrid(gridPoints, gridPoints);
z = zeros(params.spatRes,params.spatRes);
[theta, phi, r] = cart2sph(x,y,z);
params.r=r;
params.rad = 1;
params.amp = 1;
params.kVec=[1 0 0];
params.kVec = params.kVec.*params.wavNum;

rVecs = zeros(params.spatRes*params.spatRes, 3);
%rVecs(:,1) = theta(:); rVecs(:,2) = phi(:); rVecs(:,3) = r(:);
rVecs(:,1) = x(:); rVecs(:,2) = y(:); rVecs(:,3) = z(:);

ktr = (params.kVec*rVecs')';
ktr = reshape(ktr, params.spatRes, params.spatRes);

params.ordVec = 0:params.numOrd;
ordVec = (0:params.numOrd)';

params.E0 = 1;
params.NA_in = 0;
params.NA_out = 1;
alpha1 = asin(params.NA_in);
alpha2 = asin(params.NA_out);
%% focused field
%apply the partial wave expansion to get the focused field
orderEf = 200;
ordVec_Ef = (0:orderEf)';
il_Ef = 1i.^ordVec_Ef;

jl_kr = sphbesselj(orderEf,params.wavNum.*params.r,'multiple');

%normalize rVecs
normRvecs = bsxfun(@rdivide, rVecs,  sqrt(sum(rVecs.^2,2)));
normKvec = params.kVec/norm(params.kVec);
cos_theta = (normKvec*normRvecs')';
cos_theta = reshape(cos_theta, params.spatRes, params.spatRes);
Pl_costheta_ef = myLegendre(orderEf,cos_theta);
jlkr_Pcostheta = jl_kr.*Pl_costheta_ef;
il_jlkr_Pcostheta = bsxfun(@times, jlkr_Pcostheta, reshape(il_Ef,[1 1 params.numOrd+1]));

Pl_cosalpha1 = myLegendre(orderEf+1,cos(alpha1));
Pl_cosalpha2 = myLegendre(orderEf+1,cos(alpha2));

ord=0;
il_jlkr_Pcostheta(:,:,ord+1) = 1i^ord.*il_jlkr_Pcostheta(:,:,ord+1).*(P_ca1(ord+2)-P_ca2(ord+2)-P_ca1(1)+P_ca2(1));

ord=1;
il_jlkr_Pcostheta(:,:,ord+1) = 1i^ord.*il_jlkr_Pcostheta(:,:,ord+1).*(P_ca1(ord+2)-P_ca2(ord+2)-P_ca1(1)+P_ca2(1));

for ord=2:orderEf
    il_jlkr_Pcostheta(:,:,ord+1) = 1i^ord.*il_jlkr_Pcostheta(:,:,ord+1).*(P_ca1(ord+2)-P_ca2(ord+2)-P_ca1(ord)+P_ca2(ord));
end

E_f = 2*pi*params.E0.*sum(il_jlkr_Pcostheta,3);

brewer = brewermap(1000);
figure, imagesc((abs((E_f)))),title('E_f'), colorbar, axis image
colormap(brewer)

%%


%The scattered field for a single incident plane-wave k produced by
%a sphere with radius r positioned at point pf
samples = 200;

twolp1 = 2.*ordVec+1;
il = 1i.^ordVec;
twolp1_il = twolp1.*il;

params.a=1;
ka=params.wavNum*params.a;
kr = params.wavNum.*r;
kna = params.wavNum*params.n*params.a;

jl_ka = sphbesselj(params.numOrd,ka,'multiple');
jl_kna_p = derivSphBes(params.numOrd, kna);
jl_kna = sphbesselj(params.numOrd,kna,'multiple');
jl_ka_p = derivSphBes(params.numOrd, ka);

numB = jl_ka.*jl_kna_p.*params.n - jl_kna.*jl_ka_p;
hl_ka_p = derivSphHan(params.numOrd, ka);
hl_ka = shank1(params.numOrd, ka, 'multiple');

denAB = jl_kna.*hl_ka_p - hl_ka.*jl_kna_p*params.n;


B = twolp1_il.*(numB./denAB);
hl_kr = shank1(params.numOrd, kr, 'multiple');
% for i=0:params.numOrd
%    Pl_costheta(:,:,i+1)=legendreP(i,cos_theta); 
% end
hlkr_Plcostheta = hl_kr.*Pl_costheta;

k_j = monteCarlo(samples,k_vec, NA_in, NA_out);

E_s = zeros(res,res);
E_i = E_s;
h = waitbar(0, 'Monte Carlo integration...');
%compute the amplitude that makes it through the condenser
amplitude=1;
subA = 2 * pi * amplitude * ( (1 - cos(alpha2)) - (1 - cos(alpha1)) );
if samples==1
    subA=1;
end
E_f=E_i;

for i=1:samples
    E_s = E_s + sum(bsxfun(@times, hlkr_Plcostheta, reshape(B,[1 1 params.numOrd+1])),3);

end


E_s(r<params.a) = 0;%E_i(r<a);
figure, imagesc((abs((E_s)))),title('E_s'), colorbar, axis image
colormap(brewer)

%%
%The internal field specifying the field inside of a sphere
%for an incident plane wave 
knr = params.wavNum*params.n.*r;
jl_knr = sphbesselj(params.numOrd,knr,'multiple');
jlknr_Plcostheta = jl_knr.*Pl_costheta;
numA = jl_ka.*hl_ka_p.*params.n - jl_ka_p.*hl_ka;
A = twolp1_il.*(numA./denAB);

E_i = bsxfun(@times, jlknr_Plcostheta, reshape(A,[1 1 params.numOrd+1]));

E_i = sum(E_i,3);
E_i(r>params.a) = 0;
figure, imagesc((abs((E_i)))),title('E_i'), colorbar, axis image
colormap(brewer)

E_t = zeros(params.spatRes, params.spatRes);
E_t(r<params.a) = E_i(r<params.a);
E_t(r>=params.a) = E_f(r>=params.a) + E_s(r>=params.a);
figure, imagesc((abs((E_t)))),title('E_t'), colorbar, axis image
colormap(brewer)