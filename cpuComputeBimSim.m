function output = cpuComputeBimSim(params )
brewer = brewermap(1000);
output = params;

%% set up params


%specify the spatial resolution of the field plane
simRes = params.res*(2*params.padding + 1);

% compute alpha1 and alpha2 from NA_in and NA_out, respectively
alpha1 = asin(params.NA_in); alpha2 = asin(params.NA_out);
%calculate the prefix term (2l + 1)*i^l
ordVecEf=((0:params.orderEf)');
il = 1i.^ordVecEf;
il = reshape(il,[1 1 params.orderEf+1]);
%compute the amplitude that makes it through the condenser
subA = 2 * pi * params.E0 * ( (1 - cos(alpha2)) - (1 - cos(alpha1)) );


%compute the crop size needed for resampling at the detector
cropSize = params.padding*params.res;

if params.padding==0
    startIdx=1;
    endIdx=simRes;
else
    startIdx = round(simRes/2) - floor(cropSize/2);
    endIdx = startIdx + cropSize-1;
end

%direction of the incident light
[x,y,z] = sph2cart(params.theta,params.phi,1);
lightDirection  = [x y z];

%%
% matlab coordinates
% get r and rVecs
%create a grid of points representing pixel positions in the field plane
%gridPoints = linspace(-params.gridSize,params.gridSize,params.simRes)

gridSize = round(params.fov/2)*(2*params.padding + 1);

%generate grid points
gridPoints = (2*gridSize)*(0:simRes-1)/simRes - gridSize;

[x,z] = meshgrid(gridPoints, gridPoints); % field slice in the x z plane
y = ones(simRes,simRes)*(params.a);   %field plane y = 0

%%
%vectors corresponding to each point in the plane
rVecs = zeros(simRes*simRes, 3);
rVecs(:,1) = x(:); rVecs(:,2) = y(:); rVecs(:,3) = z(:); %r value at each pixel position

%vectors of each point with respect to the distance from center of the
%sphere
rVecs_ps = bsxfun(@minus, rVecs, params.ps);
%mask used later for filtering the internal and external fields
psMask=reshape(sqrt(sum(rVecs_ps.^2,2)),simRes, simRes);

%legendre polynomials
Pl_cosalpha1 = cpuLegendre(params.orderEf+1,cos(alpha1));
Pl_cosalpha2 = cpuLegendre(params.orderEf+1,cos(alpha2));

%remember the points (in vector and distance representation) with respect
%to the center of the sphere
origRVecs = rVecs;
%normalized vectors
normPMinPs = bsxfun(@rdivide, rVecs_ps,  sqrt(sum(rVecs_ps.^2,2)));
%r value at each pixel position with respect to the center of the sphere
r_ps=reshape(sqrt(sum(rVecs_ps.^2,2)),simRes, simRes);

gpu_Es = zeros(simRes,simRes);
gpu_Ei = gpu_Es;
E_t = gpu_Ei;

D_Et = zeros(params.res, params.res);
D_Ef = D_Et;
output.allEtd = zeros(params.res, params.res,params.numWav);
output.allEfd = output.allEtd ;

output.A = (zeros(params.res,params.res,params.numWav));
output.absSpec = (zeros(params.numWav,1));

%%

numRad = numel(params.rad);
%[pf_x, pf_z, pf_theta] = pointSources(params.rad, params.numPS);

%%
% figure
% for i=1:numRad
%     plot(pf_x(i,:),pf_z(i,:),'y*','linewidth',2);
%     hold on
% end
% xlim([-params.fov/2 params.fov/2])
% ylim([-params.fov/2 params.fov/2])
% axis square
% set(gca,'color','black')
%set(gcf,'color','black')
%whitebg(gcf,'k')
%%
% needed for rotation of coordinate system if using interp2
[sz1,sz2] = size(E_t);
[X,Y] = meshgrid(1:sz2, 1:sz1);
[Xi,Yi] = meshgrid(-(sz2-1)/2:(sz2-1)/2,-(sz1-1)/2:(sz1-1)/2);
% %loose interp2
% sz1New = sz1*cos(pf_theta(p))+sz2*sin(pf_theta(p));
% if (sz1New<0)
%     sz1New = -sz1New;
% end
%
% sz2New = sz2*cos(pf_theta(p))+sz1*sin(pf_theta(p));
%
% if (sz2New<0)
%     sz2New = -sz2New;
% end

%[Xi,Yi] = meshgrid(-(sz2New-1)/2:(sz2New-1)/2,-(sz1New-1)/2:(sz1New-1)/2);
%%
E_t = zeros(simRes,simRes);
h = waitbar(0, 'Per wavelength computation...');

%% compute for all wavelengths
for i=1:params.numWav
    %timing for the first wavelength
    if i==1
        tic
    end
    
    [A, B] = cpuScatteringCoefficients(params.material(i,:),params.a);
    
    lambda = params.material(i,1);
    wavNum = 2*pi/lambda;
    
    BPF = cpuBPF(gridSize, simRes, params.NA_in/lambda,  params.NA_out/lambda);
    
    kVec=(lightDirection*wavNum);
    
    %monte carlo samples
    k_j = (monteCarlo(params.s,params.samples, kVec, params.NA_in, params.NA_out));
    
    D_Et = arrayfun(@(x) 0, D_Et);
    D_Ef = arrayfun(@(x) 0, D_Ef);
    E_t = arrayfun(@(x) 0, E_t);
    
    %compute E_f and E_t for pf = [0 0 0]
    pf=[0 params.a 0];
    E_f = cpuComputeEf(params, kVec, wavNum,simRes,pf,origRVecs,Pl_cosalpha1,Pl_cosalpha2,il);
    [E_s, E_i] = cpuComputeScatteredFields(params,wavNum, r_ps, params.material(i,:), pf, simRes, subA, k_j, normPMinPs, A, B, psMask);
    E_t(psMask<params.a) = E_i(psMask<params.a);
    E_t(psMask>=params.a) = E_f(psMask>=params.a) + E_s(psMask>=params.a);
    E_t(isnan(E_t))=0; E_t(isinf(E_t))=0;
    
    [Et_bpf, Ef_bpf] = applyBandPass(E_t, E_f, BPF);
    
    [Et_ff, Ef_ff] = getFarField(Et_bpf, Ef_bpf,startIdx, endIdx);
        
    %integrate and resample
    D_Et = D_Et + (abs(Et_ff)).^2;   %this is out_i in cuda bimsim -- the measured intesity
    D_Ef = D_Ef + (abs(Ef_ff)).^2;   %this out_inc in cuda bimsim -- the measured incident field
    -log10(D_Et(params.res/2,params.res/2)/D_Ef(params.res/2,params.res/2))
    temp= -log10(D_Et./D_Ef);
    figure, imagesc(((temp))),axis image,colorbar,colormap(brewer)
    drawnow
    output.Abs(:,:,1) = -log10((D_Et)./(D_Ef));
    output.D_Et(:,:,1) = D_Et;
    output.D_Ef(:,:,1) = D_Ef;
    
    if params.numPS>=1
        for p=1:numRad
           
            numPoints = (params.numPS)*p;
            pf_theta=linspace(0,2*pi-0.01,numPoints);
            rho=ones(1,numPoints)*params.rad(p);
            [pf_x,pf_z] = pol2cart(pf_theta,rho);
            params.E0 = params.E0 - 0.02;
            %compute the amplitude that makes it through the condenser
            subA = 2 * pi * params.E0 * ( (1 - cos(alpha2)) - (1 - cos(alpha1)) );
            params.E0
            [D_Et, D_Ef, output] = cpuComputePointSources(output, X,Y,Xi,Yi,params.interpolate,startIdx, endIdx, BPF, E_t, i,params,pf_x,pf_z,pf_theta, D_Et,D_Ef,kVec,wavNum, simRes, origRVecs, Pl_cosalpha1, Pl_cosalpha2,il,r_ps,subA, k_j, normPMinPs, A,B, psMask);
             -log10(D_Et(params.res/2,params.res/2)/D_Ef(params.res/2,params.res/2))
             temp= -log10(D_Et./D_Ef);
             output.Abs(:,:,p+1) = -log10((D_Et)./(D_Ef));
             output.D_Et(:,:,p+1) = D_Et;
             output.D_Ef(:,:,p+1) = D_Ef;
             p
             imagesc(((temp))),axis image,colorbar,colormap(brewer)
             drawnow
        end
    end
    
    waitbar(i/params.numWav, h, sprintf('wavelength %f',params.material(i,1)));
    %calculate absorbance
    %save results to a structure
    
    output.A(:,:,i) = -log10((D_Et)./(D_Ef));
    output.absSpec(i) = -log10(sum(D_Et(:))/sum(D_Ef(:)));
    %output.absSpec(i)
    output.allEtd(:,:,i) = D_Et;
    output.allEfd(:,:,i) = D_Ef;
    
    if i==1
        toc
    end
    
end