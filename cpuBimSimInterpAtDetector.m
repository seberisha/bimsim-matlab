function output = cpuBimSimInterpAtDetector(params )
brewer = brewermap(1000);
output = params;

%% set up params
%specify the spatial resolution of the field plane
params.simRes = params.res*(2*params.padding + 1);

% compute alpha1 and alpha2 from NA_in and NA_out, respectively
alpha1 = asin(params.NA_in); alpha2 = asin(params.NA_out);
%calculate the prefix term (2l + 1)*i^l
ordVecEf=((0:params.orderEf)');
params.il = 1i.^ordVecEf;
params.il = reshape(params.il,[1 1 params.orderEf+1]);

%compute the crop size needed for resampling at the detector
cropSize = params.padding*params.res;

if params.padding==0
    params.startIdx=1;
    params.endIdx=params.simRes;
else
    params.startIdx = fix(params.simRes/2 + 1) - floor(cropSize/2);
    params.endIdx = params.startIdx + cropSize-1;
end

%direction of the incident light
[x,y,z] = sph2cart(params.theta,params.phi,1);
lightDirection  = [x y z];

%% matlab coordinates
% get r and rVecs
%create a grid of points representing pixel positions in the field plane
%gridPoints = linspace(-params.gridSize,params.gridSize,params.params.simRes)

halfGridSize = round(params.fov/2)*(2*params.padding + 1);

gx = linspace(-fix(halfGridSize),ceil(halfGridSize)-1,params.simRes);
gy = gx;


% gx = -fix(params.simRes/2):ceil(params.simRes/2)-1;
% gy = -fix(params.simRes/2):ceil(params.simRes/2)-1;

%generate grid points
%gridPoints = (2*gridSize)*(0:params.simRes-1)/params.simRes - gridSize;

[x,z] = meshgrid(gx, gy); % field slice in the x z plane
y = ones(params.simRes,params.simRes)*params.a;   %field plane y = 0
% dgx = gx(params.startIdx:params.endIdx);
% dgy = gy(params.startIdx:params.endIdx);
% [dx, dy] = (meshgrid(dgx(round(params.res/2 + 1):end),dgy(round(params.res/2 + 1):end)));
%
% [theta_p, r_p] = meshgrid(linspace(0,pi/4,16), 0:15);

%% vectors corresponding to each point in the plane
rVecs = zeros(params.simRes*params.simRes, 3);
rVecs(:,1) = x(:); rVecs(:,2) = y(:); rVecs(:,3) = z(:); %r value at each pixel position

%vectors of each point with respect to the distance from center of the
%sphere
rVecs_ps = bsxfun(@minus, rVecs, params.ps);
%mask used later for fparams.iltering the internal and external fields
params.psMask=reshape(sqrt(sum(rVecs_ps.^2,2)),params.simRes, params.simRes);

%legendre polynomials
params.Pl_cosalpha1 = cpuLegendre(params.orderEf+1,cos(alpha1));
params.Pl_cosalpha2 = cpuLegendre(params.orderEf+1,cos(alpha2));

%remember the points (in vector and distance representation) with respect
%to the center of the sphere
params.origRVecs = rVecs;
%normalized vectors
params.normPMinPs = bsxfun(@rdivide, rVecs_ps,  sqrt(sum(rVecs_ps.^2,2)));
%r value at each pixel position with respect to the center of the sphere
params.r_ps=reshape(sqrt(sum(rVecs_ps.^2,2)),params.simRes, params.simRes);

D_Et = zeros(params.res, params.res);
D_Ef = D_Et;
output.allEtd = zeros(params.res, params.res,params.numWav);
output.allEfd = output.allEtd ;

output.A = (zeros(params.res,params.res,params.numWav));
output.absSpec = (zeros(params.numWav,1));

%%

numRad = numel(params.rad);

%%
E_t = zeros(params.simRes,params.simRes);
%h = waitbar(0, 'Per wavelength computation...');

params.lambdas = params.material(:,1);
params.magKVectors = 2*pi./params.lambdas;
params.kVecs = lightDirection.*params.magKVectors;

%% compute for all wavelengths
for materialIdx=1:params.numWav
    %timing for the first wavelength
    if materialIdx==1
        tic
    end
    
    [params.A, params.B] = cpuScatteringCoefficients(params.material(materialIdx,:),params.a);
    
    
    params.BPF = cpuBPF(halfGridSize, params.simRes, params.NA_in/params.lambdas(materialIdx),  params.NA_out/params.lambdas(materialIdx));
    
    
    %monte carlo samples
    params.k_j = (monteCarlo(params.s,params.samples, params.kVecs, params.NA_in, params.NA_out));
    
    D_Et = arrayfun(@(x) 0, D_Et);
    D_Ef = arrayfun(@(x) 0, D_Ef);
    E_t = arrayfun(@(x) 0, E_t);
    
    %% compute E_f and E_t for pf = [0 0 0]
    params.pf=[0 0 0];
    params.subA = 2 * pi * params.E0 * ( (1 - cos(alpha2)) - (1 - cos(alpha1)) );
    
    E_f = cpuComputeEf(params, params.kVecs(materialIdx,:));
    %----do this if performing fft later
    E_f(isnan(E_f))=0; E_f(isinf(E_f))=0;
    %compute the amplitude that makes it through the condenser
    
    [E_s, E_i] = cpuComputeScatteredFields(params, params.material(materialIdx,:));
    E_t(params.psMask<params.a) = E_i(params.psMask<params.a);
    E_t(params.psMask>=params.a) = E_f(params.psMask>=params.a) + E_s(params.psMask>=params.a);
    E_t(isnan(E_t))=0; E_t(isinf(E_t))=0;
    
    [Et_params.BPF, Ef_params.BPF] = applyBandPass(E_t, E_f, params.BPF);
    
    [Et_ff, Ef_ff] = getFarField(Et_params.BPF, Ef_params.BPF,params.startIdx, params.endIdx);
    
    
    %numPoints = round(2*pi*params.rad(end));
    %params.scale = 1/numPoints;
    params.scale=1;
    %integrate and resample
    D_Et = D_Et + (params.scale)*(abs(Et_ff)).^2;   %this is out_i in cuda bimsim -- the measured intesity
    D_Ef = D_Ef + (params.scale)*(abs(Ef_ff)).^2;   %this out_inc in cuda bimsim -- the measured incident field
    
    figure
    output.R_Ef = zeros(1,numRad+1);
    output.R_Ef(1) = D_Ef(round(params.res/2+1), round(params.res/2+1));
    
    output.R_Et = zeros(1,numRad+1);
    output.R_Et(1) = D_Et(round(params.res/2+1), round(params.res/2+1));
    
    
    %     dx = -fix(params.res/2):ceil(params.res/2)-1;
    %     dy = -fix(params.res/2):ceil(params.res/2)-1;
    
    
    %[D_x,D_y] = meshgrid(dx,dy);
    D_x = x(params.startIdx:params.endIdx,params.startIdx:params.endIdx);
    D_y = z(params.startIdx:params.endIdx,params.startIdx:params.endIdx);
    [~, D_r] = cart2pol(D_x,D_y);
    
    rho = params.rad(end);
    numPoints = round(((2*pi*rho-params.spacing)/params.spacing + 1));
    maxRange = floor(2*pi*rho-params.spacing/2);
    s = linspace(0,maxRange,numPoints);
    params.pf_theta = s./rho;
    %          figure
    %          scatter(0,0,[],[1 1 0],'filled')
    %          hold on
    
    if params.numPS>1
        radiusIdx=1;
        for p=params.spacing:params.spacing:numRad
            disp(['radius = ' num2str(p)])
            temp = -log10(D_Et./D_Ef);
            dispvar('absorbance at center',temp(round(params.res/2+1), round(params.res/2+1)))
            %disp(['absorbance at center = ' num2str(temp(params.res/2+1,params.res/2+1))])
            subplot(221)
            imagesc(temp), axis image, colorbar, colormap(brewer)
            drawnow
            
            rho=params.rad(radiusIdx);
            
            %numPoints = round(((2*pi*rho-params.spacing)/params.spacing + 1));
            %numPoints = round(2*pi*rho);
            
            dispvar('num point sources',numPoints)
            
            %maxRange = floor(2*pi*rho-params.spacing/2);
            %s = linspace(0,maxRange,numPoints);
            %s = 0:numPoints;
            %params.pf_theta = s./rho;
            [params.pf_x,params.pf_z] = pol2cart(params.pf_theta,rho);
            
            %scatter(params.pf_x,params.pf_z,[],[1 1 0]* params.amplitude(p),'filled')
            
            % scale point sources by a Gaussian amplitude and average them
            %             % in each ring
            params.scale = params.amplitude(p)/numPoints;
            dispvar('gauss amplitude',params.amplitude(p))
            
            dispvar('scale', params.scale)
            dispvar('field amplitude',params.E0)
            params.subA = 2 * pi * params.E0 * ( (1 - cos(alpha2)) - (1 - cos(alpha1)) );
            
            
            
            [D_Et, D_Ef,output] = cpuComputePointSourcesInterpAtDetector(D_x, D_y, p, materialIdx,D_Et, D_Ef, E_t, params, output);
            
            D_Ef = interp1([0 params.rad(1:p)], output.R_Ef(1:p+1),D_r);
            D_Et = interp1([0 params.rad(1:p)], output.R_Et(1:p+1), D_r);
            
            subplot(222)
            imagesc(D_Ef), axis image, colorbar, colormap(brewer)
            subplot(223)
            imagesc(D_Et), axis image, colorbar, colormap(brewer)
            drawnow
            dispvar()
            radiusIdx = radiusIdx + 1;
        end
        
        %axis square
        %set(gca,'color','black')
        
        D_Ef_t = interp1([0 params.rad], output.R_Ef,D_r);
        D_Ef_t(isnan(D_Ef_t))=0; D_Ef_t(isinf(D_Ef_t))=0;
        D_Et_t = interp1([0 params.rad], output.R_Et, D_r);
        D_Et_t(isnan(D_Et_t))=0; D_Et_t(isinf(D_Et_t))=0;
        
    end
    
    
    %save results to a structure
    
    output.A(:,:,materialIdx) = -log10((D_Et_t)./(D_Ef_t));
    output.absSpec(materialIdx) = -log10(sum(D_Et_t(:))/sum(D_Ef_t(:)));
    dispvar('final abs',output.absSpec(materialIdx))
    dispvar()
    
    output.allEtd(:,:,materialIdx) = D_Et;
    output.allEfd(:,:,materialIdx) = D_Ef;
    
    if materialIdx==1
        toc
    end
    
end