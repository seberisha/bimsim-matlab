function output = gpuBimSimInterpAtDetector(params )
%% set up params
brewer = brewermap(1000);
output = params; %save input params to output structure

%spatial resolution of the field plane
params.simRes = params.res*(2*params.padding + 1);

% alpha1 and alpha2 from condenserNA_in and condenserNA_out, respectively
alpha1 = asin(params.condenserNA_in); alpha2 = asin(params.condenserNA_out);
ordVecEf=gpuArray((0:params.orderEf)'); %the prefix term (2l + 1)*i^l
params.il = 1i.^ordVecEf;
params.il = reshape(params.il,[1 1 params.orderEf+1]);

%compute the crop size needed for resampling at the detector
cropSize = params.padding*params.res;

%indices for resampling at far field
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

halfGridSize = round(params.fov/2)*(2*params.padding + 1);
gx = linspace(-fix(halfGridSize),ceil(halfGridSize)-1,params.simRes);
gy = gx;

[x,z] = meshgrid(gx, gy); % field slice in the x z plane

%field plane y at the radius of the sphere
y = ones(params.simRes,params.simRes)*params.a;
%% vectors corresponding to each point in the plane
rVecs = zeros(params.simRes*params.simRes, 3);
%r value at each pixel position
rVecs(:,1) = x(:); rVecs(:,2) = y(:); rVecs(:,3) = z(:);

%vectors of each point with respect to the distance from center of the
%sphere
rVecs_ps = bsxfun(@minus, rVecs, params.ps);

%mask used later for separating the internal and external fields in total
%field calculation
params.psMask=reshape(sqrt(sum(rVecs_ps.^2,2)),params.simRes, params.simRes);

%legendre polynomials
params.Pl_cosalpha1 = gpuLegendre(params.orderEf+1,cos(alpha1));
params.Pl_cosalpha2 = gpuLegendre(params.orderEf+1,cos(alpha2));

%remember the points (in vector and distance representation) with respect
%to the center of the sphere
params.origRVecs = rVecs;
%normalized vectors with respect to distance from sphere center
params.normPMinPs = bsxfun(@rdivide, rVecs_ps,  sqrt(sum(rVecs_ps.^2,2)));
%r value at each pixel position with respect to the center of the sphere
params.r_ps=reshape(sqrt(sum(rVecs_ps.^2,2)),params.simRes, params.simRes);

gpu_Es = zeros(params.simRes,params.simRes,'gpuArray'); %scattered field
gpu_Ei = gpu_Es;    %internal field

D_Et = zeros(params.res, params.res,'gpuArray');    %total field at the detector
D_Ef = D_Et;    %focused/incident field at the detector
%save total and focused fields at the detector to output struct
output.allEtd = zeros(params.res, params.res,params.numWav,'gpuArray');
output.allEfd = output.allEtd ;

%save absorbance image to output struct
output.A = gpuArray(zeros(params.res,params.res,params.numWav));
%save absorbance spectrum for all wavelengths to output struct
output.absSpec = gpuArray(zeros(params.numWav,1));

numRad = numel(params.rad); %number of point source samples along a radius

%h = waitbar(0, 'Per wavelength computation...');

params.lambdas = params.material(:,1);  %wavelengths in micrometers
%magnitude of k vectors for each wavelength
params.magKVectors = 2*pi./params.lambdas;
%k vectors for each wavelength
params.kVecs = lightDirection.*params.magKVectors;



    
%% compute for all wavelengths
for materialIdx=1:params.numWav
    %timing for the first wavelength
    if materialIdx==1
        tic
    end
    
    params.E0 = 1;
    params.subA = 2 * pi * params.E0 * ( (1 - cos(alpha2)) - (1 - cos(alpha1)) );
    
    
    %scattering coefficients
    [params.A, params.B] = gpuScatteringCoefficients(params.material(materialIdx,:),params.a);
    
    %band pass filter
    params.BPF = gpuBPF(halfGridSize, params.simRes, params.objectiveNA_in/params.lambdas(materialIdx),  params.objectiveNA_out/params.lambdas(materialIdx));
    
    %monte carlo samples
    params.k_j = gpuArray(monteCarlo(params.s,params.samples, params.kVecs, params.condenserNA_in, params.condenserNA_out));
    
    % compute E_f and E_t for pf = [0 0 0]
    params.pf=[0 0 0];
    
    [E_f, E_t] = gpuComputeNearFields(params, materialIdx);
       
    %angular spectrum and band pass filtering
    [Et_bpf, Ef_bpf] = applyBandPass(E_t, E_f, params.BPF);
    
    %resample to get far field
    [Et_ff, Ef_ff] = getFarField(Et_bpf, Ef_bpf,params.startIdx, params.endIdx);
    
    %no scaling (averaging and Gaussian amplitude) for the first point
    %source at [0 0 0]
    params.scale=1;
    
    % detector x, z coordinates
    x_fp = x - params.pf(1);
    z_fp = z - params.pf(3);
    D_x_fp = x_fp(params.startIdx:params.endIdx,params.startIdx:params.endIdx);
    D_z_fp = z_fp(params.startIdx:params.endIdx,params.startIdx:params.endIdx);
    [~, D_r_fp] = cart2pol(D_x_fp,D_z_fp);   %rho polar values at the above cartesian coordinates
    
    D_x_ps = x(params.startIdx:params.endIdx,params.startIdx:params.endIdx);
    D_z_ps = z(params.startIdx:params.endIdx,params.startIdx:params.endIdx);
    [~, D_r_ps] = cart2pol(D_x_ps,D_z_ps);   %rho polar values at the above cartesian coordinates
    
    
    % focused/incident and total field values at the detector for each ring
    output.R_Ef = zeros(1,numRad+1,'gpuArray');
    output.R_Et = zeros(1,numRad+1,'gpuArray');
    
    %integrate
    D_Et = (params.scale)*(abs(Et_ff)).^2;   %this is out_i in cuda bimsim -- the measured intesity
    D_Ef = (params.scale)*(abs(Ef_ff)).^2;   %this out_inc in cuda bimsim -- the measured incident field
    
    % Ef and Et values at [0 0  0] taken from the above full field computation
    output.R_Ef(1) = D_Ef(round(params.res/2+1), round(params.res/2+1));
    output.R_Et(1) = D_Et(round(params.res/2+1), round(params.res/2+1));
        
    rho = params.rad(end);  %largest radius
    
    % number of point source samples at the outermost radius -- chosen such
    % that spacing between sectors in the circle is 1
    numPoints = ceil(((2*pi*rho-params.spacing)/params.spacing + 1));
    numPoints=1;
    maxRange = floor(2*pi*rho-params.spacing/2);
    s = linspace(0,maxRange,numPoints); %sectors
    
    params.pf_theta = 2*pi/numPoints:2*pi/numPoints:2*pi;
    R = numel(params.rad)+1;
    radius = params.rad(end);   
  
    figure
    % get the 1d profiles for the first point source located at [0 0 0]
    totalEf_profile = (simRing(numPoints, R, radius, D_Ef, D_x_fp, D_z_fp));
    
    %totalEf_profile(isnan(totalEf_profile))=0;
    
    %totalEf_profile = totalEf_profile./max(totalEf_profile(:));
    
    totalEt_profile = (simRing(numPoints, R, radius, D_Et, D_x_ps, D_z_ps));
    
    %totalEt_profile(isnan(totalEt_profile))=0;
   % totalEt_profile = totalEt_profile./max(totalEt_profile(:));
    subplot(121)
    plot(totalEf_profile)
    hold on
    
    subplot(122)
    plot(totalEf_profile); hold on
    
    % obtain and sum the 1d profiles of the rest of the point sources
    for i=1:numRad

        dispvar('ring',i)
        %focal point for the first point source in the ring
        params.pf(1) = params.rad(i);
        params.pf(2) = 0;
        params.pf(3) = 0;
                
        %params.E0 = numPoints;
        params.E0 = 1;
        %params.E0 = 2*pi*params.rad(i);

        params.subA = 2 * pi * params.E0 * ( (1 - cos(alpha2)) - (1 - cos(alpha1)) );
        
        %compute E_f, E_t for the first point source
        [E_f, E_t] = gpuComputeNearFields(params, materialIdx);
           
        %apply bandpass filtering in the frequency domain
        [Et_bpf, Ef_bpf] = applyBandPass(E_t, E_f, params.BPF);
        
        %resample to get far field
        [Et_ff, Ef_ff] = getFarField(Et_bpf, Ef_bpf,params.startIdx, params.endIdx);
        
        %save intensities at the detector for the first point source
        D_Et_psi = (abs(Et_ff)).^2;
        D_Ef_psi = (abs(Ef_ff)).^2;
        
        % detector x, z coordinates
        x_fp = x - params.pf(1);
        z_fp = z - params.pf(3);
        D_x_fp = x_fp(params.startIdx:params.endIdx,params.startIdx:params.endIdx);
        D_z_fp = z_fp(params.startIdx:params.endIdx,params.startIdx:params.endIdx);
        [~, D_r_fp] = cart2pol(D_x_fp,D_z_fp);   %rho polar values at the above cartesian coordinates
        
        rho=params.rad(i);
        numPoints = ceil(((2*pi*rho-params.spacing)/params.spacing + 1));
        params.pf_theta = 2*pi/numPoints:2*pi/numPoints:2*pi;

        temp_Ef = (simRing(numPoints, R, radius, D_Ef_psi, D_x_ps, D_z_ps));
        %temp_Ef = temp_Ef./max(temp_Ef(:));
        %temp_Ef(isnan(temp_Ef))=0;
        subplot(121)
        plot(temp_Ef)
        hold on
        drawnow
        
        totalEf_profile = totalEf_profile + temp_Ef;
        subplot(122)
        plot(totalEf_profile)
        drawnow
        
        temp_Et = (simRing(numPoints, R, radius, D_Et_psi, D_x_ps, D_z_ps));
        
        %temp_Et(isnan(temp_Et))=0;
        %temp_Et = temp_Et./max(temp_Et(:));
        
        
        totalEt_profile = totalEt_profile + temp_Et;
            
  
        dispvar()
        
    end
    
    %scaling?    
    %totalEf_profile(2:end) = totalEf_profile(2:end).*params.amplitude;
    %totalEt_profile(2:end) = totalEt_profile(2:end).*params.amplitude;
    
    D_Ef_t = interp1([0 params.rad], totalEf_profile,D_r_ps);
    D_Ef_t(isnan(D_Ef_t))=0; D_Ef_t(isinf(D_Ef_t))=0;
    D_Et_t = interp1([0 params.rad], totalEt_profile, D_r_ps);
    D_Et_t(isnan(D_Et_t))=0; D_Et_t(isinf(D_Et_t))=0;
    
    
    

    
    
    
    
    %figure
    %%
    %extended source computation
    %     if params.numPS>1
    %         % for all point source rings
    %         for p=1:numRad
    %             temp = -log10(D_Et./D_Ef);
    %             dispvar('radius', p)
    %             dispvar('absorbance at center',temp(round(params.res/2+1), round(params.res/2+1)))
    % %             subplot(221)
    % %             imagesc(temp), axis image, colorbar, colormap(brewer)
    % %             drawnow
    %
    %             rho=params.rad(p);  %radius
    %
    %             dispvar('num point sources',numPoints)
    %
    %             %maxRange = floor(2*pi*rho-params.spacing/2);
    %             %s = linspace(0,maxRange,numPoints);
    %             %s = 0:numPoints;
    %             %params.pf_theta = s./rho;
    %
    %             %x and z coordinates for the point source ring
    %             [params.pf_x,params.pf_z] = pol2cart(params.pf_theta,rho);
    %
    %             % scale point sources by a Gaussian amplitude and average them
    %             % in each ring
    %             params.scale = params.amplitude(p)/numPoints;
    %
    %             dispvar('gauss amplitude',params.amplitude(p))
    %
    %             dispvar('scale', params.scale)
    %             dispvar('field amplitude',params.E0)
    %             dispvar()
    %
    %             % compute the amplitude that makes it through the condenser
    %             %   -for a single plane wave this should be equal to E0
    %             %   -no need to recompute if amplitude doesn't change
    %             %   -when computing the absorbance image amplitude gets
    %             %   cancelled
    %             params.subA = 2 * pi * params.E0 * ( (1 - cos(alpha2)) - (1 - cos(alpha1)) );
    %
    %             [psProfileVec_Et, psProfileVec_Ef,output] = gpuComputePointSourcesInterpAtDetector(numPoints, D_x, D_z, materialIdx,D_Et, D_Ef, E_t, params, output);
    %
    %             totalEf_profile = totalEf_profile + psProfileVec_Ef;
    %             totalEt_profile = totalEt_profile + psProfileVec_Et;
    %
    %         end
    %         totalEf_profile(2:end) = totalEf_profile(2:end);%.*params.amplitude/numRad;
    %         totalEt_profile(2:end) = totalEt_profile(2:end);%.*params.amplitude/numRad;
    %
    %         D_Ef_t = interp1([0 params.rad], totalEf_profile,D_r);
    %         D_Ef_t(isnan(D_Ef_t))=0; D_Ef_t(isinf(D_Ef_t))=0;
    %         D_Et_t = interp1([0 params.rad], totalEt_profile, D_r);
    %         D_Et_t(isnan(D_Et_t))=0; D_Et_t(isinf(D_Et_t))=0;
    %     end
    
    
    %save results to a structure
    output.A(:,:,materialIdx) = -log10((D_Et_t)./(D_Ef_t));
    output.absSpec(materialIdx) = -log10(sum(D_Et_t(:))/sum(D_Ef_t(:)));
    
    %output.A(:,:,materialIdx) = -log10((D_Et)./(D_Ef));
    %output.absSpec(materialIdx) = -log10(sum(D_Et(:))/sum(D_Ef(:)));
    output.allEtd(:,:,materialIdx) = D_Et_t;
    output.allEfd(:,:,materialIdx) = D_Ef_t;
    
    if materialIdx==1
        toc
    end
    
end