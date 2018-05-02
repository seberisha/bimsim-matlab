function params = bimSimParamsStruct(input)

params = input;

if input.highMag
    pixRes = 1.1;   %pixel resolution in um for high mag
else
    pixRes = 5.5;  %pixel resolution in um for low mag
end

if isempty(input.fov)
    params.fov = round(input.resolution*pixRes);  %FOV in um
end

% this is ok for low mage mode, e.g. for resolution of 128 (at the detector) the field of
% view will be 141 um => sampling radius for interpolation will go from 0
% to 63; actually to cover all of FOV it should go up to 69
% However, this gets worse in high mag mode. In this case the FOV would be
% 704 um while sampling radius range is [0 63].
% Possible solution: increase the size of the detector image at far field,
% or don't resample at all at far field. Do the resampling at the end after
% interpolation is completely done.
params.rad = 0:round(params.fov/2-1);

params.s=rng;   %control random number generation

% point sources
params.numRad = numel(params.rad);

params.gaussianAmplitude = exp(-(([input.spacing:input.spacing:params.numRad].^2))/(2*(input.sigma^2)))/(2*pi*(input.sigma^2));
params.gaussianAmplitude = params.gaussianAmplitude./max(params.gaussianAmplitude(:)); %normalize the amplitudes

params.gaussianAmplitude = params.gaussianAmplitude(1:input.spacing:end);

%spatial resolution of the field plane
params.simRes = input.resolution*(2*input.padding + 1);

% params.alpha1 and params.alpha2 from condenserNA_in and condenserNA_out, respectively
params.alpha1 = asin(input.condenserNA_in); 
params.alpha2 = asin(input.condenserNA_out);
ordVecEf=gpuArray((0:input.orderEf)'); %the prefix term (2l + 1)*i^l
params.il = 1i.^ordVecEf;
params.il = reshape(params.il,[1 1 input.orderEf+1]);

%compute the crop size needed for resampling at the detector
cropSize = input.padding*input.resolution;

%indices for resampling at far field
if input.padding==0
    params.startIdx=1;
    params.endIdx=params.simRes;
else
    params.startIdx = fix(params.simRes/2 + 1) - floor(cropSize/2);
    params.endIdx = params.startIdx + cropSize-1;
end

%direction of the incident light
[x,y,z] = sph2cart(input.theta,input.phi,1);
lightDirection  = [x y z];

%% matlab coordinates
% get r and rVecs
%create a grid of points representing pixel positions in the field plane

halfGridSize = round(params.fov/2)*(2*input.padding + 1);
gx = linspace(-fix(halfGridSize),ceil(halfGridSize)-1,params.simRes);
gy = gx;

[x,z] = meshgrid(gx, gy); % field slice in the x z plane

%field plane y at the radius of the sphere
y = ones(params.simRes,params.simRes)*input.sphereRadius;
params.x = x;
params.y = y;
params.z = z;

%% vectors corresponding to each point in the plane
rVecs = zeros(params.simRes*params.simRes, 3);
%r value at each pixel position
rVecs(:,1) = x(:); rVecs(:,2) = y(:); rVecs(:,3) = z(:);

%vectors of each point with respect to the distance from center of the
%sphere
rVecs_ps = bsxfun(@minus, rVecs, input.ps);

%mask used later for separating the internal and external fields in total
%field calculation
params.psMask=reshape(sqrt(sum(rVecs_ps.^2,2)),params.simRes, params.simRes);

%legendre polynomials
params.Pl_cosalpha1 = gpuLegendre(input.orderEf+1,cos(params.alpha1));
params.Pl_cosalpha2 = gpuLegendre(input.orderEf+1,cos(params.alpha2));

%remember the points (in vector and distance representation) with respect
%to the center of the sphere
params.origRVecs = rVecs;
%normalized vectors with respect to distance from sphere center
params.normPMinPs = bsxfun(@rdivide, rVecs_ps,  sqrt(sum(rVecs_ps.^2,2)));
%r value at each pixel position with respect to the center of the sphere
params.r_ps=reshape(sqrt(sum(rVecs_ps.^2,2)),params.simRes, params.simRes);

params.lambda = input.material(1);  %wavelengths in micrometers
%magnitude of k vectors for each wavelength
params.magKVector = 2*pi./params.lambda;
%k vectors for each wavelength
params.kVec = lightDirection.*params.magKVector;

params.subA = 2 * pi * input.E0 * ( (1 - cos(params.alpha2)) - (1 - cos(params.alpha1)) );

%band pass filter
params.BPF = gpuBPF(halfGridSize, params.simRes, input.objectiveNA_in/params.lambda,  input.objectiveNA_out/params.lambda);
params.scale = 1;


