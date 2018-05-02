
%% set up and input parameters
clear
gpuDevice(2)
addpath(genpath('~/source/stim-matlab/'));
datDir = '~/data/bimsim/matlab/bulkPMMA_6.5'
addpath(genpath('~/source/stim-matlab/'))

input.resolution = 32;   %resolution at the detector

% standard deviation for a Gaussian extended source --to be estimated
% from FTIR interferogram data - for Agilent FTIR select show raw data in
% lancer control before start of data acquisition (when save dialog box pops up)
input.sigma = 11;

input.interpolate=0;   %option to interpolate for extended sources
input.sphereRadius = 2.25;    %radius of the sphere
input.highMag = 1;    %for high mag mode simulation choose 1, 0 for low mag

%material properties: [lambda(um) real(n) imag(n)]
input.material = [1e4/3502.09657 1.47634 0.000258604];
input.fov = [];
input.spacing=1;   %spacing between point sources in the sampling radius

input.numWav=1;    %number of wavelengths to simulate
input.ps = [0 0 0];    %position of the sphere in x,y,z

input.samples=400; %number of Monte Carlo samples
input.orderEf=50;  %number of orders for the focused/incident field

input.padding = 1; %specify padding
input.E0 = 1;  %amplitude of the field

% condenser NA pair of values specifying an inner obscuration
input.condenserNA_in = .34;
input.condenserNA_out = 0.62;

% objective NA pair of values specifying an inner obscuration
input.objectiveNA_in = .34;
input.objectiveNA_out = 0.62;

%this gives the kVec direction from the y axis
input.theta=1.5708;
input.phi=0;

%focal point
input.pf = [0 0 0];

% create a struct with input parameters for simulation
params = bimSimParamsStruct(input);


%% run forward BimSim
figure
if (input.interpolate)
    
   [D_Et, D_Ef] = gpuInterpolateDetector(params);
      
else
    % run the near-field simulation
    [E_t, E_f] = gpuSimulateScattering(params);
    
    % run the far-field simulation
    [D_Et, D_Ef] = gpuSimulateImaging(E_t, E_f, params);
end

%% save output

output.A = -log10((D_Et)./(D_Ef));
output.absSpec = -log10(sum(D_Et(:))/sum(D_Ef(:)));
output.D_Et = D_Et;
output.D_Ef = D_Ef;
