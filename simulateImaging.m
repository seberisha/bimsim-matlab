%% simulate imaging
%apply the band pass filter
%compute the fft of the field

D = params.gridSize*2;
% Set up range of variables.
u = 0:(params.simRes-1);
v = 0:(params.simRes-1);

% Compute the indices for use in meshgrid
idx = find(u > params.simRes/2);
u(idx) = u(idx) - params.simRes;
idy = find(v > params.simRes/2);
v(idy) = v(idy) - params.simRes;

% Compute the meshgrid arrays
[v, u] = meshgrid(v, u);

df = 1/D;
u=u.*df; v=v.*df;
% this is not correct if the slice goes through the sphere
E_d = fft2(E_t); 
params.objectiveMin = 0.6; params.objectiveMax=1;

fmag = sqrt(u.*u + v.*v);
E_d = fftshift((E_d));

BPF = fftshift((fmag < params.objectiveMin/params.lambda | fmag > params.objectiveMax/params.lambda));
E_d(BPF) = 0;

iftE_d = ifft2(ifftshift(E_d));


%first crop the filtered near-field image of the source and scattered fields
cropSize = padding*params.res;
startIdx = round((params.simRes  - cropSize)/2);
endIdx = startIdx + cropSize-1;
cropEd=iftE_d(startIdx:endIdx,startIdx:endIdx);
cropEf=E_f(startIdx:endIdx, startIdx:endIdx);
%integrate and resample
supersample=1;
I_ed = zeros(params.res/supersample, params.res/supersample);
I_ef = I_ed;

for i=1:supersample
  I_ed = I_ed + (abs(cropEd)).^2;
  I_ef = I_ef + (abs(cropEf)).^2;
end
DE_d = I_ed./(supersample*supersample);
DE_f = I_ef./(supersample*supersample);
figure, imagesc((DE_d)), title('Et_d'),axis image, colormap(brewer), colorbar

%calculate absorbance
A = -log10(DE_d./DE_f);
figure, imagesc((A)), title('A'),axis image, colormap(brewer), colorbar

