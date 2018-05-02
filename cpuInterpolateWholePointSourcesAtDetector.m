function [D_Et, D_Ef] = cpuInterpolateWholePointSourcesAtDetector(params)
% This version computes one points source at a certain distance (ring)
% from the location of the sphere. Then it approximates the rest of the
% points sources in the same ring by interpolating the first point source
% at different focal points in the ring.


brewer = brewermap(1000);colormap(brewer)

[X, Y] = meshgrid(1:params.simRes,1:params.simRes);

uc = round(params.simRes/2+1);
vc = round(params.simRes/2+1);

D_Et = 0;
D_Ef = 0;

% X_et = params.x + params.ps(1);
% Y_et = params.z + params.ps(3);

showPointSourceSamples(params.rad, params.numPoints, ones(size(params.rad)));


%compute point source at (0,0,0)

D_Et_ring = 0;
D_Ef_ring = 0;

%focal point for the first point source in the ring
params.pf(1) = params.rad(1);
params.pf(2) = 0;
params.pf(3) = 0;

% run the near-field simulation
[E_t, E_f] = cpuSimulateScattering(params);

%angular spectrum and band pass filtering
[Et_bpf, Ef_bpf] = applyBandPass(E_t, E_f, params.BPF);

[D_Et_ps1, D_Ef_ps1] = cpuIntegrateDetector(params, Et_bpf, Ef_bpf);

D_Et = D_Et + D_Et_ps1;
D_Ef = D_Ef + D_Ef_ps1;
figure

%resample to get far field
[temp_Et, temp_Ef] = getFarField(D_Et, D_Ef,params.startIdx, params.endIdx);
temp_A = -log10((temp_Et)./(temp_Ef));
imagesc(temp_A), axis image,colorbar,colormap(brewer)
drawnow


% subplot(122)
% plot(D_Ef(round(params.simRes/2)+1,round(params.simRes/2)+1:end))
% hold on
% drawnow

for r = 2:params.numRad
    
%     dispvar('ring',r)
    
    theta = 2*pi/params.numPoints(r):2*pi/params.numPoints(r):2*pi;
    
    %focal point for the first point source in the ring
    params.pf(1) = params.rad(r);
    params.pf(2) = 0;
    params.pf(3) = 0;
    
    %     origX_ef = params.x + params.pf(1);
    %     origY_ef = params.z + params.pf(3);
    %
    %     [origTheta, ~, origR] = cart2sph(origX_ef,params.y, origY_ef);
    
    % run the near-field simulation
    [E_t, E_f] = cpuSimulateScattering(params);
    
    %angular spectrum and band pass filtering
    [Et_bpf, Ef_bpf] = applyBandPass(E_t, E_f, params.BPF);
    
    [D_Et_ps1, D_Ef_ps1] = cpuIntegrateDetector(params, Et_bpf, Ef_bpf);
    
    D_Et_ring = D_Et_ring + D_Et_ps1;
    D_Ef_ring = D_Ef_ring + D_Ef_ps1;
    
    for p = 1:params.numPoints(r)        
        
        %use interp2
        newx = uc + (X-uc)*cos(theta(p)) + (Y-vc)*sin(theta(p)) ;
        newy = vc - (X-uc)*sin(theta(p)) + (Y-vc)*cos(theta(p));
        %         newx = params.x + pf_x;
        %         newy = params.z + pf_y;
        
        interpD_Ef = interp2(X, Y, D_Ef_ps1, newx, newy);
        
        interpD_Ef(isnan(interpD_Ef))=0; interpD_Ef(isinf(interpD_Ef))=0;
      
        
        interpD_Et = interp2(X, Y, (D_Et_ps1), newx, newy);
        
        interpD_Et(isnan(interpD_Et))=0; interpD_Et(isinf(interpD_Et))=0;
        
        D_Et_ring = D_Et_ring + params.scale*interpD_Et;   %this is out_i in cuda bimsim -- the measured intesity
        D_Ef_ring = D_Ef_ring + params.scale*interpD_Ef;   %this out_inc in cuda bimsim -- the measured incident field
        
        
    end
    
%     subplot(121)
%     plot(D_Ef_ring(round(params.simRes/2)+1,round(params.simRes/2)+1:end))
%     hold on
%     drawnow
    
    D_Et = D_Et + D_Et_ring;
    %subplot(121), imagesc(D_Et), axis image,colorbar,colormap(brewer)
    %drawnow
    D_Ef = D_Ef + D_Ef_ring;
    
%     subplot(122)
%     plot(D_Ef(round(params.simRes/2)+1,round(params.simRes/2)+1:end))
%     hold on
%     drawnow
        
    %subplot(122), imagesc(D_Ef), axis image,colorbar,colormap(brewer)
    %drawnow
    %resample to get far field
    [temp_Et, temp_Ef] = getFarField(D_Et, D_Ef,params.startIdx, params.endIdx);
    temp_A = -log10((temp_Et)./(temp_Ef));
    imagesc(temp_A), axis image,colorbar,colormap(brewer)
    drawnow
    

end

%resample to get far field
[D_Et, D_Ef] = getFarField(D_Et, D_Ef,params.startIdx, params.endIdx);



