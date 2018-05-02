function [D_Et, D_Ef,output] = computePointSourcesInterpAtDetector(materialIdx,D_Et, D_Ef, E_t, params, output)

[X, Y] = meshgrid(1:params.res,1:params.res);
uc = round(params.res/2+1);
vc = round(params.res/2+1);

if(params.interpolate)
    
    params.pf(2) = 0;
    params.pf(1) = params.pf_x(1);
    params.pf(3) = params.pf_z(1);
    
    %compute E_f, E_t for one point source
    
    E_f = gpuComputeEf(params, params.kVecs(materialIdx,:));
    E_f(isnan(E_f))=0; E_f(isinf(E_f))=0;
    
    [E_s, E_i] = gpuComputeScatteredFields(params, params.material(materialIdx,:));
    E_t(params.psMask<params.a) = E_i(params.psMask<params.a);
    E_t(params.psMask>=params.a) = E_f(params.psMask>=params.a) + E_s(params.psMask>=params.a);
    E_t(isnan(E_t))=0; E_t(isinf(E_t))=0;
    
    %apply bandpass fparams.iltering
    [Et_params.BPF, Ef_params.BPF] = applyBandPass(E_t, E_f, params.BPF);
    
    %resample to get far field
    [Et_ff, Ef_ff] = getFarField(Et_params.BPF, Ef_params.BPF,params.startIdx, params.endIdx);
    
    %save intensities for the first point source
    D_Et_ps1 = (abs(Et_ff)).^2;
    D_Ef_ps1 = (abs(Ef_ff)).^2;
    
    %integrate to get intensity at the detector
    D_Et_ring = params.scale*D_Et_ps1;
    %D_Et = D_Et + params.scale*D_Et_ps1;   %this is out_i in cuda bimsim -- the measured intesity
    %D_Ef = D_Ef + params.scale*D_Ef_ps1;   %this out_inc in cuda bimsim -- the measured incident field
    D_Ef_ring = params.scale*D_Ef_ps1;
    
    numPS = numel(params.pf_x);
    
    % interpolate the intensities at the detector for the rest of the point
    % sources in the same ring
    
    for p=2:numPS
        
        params.pf(1) = params.pf_x(p);
        params.pf(3) = params.pf_z(p);
      
        %use interp2
        newx = uc + (X-uc)*cos(params.pf_theta(p)) + (Y-vc)*sin(params.pf_theta(p)) ;
        newy = vc - (X-uc)*sin(params.pf_theta(p)) + (Y-vc)*cos(params.pf_theta(p));
        
        interpD_Ef = interp2(X, Y, (D_Ef_ps1), newx, newy);
       
        interpD_Ef(isnan(interpD_Ef))=0; interpD_Ef(isinf(interpD_Ef))=0;
        
        %use interp2
        interpD_Et = interp2(X, Y, (D_Et_ps1), newx, newy);
        
        interpD_Et(isnan(interpD_Et))=0; interpD_Et(isinf(interpD_Et))=0;
        
        %integrate and resample
        %D_Et = D_Et + params.scale*interpD_Et;   %this is out_i in cuda bimsim -- the measured intesity
        %D_Ef = D_Ef + params.scale*interpD_Ef;   %this out_inc in cuda bimsim -- the measured incident field
        
        D_Et_ring = D_Et_ring + params.scale*interpD_Et;   %this is out_i in cuda bimsim -- the measured intesity
        D_Ef_ring = D_Ef_ring + params.scale*interpD_Ef;   %this out_inc in cuda bimsim -- the measured incident field
        
    end
    
    D_Et = D_Et + D_Et_ring;
    D_Ef = D_Ef + D_Ef_ring;
    
else
    numPS = numel(params.pf_x);
    %all point sources
    for p=1:numPS
        
        pf(1) = params.pf_x(p);
        pf(2) = 0;
        pf(3) = params.pf_z(p);
        
        E_f = gpuComputeEf(params);
        E_f(isnan(E_f))=0; E_f(isinf(E_f))=0;
        
        [E_s, E_i] = gpuComputeScatteredFields(params, params.material(materialIdx,:));
        E_t(params.psMask<params.a) = E_i(params.psMask<params.a);
        E_t(params.psMask>=params.a) = E_f(params.psMask>=params.a) + E_s(params.psMask>=params.a);
        
        E_t(isnan(E_t))=0; E_t(isinf(E_t))=0;
        
        [Et_params.BPF, Ef_params.BPF] = applyBandPass(E_t, E_f, params.BPF);
        
        [Et_ff, Ef_ff] = getFarField(Et_params.BPF, Ef_params.BPF,params.startIdx, params.endIdx);
        
        %integrate and resample
        D_Et = D_Et + (abs(Et_ff)).^2;   %this is out_i in cuda bimsim -- the measured intesity
        D_Ef = D_Ef + (abs(Ef_ff)).^2;   %this out_inc in cuda bimsim -- the measured incident field
        
    end
    
end