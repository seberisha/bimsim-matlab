function [D_Et, D_Ef,output] = computePointSources(X,Y,Xi,Yi,interpolate, ...
    startIdx, endIdx, BPF, E_t, i,params,pf_x,pf_z,pf_theta, D_Et,D_Ef,kVec,wavNum, ...
    simRes, origRVecs, Pl_cosalpha1, Pl_cosalpha2,il,r_ps,subA, k_j, normPMinPs, ...
    A,B, psMask,output)


if(interpolate)
    %compute E_f, E_t for one point source
    %pf(2) = ceil(params.a);
    pf(2) = 0;
    pf(1) = pf_x(1);
    pf(3) = pf_z(1);
    
    E_f = gpuComputeEf(params, kVec, wavNum,simRes, pf,origRVecs,Pl_cosalpha1,Pl_cosalpha2,il);
    E_f(isnan(E_f))=0; E_f(isinf(E_f))=0;
    Ef_ps1 = E_f;
    [E_s, E_i] = gpuComputeScatteredFields(params,wavNum, r_ps, params.material(i,:), pf, simRes, subA, k_j, normPMinPs, A, B, psMask);
    E_t(psMask<params.a) = E_i(psMask<params.a);
    E_t(psMask>=params.a) = E_f(psMask>=params.a) + E_s(psMask>=params.a);
    E_t(isnan(E_t))=0; E_t(isinf(E_t))=0;
    Et_ps1 = E_t;
    
    [Et_bpf, Ef_bpf] = applyBandPass(E_t, E_f, BPF);
    
    [Et_ff, Ef_ff] = getFarField(Et_bpf, Ef_bpf,startIdx, endIdx);

    
    %integrate and resample
    D_Et = D_Et + (abs(Et_ff)).^2;   %this is out_i in cuda bimsim -- the measured intesity
    D_Ef = D_Ef + (abs(Ef_ff)).^2;   %this out_inc in cuda bimsim -- the measured incident field
    
    
%     output.Et_ff(:,:,psIdx) = (abs(Et_ff)).^2;
%     output.Ef_ff(:,:,psIdx) = (abs(Ef_ff)).^2;
%     psIdx = psIdx+1;
    
    numPS = numel(pf_x);
    %all point sources
    for p=2:numPS
        pf(1) = pf_x(p);
        pf(3) = pf_z(p);
        
        
        %% use full field simulation
        %         E_f = gpuComputeEf(params, kVec, wavNum,simRes, pf,origRVecs,Pl_cosalpha1,Pl_cosalpha2,il);
        %         E_f(isnan(E_f))=0; E_f(isinf(E_f))=0;
        
        %interpolate E_f
        %E_f = imrotate(single(Ef_ps1), -pf_theta(p)*180/pi,'crop');
        
        % use interp2
        %nodata=find(Xintp<1 | Xintp>maxX | Yintp < 0 | Yintp>maxY);
        ct = cos(pf_theta(p)); st = sin(pf_theta(p));
        Xii = (1+128)/2 + ct*Xi + st*Yi;
        Yii = (1+128)/2 - st*Xi + ct*Yi;
        E_f = interp2(X, Y, (Ef_ps1), Xii, Yii);
        
        E_f(isnan(E_f))=0; E_f(isinf(E_f))=0;
        
        %E_t = imrotate(single(Et_ps1), -pf_theta(p)*180/pi,'crop');
        
        % use interp2
        E_t = interp2(X, Y, (Et_ps1), Xii, Yii);
        
        E_t(isnan(E_t))=0; E_t(isinf(E_t))=0;
        %         [E_s, E_i] = gpuComputeScatteredFields(params,wavNum, r_ps, params.material(i,:), pf, simRes, subA, k_j, normPMinPs, A, B, psMask);
        %         E_t(psMask<params.a) = E_i(psMask<params.a);
        %         E_t(psMask>=params.a) = E_f(psMask>=params.a) + E_s(psMask>=params.a);
        %         E_t(isnan(E_t))=0; E_t(isinf(E_t))=0;
        
        [Et_bpf, Ef_bpf] = applyBandPass(E_t, E_f, BPF);
        
        [Et_ff, Ef_ff] = getFarField(Et_bpf, Ef_bpf,startIdx, endIdx);
        
    
        
        %integrate and resample
        D_Et = D_Et + (abs(Et_ff)).^2;   %this is out_i in cuda bimsim -- the measured intesity
        D_Ef = D_Ef + (abs(Ef_ff)).^2;   %this out_inc in cuda bimsim -- the measured incident field
        
%         output.Et_ff(:,:,psIdx) = (abs(Et_ff)).^2;
%         output.Ef_ff(:,:,psIdx) = (abs(Ef_ff)).^2;
%         psIdx = psIdx+1;
        
        
        
        
    end
else
    numPS = numel(pf_x);
    %all point sources
    for p=1:numPS
        
        pf(1) = pf_x(p);
        pf(2) = 0;
        pf(3) = pf_z(p);
        
        E_f = gpuComputeEf(params, kVec, wavNum,simRes,pf,origRVecs,Pl_cosalpha1,Pl_cosalpha2,il);
        E_f(isnan(E_f))=0; E_f(isinf(E_f))=0;
        
        [E_s, E_i] = gpuComputeScatteredFields(params,wavNum, r_ps, params.material(i,:), pf, simRes, subA, k_j, normPMinPs, A, B, psMask);
        E_t(psMask<params.a) = E_i(psMask<params.a);
        E_t(psMask>=params.a) = E_f(psMask>=params.a) + E_s(psMask>=params.a);
        
        E_t(isnan(E_t))=0; E_t(isinf(E_t))=0;
        
        [Et_bpf, Ef_bpf] = applyBandPass(E_t, E_f, BPF);
        
        [Et_ff, Ef_ff] = getFarField(Et_bpf, Ef_bpf,startIdx, endIdx);
        
        %integrate and resample
        D_Et = D_Et + (abs(Et_ff)).^2;   %this is out_i in cuda bimsim -- the measured intesity
        D_Ef = D_Ef + (abs(Ef_ff)).^2;   %this out_inc in cuda bimsim -- the measured incident field
        
    end
    
end