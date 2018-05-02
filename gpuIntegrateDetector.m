function [D_Et, D_Ef] = gpuIntegrateDetector(params, Et_ff, Ef_ff)

%integrate
D_Et = (params.scale)*(abs(Et_ff)).^2;   %this is out_i in cuda bimsim -- the measured intesity
D_Ef = (params.scale)*(abs(Ef_ff)).^2;   %this out_inc in cuda bimsim -- the measured incident field