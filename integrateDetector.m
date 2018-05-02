function integrateDetector()

%integrate and resample
D_Et = D_Et + (abs(cropEt_d)).^2;   %this is out_i in cuda bimsim -- the measured intesity
D_Ef = D_Ef + (abs(cropEf_d)).^2;   %this out_inc in cuda bimsim -- the measured incident field
