%%
%generate random focal point vectors
function genFocPointVecs(numFP, fov, a)

fpv = -fov + rand([numFP 3])*fov*2;
fpv(:,2) = a/2;
fileName = sprintf('fpv%i_%.1f.mat', numFP, a);
save(fileName, 'fpv')