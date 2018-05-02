function BPF = cpuBPF(gridSize, simRes, objMin, objMax)

df = 1/(gridSize*2);

[iv, iu] = meshgrid(0:simRes-1, 0:simRes-1);
%[iv, iu] = meshgrid(gridPoints, gridPoints);
u=zeros(size(iu)); v = u;
idx = find(iu <= simRes/2);
u(idx) = iu(idx);
idx = find(iu > simRes/2);
u(idx) = (iu(idx) - simRes+1);
u=u.*df;
idx = find(iv <= simRes/2);
v(idx) = iv(idx);
idx = find(iv > simRes/2);
v(idx) = iv(idx) - simRes + 1;
v=v.*df;
fmag = sqrt(u.*u + v.*v);    %compute the magnitude of the frequencies

BPF = (fmag < objMin | fmag > objMax);   %compute the band pass filter

