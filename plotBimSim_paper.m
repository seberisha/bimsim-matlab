%%
%plot absorbance spectra comparison for different number of mc samples
%field order = 100, number of point sources = 1


dir = '/home/sberisha/build/bimsim-bld/bimsimData/pmma/a6.5_fo_100_ps1/mc';

mc = [1 10 100 200 400 1000];
numComp = numel(mc);

h=axes;
lw = 1;
ms = 2;
for ii=1:numComp
    [wav, spec] = loadCudaBimSim(sprintf('%s%i%s',dir,mc(ii),'/absSpec'));
    plot(wav,spec,sprintf('-%c',getMarker(ii)), 'MarkerSize',2)
    legendInfo{ii} = [num2str(mc(ii)) ' mc'];
    hold on
end

set(gca,'FontSize',14)
xlabel('wavenumber')
ylabel('absorbance')
xlim([800 4000])
set(h, 'Xdir','reverse')
legend(legendInfo,'Location','Best')

%%

dir = '/home/sberisha/data/bimsim/matlab/100ps_mc';

mc = [1 10 100];% 200 400 1000];
numComp = numel(mc);

ha=axes;
lw = 1;
ms = 2;
for ii=1:numComp
    matFile = sprintf('%s%i.mat',dir,mc(ii));
    load(matFile)
    plot(wavenumbers,absSpec,sprintf('-%c',getMarker(ii)), 'MarkerSize',2)
    legendInfo{ii} = [num2str(mc(ii)) ' mc'];
    hold on
end

set(ha,'FontSize',14)
xlabel('wavenumber')
ylabel('absorbance')
xlim([800 4000])
set(ha, 'Xdir','reverse')
legend(legendInfo,'Location','Best')
