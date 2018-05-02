%% paths
addpath(genpath('~/source/stim-matlab/'))
addpath(genpath('~/source/bimsim-matlab/'))
addpath(genpath('~/source/stimlib/stim/matlab/'))

%%
brewer = brewermap(1000);
%dataPath = '/Users/sbstn/data/ftir/all_ftir/ftir/pmma/5-20um/sphere5_20_middle'
%dataPath = '/Users/sbstn/data/ftir/all_ftir/ftir/pmma/5-27um/sphere1_middle'
%dataPath = '/Users/sbstn/data/ftir/all_ftir/ftir/pmma/6_5/pmma6_5_middle'
%dataPath = '/Users/sbstn/data/ftir/all_ftir/ftir/pmma/d90_sd2_15/middle_8r/spec'
%dataPath = '/Users/sbstn/data/ftir/all_ftir/ftir/oneSphere/spec_128ca_8r_sphere2'
dataPath = '/Users/sbstn/data/ftir/all_ftir/ftir/oneSphere/spec_128ca_8r_sphere1'
%dataPath = '/Users/sbstn/data/ftir/all_ftir/ftir/pmma_20-27/high_mag/sphere_1'

[wav, absImages, spec] = loadFTIR([dataPath '.hdr'], [dataPath '.dat'], 0, [], []);
figure, imagesc(absImages(:,:,end)), axis image, colormap(brewer)

%% get spectrum for a pixel
%r_idx = 59; c_idx = 64;
%r_idx = 61; c_idx = 63;
%r_idx = 60; c_idx = 65;
%r_idx = 62; c_idx = 64;
%r_idx = 66; c_idx = 64;
r_idx = 66; c_idx = 60;
%r_idx = 65; c_idx = 63;

spec = squeeze(absImages(r_idx, c_idx,:));

%% plot
figure
ax=axes('FontSize',14);
plot(wav,spec,'LineWidth',2)
set(ax, 'Xdir','reverse')
xlabel('Wavenumber (cm^{-1})')
ylabel('Absorbance')
set(gca,'FontSize',14)
xlim([1000 3500])

%% plot the subplot of light scattering example
load ~/data/refIdx/etaPmma.mat
figure
subplot(3,1,1)
plot(etaPmma(:,1),etaPmma(:,3),'LineWidth',2)
set(gca, 'Xdir','reverse')
xlabel('Wavenumber (cm^{-1})')
ylabel('Absorbance')
set(gca,'FontSize',14)
xlim([1000 3500])

dataPath = '/Users/sbstn/data/ftir/all_ftir/ftir/oneSphere/spec_128ca_8r_sphere1'
[wav, absImages, spec] = loadFTIR([dataPath '.hdr'], [dataPath '.dat'], 0, [], []);
r_idx = 66; c_idx = 60;
spec = squeeze(absImages(r_idx, c_idx,:));

subplot(3,1,2)
plot(wav,spec,'LineWidth',2)
set(gca, 'Xdir','reverse')
set(gca,'xtick',[])
xlabel('Wavenumber (cm^{-1})')
ylabel('Absorbance')
set(gca,'FontSize',14)
xlim([1000 3500])

dataPath = '/Users/sbstn/data/ftir/all_ftir/ftir/pmma/5-20um/sphere5_20_middle'
[wav, absImages, spec] = loadFTIR([dataPath '.hdr'], [dataPath '.dat'], 0, [], []);
r_idx = 59; c_idx = 64;
spec = squeeze(absImages(r_idx, c_idx,:));

subplot(3,1,3)
plot(wav,spec,'LineWidth',2)
set(gca, 'Xdir','reverse')

xlabel('Wavenumber (cm^{-1})')
ylabel('Absorbance')
set(gca,'FontSize',14)
xlim([1000 3500])

saveFigMF(gcf, '~/Desktop/conferences/2016/MM/Poster/tes/Figures/scatteringSpectrumExample', 'm', 'eps', 'fig')


%% plot complex refrective index for PMMA
fontSize = 24;
figure
pu = get(gcf,'PaperUnits');
pp = get(gcf,'PaperPosition');
set(gcf,'Units',pu,'Position',pp)

ppi = get(0,'ScreenPixelsPerInch');
set(0,'ScreenPixelsPerInch',ppi*1.7)
etaPMMA = csvread('~/data/refIdx/etaPMMA_forMatlab.txt');
[hAx,hLine1,hLine2] = plotyy(etaPMMA(:,1),etaPMMA(:,2),etaPMMA(:,1),etaPMMA(:,3));
set(hAx(1), 'Xdir','reverse')
set(hAx(2), 'Xdir','reverse')
set(hLine1,'linewidth',2)
set(hLine2,'linewidth',2)

% plot(etaPMMA(:,1), etaPMMA(:,2))
% hold on
% plot(etaPMMA(:,1), etaPMMA(:,3))
%yyaxis left

xlabel('Wavenumber (cm^{-1})','interpreter','tex')
ylabel(hAx(1),'Effective Refractive Index') % left y-axis

ylabel(hAx(2),'Effective Refractive Index') % right y-axis
set(hAx(1),'FontSize',fontSize)
set(hAx(2),'FontSize',fontSize)

xlim(hAx(1), [1000 3500])
xlim(hAx(2), [1000 3500])

ylim(hAx(1),[1.3 1.55])
ylim(hAx(2),[0 0.2])
h_legend = legend('\eta, optical path length', '\kappa, absorption coefficient', 'Location', 'Best')
set(h_legend,'FontSize',fontSize)
grid on
%set(gca,'LooseInset',get(gca,'TightInset'));
%set(gcf, 'Renderer', 'opengl');
%set(gcf,'PaperPositionMode','auto')


%% plot complex refrective index for polystyrene
figure
etaPS = csvread('~/data/refIdx/etaPolystyrene_forMatlab.txt');
[hAx,hLine1,hLine2] = plotyy(etaPS(:,1),etaPS(:,2),etaPS(:,1),etaPS(:,3));
set(hAx(1), 'Xdir','reverse')
set(hAx(2), 'Xdir','reverse')
set(hLine1,'linewidth',4)
set(hLine2,'linewidth',4)

% plot(etaPMMA(:,1), etaPMMA(:,2))
% hold on
% plot(etaPMMA(:,1), etaPMMA(:,3))
%yyaxis left

xlabel('Wavenumber (cm^{-1})','fontweight','bold') 
ylabel(hAx(1),'Effective Refractive Index','fontweight','bold') % left y-axis

ylabel(hAx(2),'Effective Refractive Index','fontweight','bold') % right y-axis
set(hAx(1),'FontSize',fontSize,'fontweight','bold')
set(hAx(2),'FontSize',fontSize,'fontweight','bold')

xlim(hAx(1), [1000 3500])
xlim(hAx(2), [1000 3500])

ylim(hAx(1),[0.9 1.8])
ylim(hAx(2),[-0.02 1])
h_legend = legend('\eta, optical path length', '\kappa, absorption coefficient', 'Location', 'Best')
set(h_legend,'FontSize',fontSize, 'fontweight','bold')
grid on

%%
fontsize = 14;
set(0,'defaultaxesfontsize',fontsize);
set(0,'defaulttextfontsize',fontsize);
set(0,'defaultaxesfontweight','bold');
set(0,'defaulttextfontweight','bold');

%% plot results of pmma d = 3_36, simulation a = 1.5
dataPath = '~/data/ftir/pmma/3_36/pmma3_36_middle_highmag';
[wav, absImages, spec] = loadFTIR([dataPath '.hdr'], [dataPath '.dat'], 2, [], []);
spec = squeeze(absImages(60,63,:));
[wav_b, spec_b, A_b, int_b, inc_b] = loadCudaBimSim('~/data/bimsim/cuda/a1_5_fo_32_mc1024_highMag_fpAta_ps48_141_nBP_allWav_s0_2_n1_4_k_pmma/','absSpec', 'out_a','out_i', 'out_inc', 128,128, 2);
absImages = absImages(:,:,9:8:end);
wav_b = wav_b(2:end); wav = wav(9:8:end);
spec_b = spec_b(2:end); spec = spec(9:8:end);
figure
plot(wav,spec, 'LineWidth',2); hold on; plot(wav_b,spec_b, 'LineWidth',2)
xlim(gca, [1000 3500])
ylim(gca,[0 1])
xlabel('\textbf{Wavenumber ($\mathbf{cm^{-1}}$)}', 'Interpreter', 'latex') 
ylabel(gca,'Absorbance','fontweight','bold') 
set(gca,'FontSize',14,'fontweight','bold')
h_legend = legend('measured', 'predicted', 'Location', 'Best');
set(h_legend,'FontSize',14, 'fontweight','bold')
grid on
%% save spectra results
saveas(gcf,'~/source/latex/bimsim/figures/specComp_pmma_1_5','epsc')
saveFigMF(gcf, '~/source/latex/bimsim/figures/source/specComp_pmma_1_5','m','fig')

%% crop and save absorbance images
cr = 64; sz = 128;
idx = 30;
t_m = absImages(:,:,idx);
%t_m(t_m<0) = 0;
t_m = t_m(round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1,round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1);
t_p = A_b(:,:,30);
%t_p(t_p<0) = 0;
t_p = t_p(round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1,round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1);
figure; imagesc(t_m),axis image, colormap(brewer),colorbar, axis off
%xlabel(['{\lambda} = ', num2str(1e4/wav(idx))],'Interpreter','latex')
set(gca,'FontSize',14,'fontweight','bold')
hold on
plot([60; 50], [60; 60], '-k', 'LineWidth', 5)
text(55,57, '11$\mu m$', 'HorizontalAlignment','center','Interpreter','latex','FontWeight','bold','FontSize',25,'Color','k')
saveFigMF(gcf, '~/source/latex/bimsim/figures/source/m_pmma_a1_5_wav_6_8','m','fig')

figure; imagesc(t_p),axis image, colormap(brewer),colorbar, axis off
%xlabel(['\lambda = ', num2str(1e4/wav(idx))],'Interpreter','latex')
set(gca,'FontSize',14,'fontweight','bold')
hold on
plot([60; 50], [60; 60], '-w', 'LineWidth', 5)
text(55,57, '11$\mu m$', 'HorizontalAlignment','center','Interpreter','latex','FontWeight','bold','FontSize',25,'Color','w')
saveFigMF(gcf, '~/source/latex/bimsim/figures/source/p_pmma_a1_5_wav_6_8','m','fig')

%% crop and save absorbance images
cr = 64; sz = 128;
idx = 60;
t_m = absImages(:,:,idx);
%t_m(t_m<0) = 0;
t_m = t_m(round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1,round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1);
t_p = A_b(:,:,idx);
%t_p(t_p<0) = 0;
t_p = t_p(round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1,round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1);
figure; imagesc(t_m),axis image, colormap(brewer),colorbar, axis off
%xlabel(['{\lambda} = ', num2str(1e4/wav(idx))],'Interpreter','latex')
set(gca,'FontSize',14,'fontweight','bold')
hold on
plot([60; 50], [60; 60], '-k', 'LineWidth', 5)
text(55,57, '11$\mu m$', 'HorizontalAlignment','center','Interpreter','latex','FontWeight','bold','FontSize',25,'Color','k')
saveFigMF(gcf, '~/source/latex/bimsim/figures/source/m_pmma_a1_5_wav_5_2','m','fig')

figure; imagesc(t_p),axis image, colormap(brewer),colorbar, axis off
%xlabel(['\lambda = ', num2str(1e4/wav(idx))],'Interpreter','latex')
set(gca,'FontSize',14,'fontweight','bold')
hold on
plot([60; 50], [60; 60], '-w', 'LineWidth', 5)
text(55,57, '11$\mu m$', 'HorizontalAlignment','center','Interpreter','latex','FontWeight','bold','FontSize',25,'Color','w')
saveFigMF(gcf, '~/source/latex/bimsim/figures/source/p_pmma_a1_5_wav_5_2','m','fig')

%% crop and save absorbance images
cr = 64; sz = 128;
idx = 90;
t_m = absImages(:,:,idx);
%t_m(t_m<0) = 0;
t_m = t_m(round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1,round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1);
t_p = A_b(:,:,idx);
%t_p(t_p<0) = 0;
t_p = t_p(round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1,round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1);
figure; imagesc(t_m),axis image, colormap(brewer),colorbar, axis off
%xlabel(['{\lambda} = ', num2str(1e4/wav(idx))],'Interpreter','latex')
set(gca,'FontSize',14,'fontweight','bold')
hold on
plot([60; 50], [60; 60], '-k', 'LineWidth', 5)
text(55,57, '11$\mu m$', 'HorizontalAlignment','center','Interpreter','latex','FontWeight','bold','FontSize',25,'Color','k')
saveFigMF(gcf, '~/source/latex/bimsim/figures/source/m_pmma_a1_5_wav_4_2','m','fig')

figure; imagesc(t_p),axis image, colormap(brewer),colorbar, axis off
%xlabel(['\lambda = ', num2str(1e4/wav(idx))],'Interpreter','latex')
set(gca,'FontSize',14,'fontweight','bold')
hold on
plot([60; 50], [60; 60], '-w', 'LineWidth', 5)
text(55,57, '11$\mu m$', 'HorizontalAlignment','center','Interpreter','latex','FontWeight','bold','FontSize',25,'Color','w')
saveFigMF(gcf, '~/source/latex/bimsim/figures/source/p_pmma_a1_5_wav_4_2','m','fig')

%% crop and save absorbance images
cr = 64; sz = 128;
idx = 120;
t_m = absImages(:,:,idx);
%t_m(t_m<0) = 0;
t_m = t_m(round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1,round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1);
t_p = A_b(:,:,idx);
%t_p(t_p<0) = 0;
t_p = t_p(round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1,round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1);
figure; imagesc(t_m),axis image, colormap(brewer),colorbar, axis off
%xlabel(['{\lambda} = ', num2str(1e4/wav(idx))],'Interpreter','latex')
set(gca,'FontSize',14,'fontweight','bold')
hold on
plot([60; 50], [60; 60], '-k', 'LineWidth', 5)
text(55,57, '11$\mu m$', 'HorizontalAlignment','center','Interpreter','latex','FontWeight','bold','FontSize',25,'Color','k')
saveFigMF(gcf, '~/source/latex/bimsim/figures/source/m_pmma_a1_5_wav_3_5','m','fig')

figure; imagesc(t_p),axis image, colormap(brewer),colorbar, axis off
%xlabel(['\lambda = ', num2str(1e4/wav(idx))],'Interpreter','latex')
set(gca,'FontSize',14,'fontweight','bold')
hold on
plot([60; 50], [60; 60], '-w', 'LineWidth', 5)
text(55,57, '11$\mu m$', 'HorizontalAlignment','center','Interpreter','latex','FontWeight','bold','FontSize',25,'Color','w')
saveFigMF(gcf, '~/source/latex/bimsim/figures/source/p_pmma_a1_5_wav_3_5','m','fig')



%% plot results of pmma d = spec_128ca_8r_sphere1, simulation a = 2.2 -- pmma with d = 6.5
dataPath = '~/data/ftir/oneSphere/spec_128ca_8r_sphere1';
dataGenCudaPath = '~/data/bimsim/cuda/';
[wav, absImages, spec] = loadFTIR([dataPath '.hdr'], [dataPath '.dat'], 2, [], []);
spec = squeeze(absImages(65,60,:));
[wav_b, spec_b, A_b, int_b, inc_b] = loadCudaBimSim([dataGenCudaPath 'pmma_a2_2_fo_32_mc1024_highMag_fpAta_ps64_141_nBP_allWav_s0_15_n1_48_k/'],'absSpec', 'out_a','out_i', 'out_inc', 128,128, 2);
absImages = absImages(:,:,15:4:end);
wav = wav(15:4:end);
spec = spec(15:4:end);
figure
plot(wav,spec, 'LineWidth',2); hold on; plot(wav_b,spec_b, 'LineWidth',2)
xlim(gca, [1000 3500])
ylim(gca,[0 1])
xlabel('\textbf{Wavenumber ($\mathbf{cm^{-1}}$)}', 'Interpreter', 'latex') 
ylabel(gca,'Absorbance','fontweight','bold') 
set(gca,'FontSize',14,'fontweight','bold')
h_legend = legend('measured', 'predicted', 'Location', 'Best');
set(h_legend,'FontSize',14, 'fontweight','bold')
grid on
%% save spectra results
saveas(gcf,'~/source/latex/bimsim/figures/specComp_pmma_2_2','epsc')
saveFigMF(gcf, '~/source/latex/bimsim/figures/source/specComp_pmma_2_2','m','fig')

%% crop and save absorbance images
cr = 64; sz = 128;
idx = 5;
t_m = absImages(:,:,idx);
%t_m(t_m<0) = 0;
t_m = t_m(round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1,round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1);
t_p = A_b(:,:,idx);
%t_p(t_p<0) = 0;
t_p = t_p(round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1,round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1);
figure; imagesc(t_m),axis image, colormap(brewer),colorbar, axis off
%xlabel(['{\lambda} = ', num2str(1e4/wav(idx))],'Interpreter','latex')
set(gca,'FontSize',14,'fontweight','bold')
hold on
plot([60; 50], [60; 60], '-w', 'LineWidth', 5)
text(55,57, '11$\mu m$', 'HorizontalAlignment','center','Interpreter','latex','FontWeight','bold','FontSize',25,'Color','w')
%saveFigMF(gcf, '~/source/latex/bimsim/figures/source/m_pmma_a2_2_wav_6_9','m','fig')

figure; imagesc(t_p),axis image, colormap(brewer),colorbar, axis off
%xlabel(['\lambda = ', num2str(1e4/wav(idx))],'Interpreter','latex')
set(gca,'FontSize',14,'fontweight','bold')
hold on
plot([60; 50], [60; 60], '-w', 'LineWidth', 5)
text(55,57, '11$\mu m$', 'HorizontalAlignment','center','Interpreter','latex','FontWeight','bold','FontSize',25,'Color','w')
%saveFigMF(gcf, '~/source/latex/bimsim/figures/source/p_pmma_a1_5_wav_6_8','m','fig')


%% crop and save absorbance images
cr = 64; sz = 128;
idx = 15;
t_m = absImages(:,:,idx);
%t_m(t_m<0) = 0;
t_m = t_m(round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1,round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1);
t_p = A_b(:,:,idx);
%t_p(t_p<0) = 0;
t_p = t_p(round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1,round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1);
figure; imagesc(t_m),axis image, colormap(brewer),colorbar, axis off
%xlabel(['{\lambda} = ', num2str(1e4/wav(idx))],'Interpreter','latex')
set(gca,'FontSize',14,'fontweight','bold')
hold on
plot([60; 50], [60; 60], '-w', 'LineWidth', 5)
text(55,57, '11$\mu m$', 'HorizontalAlignment','center','Interpreter','latex','FontWeight','bold','FontSize',25,'Color','w')
%saveFigMF(gcf, '~/source/latex/bimsim/figures/source/m_pmma_a2_2_wav_6_9','m','fig')

figure; imagesc(t_p),axis image, colormap(brewer),colorbar, axis off
%xlabel(['\lambda = ', num2str(1e4/wav(idx))],'Interpreter','latex')
set(gca,'FontSize',14,'fontweight','bold')
hold on
plot([60; 50], [60; 60], '-w', 'LineWidth', 5)
text(55,57, '11$\mu m$', 'HorizontalAlignment','center','Interpreter','latex','FontWeight','bold','FontSize',25,'Color','w')
%saveFigMF(gcf, '~/source/latex/bimsim/figures/source/p_pmma_a1_5_wav_6_8','m','fig')

%% crop and save absorbance images
cr = 64; sz = 128;
idx = 25;
t_m = absImages(:,:,idx);
%t_m(t_m<0) = 0;
t_m = t_m(round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1,round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1);
t_p = A_b(:,:,idx);
%t_p(t_p<0) = 0;
t_p = t_p(round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1,round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1);
figure; imagesc(t_m),axis image, colormap(brewer),colorbar, axis off
%xlabel(['{\lambda} = ', num2str(1e4/wav(idx))],'Interpreter','latex')
set(gca,'FontSize',14,'fontweight','bold')
hold on
plot([60; 50], [60; 60], '-w', 'LineWidth', 5)
text(55,57, '11$\mu m$', 'HorizontalAlignment','center','Interpreter','latex','FontWeight','bold','FontSize',25,'Color','w')
%saveFigMF(gcf, '~/source/latex/bimsim/figures/source/m_pmma_a2_2_wav_6_9','m','fig')

figure; imagesc(t_p),axis image, colormap(brewer),colorbar, axis off
%xlabel(['\lambda = ', num2str(1e4/wav(idx))],'Interpreter','latex')
set(gca,'FontSize',14,'fontweight','bold')
hold on
plot([60; 50], [60; 60], '-w', 'LineWidth', 5)
text(55,57, '11$\mu m$', 'HorizontalAlignment','center','Interpreter','latex','FontWeight','bold','FontSize',25,'Color','w')
%saveFigMF(gcf, '~/source/latex/bimsim/figures/source/p_pmma_a1_5_wav_6_8','m','fig')

%% plot cluster of spectra around center pmma 6.5dataPath = '~/data/ftir/oneSphere/spec_128ca_8r_sphere1';
dataGenCudaPath = '~/data/bimsim/cuda/';
[wav, absImages, spec] = loadFTIR([dataPath '.hdr'], [dataPath '.dat'], 2, [], []);
spec = squeeze(absImages(65,60,:));
[wav_b, spec_b, A_b, int_b, inc_b] = loadCudaBimSim([dataGenCudaPath 'pmma_a2_2_fo_32_mc1024_highMag_fpAta_ps64_141_nBP_allWav_s0_15_n1_48_k/'],'absSpec', 'out_a','out_i', 'out_inc', 128,128, 2);
absImages = absImages(:,:,15:4:end);
wav = wav(15:4:end);
spec = squeeze(absImages(65,59,:));
figure
p1 = plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(65,60,:));
hold on
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(65,61,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(66,59,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(66,60,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(66,61,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(67,59,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(67,60,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(67,61,:));
plot(wav,spec, 'b', 'LineWidth', 2)

spec_b = squeeze(A_b(64,64,:));
p2 = plot(wav_b,spec_b, 'r--','LineWidth',2)
hold on
spec_b = squeeze(A_b(64,65,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(64,66,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(65,64,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(65,65,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(65,66,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(66,64,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(66,65,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(66,66,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)

xlim(gca, [1000 3500])
ylim(gca,[0 1])
xlabel('\textbf{Wavenumber ($\mathbf{cm^{-1}}$)}', 'Interpreter', 'latex') 
ylabel(gca,'Absorbance','fontweight','bold') 
set(gca,'FontSize',14,'fontweight','bold')
h_legend = legend([p1 p2],'measured', 'predicted', 'Location', 'Best');
set(h_legend,'FontSize',14, 'fontweight','bold')
grid on


%% plot cluster of spectra for pmma of d = 3.36
dataPath = '~/data/ftir/pmma/3_36/pmma3_36_middle_highmag';
[wav, absImages, spec] = loadFTIR([dataPath '.hdr'], [dataPath '.dat'], 2, [], []);
[wav_b, spec_b, A_b, int_b, inc_b] = loadCudaBimSim('~/data/bimsim/cuda/a1_5_fo_32_mc1024_highMag_fpAta_ps48_141_nBP_allWav_s0_2_n1_4_k_pmma/','absSpec', 'out_a','out_i', 'out_inc', 128,128, 2);
absImages = absImages(:,:,9:8:end);
wav_b = wav_b(2:end); wav = wav(9:8:end);
A_b = A_b(:,:,2:end);

spec = squeeze(absImages(59,62,:));
figure
p1 = plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(59,63,:));
hold on
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(59,64,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(60,62,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(60,63,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(60,64,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(61,62,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(61,63,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(61,64,:));
plot(wav,spec, 'b', 'LineWidth', 2)

spec_b = squeeze(A_b(64,64,:));
p2 = plot(wav_b,spec_b, 'r--','LineWidth',2)
hold on
spec_b = squeeze(A_b(64,65,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(64,66,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(65,64,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(65,65,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(65,66,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(66,64,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(66,65,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(66,66,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)

xlim(gca, [1000 3500])
ylim(gca,[0 1])
xlabel('\textbf{Wavenumber ($\mathbf{cm^{-1}}$)}', 'Interpreter', 'latex') 
ylabel(gca,'Absorbance','fontweight','bold') 
set(gca,'FontSize',14,'fontweight','bold')
h_legend = legend([p1 p2],'measured', 'predicted', 'Location', 'Best');
set(h_legend,'FontSize',14, 'fontweight','bold')
grid on


%% plot results of pmma d = spec_128ca_8r_sphere1, simulation a = 2.2 -- ps with d = 4.52+-0.15
[wav, absImages, spec] = loadFTIR('/home/sberisha/data/ftir/polystyrene/polystyrene_4_52/polystyrene_4_52.hdr', '/home/sberisha/data/ftir/polystyrene/polystyrene_4_52/polystyrene_4_52',3,[],[]);

spec = squeeze(absImages(66,64,:));
[wav_b, spec_b, A_b, int_b, inc_b] = loadCudaBimSim([dataGenCudaPath 'ps_a2_2_fo_32_mc1024_highMag_fpAta_ps76_141_nBP_allWav_s5_n1_4/'],'absSpec', 'out_a','out_i', 'out_inc', 128,128, 2);

absImages = absImages(:,:,2:2:end);
wav = wav(2:2:end);
spec = spec(2:2:end);
figure
plot(wav,spec, 'LineWidth',2); hold on; plot(wav_b,spec_b, 'LineWidth',2)
xlim(gca, [1000 3500])
ylim(gca,[min(spec) 1])
xlabel('\textbf{Wavenumber ($\mathbf{cm^{-1}}$)}', 'Interpreter', 'latex') 
ylabel(gca,'Absorbance','fontweight','bold') 
set(gca,'FontSize',14,'fontweight','bold')
h_legend = legend('measured', 'predicted', 'Location', 'Best');
set(h_legend,'FontSize',14, 'fontweight','bold')
grid on

%% crop and save absorbance images
cr = 64; sz = 128;
idx = 5;
t_m = absImages(:,:,idx);
%t_m(t_m<0) = 0;
t_m = t_m(round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1,round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1);
t_p = A_b(:,:,idx);
%t_p(t_p<0) = 0;
t_p = t_p(round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1,round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1);
figure; imagesc(t_m),axis image, colormap(brewer),colorbar, axis off
%xlabel(['{\lambda} = ', num2str(1e4/wav(idx))],'Interpreter','latex')
set(gca,'FontSize',14,'fontweight','bold')
hold on
plot([60; 50], [60; 60], '-w', 'LineWidth', 5)
text(55,57, '11$\mu m$', 'HorizontalAlignment','center','Interpreter','latex','FontWeight','bold','FontSize',25,'Color','w')
%saveFigMF(gcf, '~/source/latex/bimsim/figures/source/m_pmma_a2_2_wav_6_9','m','fig')

figure; imagesc(t_p),axis image, colormap(brewer),colorbar, axis off
%xlabel(['\lambda = ', num2str(1e4/wav(idx))],'Interpreter','latex')
set(gca,'FontSize',14,'fontweight','bold')
hold on
plot([60; 50], [60; 60], '-w', 'LineWidth', 5)
text(55,57, '11$\mu m$', 'HorizontalAlignment','center','Interpreter','latex','FontWeight','bold','FontSize',25,'Color','w')
%saveFigMF(gcf, '~/source/latex/bimsim/figures/source/p_pmma_a1_5_wav_6_8','m','fig')


%% crop and save absorbance images
cr = 64; sz = 128;
idx = 15;
t_m = absImages(:,:,idx);
%t_m(t_m<0) = 0;
t_m = t_m(round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1,round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1);
t_p = A_b(:,:,idx);
%t_p(t_p<0) = 0;
t_p = t_p(round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1,round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1);
figure; imagesc(t_m),axis image, colormap(brewer),colorbar, axis off
%xlabel(['{\lambda} = ', num2str(1e4/wav(idx))],'Interpreter','latex')
set(gca,'FontSize',14,'fontweight','bold')
hold on
plot([60; 50], [60; 60], '-w', 'LineWidth', 5)
text(55,57, '11$\mu m$', 'HorizontalAlignment','center','Interpreter','latex','FontWeight','bold','FontSize',25,'Color','w')
%saveFigMF(gcf, '~/source/latex/bimsim/figures/source/m_pmma_a2_2_wav_6_9','m','fig')

figure; imagesc(t_p),axis image, colormap(brewer),colorbar, axis off
%xlabel(['\lambda = ', num2str(1e4/wav(idx))],'Interpreter','latex')
set(gca,'FontSize',14,'fontweight','bold')
hold on
plot([60; 50], [60; 60], '-w', 'LineWidth', 5)
text(55,57, '11$\mu m$', 'HorizontalAlignment','center','Interpreter','latex','FontWeight','bold','FontSize',25,'Color','w')
%saveFigMF(gcf, '~/source/latex/bimsim/figures/source/p_pmma_a1_5_wav_6_8','m','fig')

%% crop and save absorbance images
cr = 64; sz = 128;
idx = 25;
t_m = absImages(:,:,idx);
%t_m(t_m<0) = 0;
t_m = t_m(round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1,round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1);
t_p = A_b(:,:,idx);
%t_p(t_p<0) = 0;
t_p = t_p(round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1,round(sz/2)- round(cr/2):round(sz/2) + round(cr/2) - 1);
figure; imagesc(t_m),axis image, colormap(brewer),colorbar, axis off
%xlabel(['{\lambda} = ', num2str(1e4/wav(idx))],'Interpreter','latex')
set(gca,'FontSize',14,'fontweight','bold')
hold on
plot([60; 50], [60; 60], '-w', 'LineWidth', 5)
text(55,57, '11$\mu m$', 'HorizontalAlignment','center','Interpreter','latex','FontWeight','bold','FontSize',25,'Color','w')
%saveFigMF(gcf, '~/source/latex/bimsim/figures/source/m_pmma_a2_2_wav_6_9','m','fig')

figure; imagesc(t_p),axis image, colormap(brewer),colorbar, axis off
%xlabel(['\lambda = ', num2str(1e4/wav(idx))],'Interpreter','latex')
set(gca,'FontSize',14,'fontweight','bold')
hold on
plot([60; 50], [60; 60], '-w', 'LineWidth', 5)
text(55,57, '11$\mu m$', 'HorizontalAlignment','center','Interpreter','latex','FontWeight','bold','FontSize',25,'Color','w')

%% plot cluster of spectra for pmma of d = 4.52+-0.15
[wav, absImages, spec] = loadFTIR('/home/sberisha/data/ftir/polystyrene/polystyrene_4_52/polystyrene_4_52.hdr', '/home/sberisha/data/ftir/polystyrene/polystyrene_4_52/polystyrene_4_52',3,[],[]);

[wav_b, spec_b, A_b, int_b, inc_b] = loadCudaBimSim([dataGenCudaPath 'ps_a2_2_fo_32_mc1024_highMag_fpAta_ps76_141_nBP_allWav_s5_n1_4/'],'absSpec', 'out_a','out_i', 'out_inc', 128,128, 2);

absImages = absImages(:,:,2:2:end);
wav = wav(2:2:end);
spec = spec(2:2:end);


spec = squeeze(absImages(65,63,:));
figure
p1 = plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(65,64,:));
hold on
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(65,65,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(66,63,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(66,64,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(66,65,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(67,63,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(67,64,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(67,65,:));
plot(wav,spec, 'b', 'LineWidth', 2)

spec_b = squeeze(A_b(64,64,:));
p2 = plot(wav_b,spec_b, 'r--','LineWidth',2)
hold on
spec_b = squeeze(A_b(64,65,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(64,66,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(65,64,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(65,65,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(65,66,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(66,64,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(66,65,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(66,66,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)

xlim(gca, [1000 3500])
ylim(gca,[0 1])
xlabel('\textbf{Wavenumber ($\mathbf{cm^{-1}}$)}', 'Interpreter', 'latex') 
ylabel(gca,'Absorbance','fontweight','bold') 
set(gca,'FontSize',14,'fontweight','bold')
h_legend = legend([p1 p2],'measured', 'predicted', 'Location', 'Best');
set(h_legend,'FontSize',14, 'fontweight','bold')
grid on



%% plot results of pmma d = spec_128ca_8r_sphere1, simulation a = 2.2 -- ps with d = 4.52+-0.15

[wav, absImages, spec] = loadFTIR('/home/sberisha/data/ftir/polystyrene/polystyrene_5_43/polystyrene_5_43.hdr', '/home/sberisha/data/ftir/polystyrene/polystyrene_5_43/polystyrene_5_43.dat',3,[],[]);

spec = squeeze(absImages(62,65,:));

[wav_b, spec_b, A_b, int_b, inc_b] = loadCudaBimSim([dataGenCudaPath 'ps_a2_5_fo_32_mc1024_highMag_fpAta_ps88_141_nBP_allWav_s5_n1_45/'],'absSpec', 'out_a','out_i', 'out_inc', 128,128, 2);


absImages = absImages(:,:,2:2:end);
wav = wav(2:2:end);
spec = spec(2:2:end);
figure
plot(wav,spec, 'LineWidth',2); hold on; plot(wav_b,spec_b, 'LineWidth',2)
xlim(gca, [1000 3500])
ylim(gca,[min(spec) 1])
xlabel('\textbf{Wavenumber ($\mathbf{cm^{-1}}$)}', 'Interpreter', 'latex') 
ylabel(gca,'Absorbance','fontweight','bold') 
set(gca,'FontSize',14,'fontweight','bold')
h_legend = legend('measured', 'predicted', 'Location', 'Best');
set(h_legend,'FontSize',14, 'fontweight','bold')
grid on


%% plot cluster of spectra for pmma of d = 5.43
[wav, absImages, spec] = loadFTIR('/home/sberisha/data/ftir/polystyrene/polystyrene_5_43/polystyrene_5_43.hdr', '/home/sberisha/data/ftir/polystyrene/polystyrene_5_43/polystyrene_5_43.dat',3,[],[]);

[wav_b, spec_b, A_b, int_b, inc_b] = loadCudaBimSim([dataGenCudaPath 'ps_a2_5_fo_32_mc1024_highMag_fpAta_ps88_141_nBP_allWav_s5_n1_45/'],'absSpec', 'out_a','out_i', 'out_inc', 128,128, 2);

absImages = absImages(:,:,2:2:end);
wav = wav(2:2:end);
spec = spec(2:2:end);


spec = squeeze(absImages(61,64,:));
figure
p1 = plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(61,65,:));
hold on
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(61,66,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(62,64,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(62,65,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(62,66,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(63,64,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(63,65,:));
plot(wav,spec, 'b', 'LineWidth', 2)
spec = squeeze(absImages(63,66,:));
plot(wav,spec, 'b', 'LineWidth', 2)

spec_b = squeeze(A_b(64,64,:));
p2 = plot(wav_b,spec_b, 'r--','LineWidth',2)
hold on
spec_b = squeeze(A_b(64,65,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(64,66,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(65,64,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(65,65,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(65,66,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(66,64,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(66,65,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)
spec_b = squeeze(A_b(66,66,:));
plot(wav_b,spec_b, 'r--','LineWidth',2)

xlim(gca, [1000 3500])
ylim(gca,[0 1])
xlabel('\textbf{Wavenumber ($\mathbf{cm^{-1}}$)}', 'Interpreter', 'latex') 
ylabel(gca,'Absorbance','fontweight','bold') 
set(gca,'FontSize',14,'fontweight','bold')
h_legend = legend([p1 p2],'measured', 'predicted', 'Location', 'Best');
set(h_legend,'FontSize',14, 'fontweight','bold')
grid on


