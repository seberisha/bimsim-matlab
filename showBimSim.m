
function showBimSim(p, E_f, E_s, E_i,E_t, Et_d, Ef_d, A)

brewer = brewermap(1000);
colormap(brewer)

figure, 
subplot(3,3,1)
imagesc((abs((E_f)))),title(sprintf('abs(E_f)')), colorbar, axis image
subplot(3,3,2)
imagesc((abs((E_s)))),title(sprintf('abs(E_s)')), colorbar, axis image
subplot(3,3,3)
imagesc((abs((E_i)))),title(sprintf('abs(E_i)')), colorbar, axis image
subplot(3,3,4)
imagesc((abs((E_t)))),title(sprintf('abs(E_t)')), colorbar, axis image
subplot(3,3,5)
imagesc(abs(Et_d)), title(sprintf('D_{Et}')),axis image, colormap(brewer), colorbar
subplot(3,3,6)
imagesc(abs(Ef_d)), title(sprintf('D_{Ef}')),axis image, colormap(brewer), colorbar
subplot(3,3,7)
imagesc((A)), title(sprintf('A ')),axis image, colormap(brewer), colorbar
%suptitle(sprintf('Results for point source %i',p))

axes('Position',[0 0 1 1],'Visible','off');
text(0.45,0.97,sprintf('Results for point source %i',p))