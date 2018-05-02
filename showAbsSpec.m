function showAbsSpec(wav,spec)

h=axes('FontSize',14);
plot(wav,spec,'LineWidth',2)
set(h, 'Xdir','reverse')