[bic,waxis] = bicoherence(lfp(1, :), 'nsamp', 1024);

freq = waxis.*1000;
[c, f1] = min(abs(freq - 0));
[c, f2] = min(abs(freq - 150));
bico = bic(f1:f2, f1:f2);
freq = freq(f1:f2);


%%
bico = triu(bico); bico = flip(bico);
bico = tril(bico); bico = flip(bico);
bico(bico == 0) = NaN;

contourf(freq, freq, bico, 100,'linecolor','none')

hold on
axis xy
Reds=cbrewer('seq', 'Reds', 201); %caxis([0 0.03]);
colormap(Reds);
cb = colorbar;
set(cb,'Position',[0.85 0.20 0.02 0.3])
xlabel('Frequency[Hz]'); ylabel('Frequency[Hz]')
ylim([0 75])
pbaspect([1 1 1])
set(gca, 'box', 'off', 'tickdir', 'out','LineWidth', 1.5)

set(gca,'ytick', [10 30 50 75])
set(gca,'xtick', [10 30 50 100 150])
