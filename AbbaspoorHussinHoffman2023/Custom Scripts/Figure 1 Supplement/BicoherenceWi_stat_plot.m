clear all, clc

load('D:\HiTGun\Wi TH Data\Analyses\Bicoherence\WIBicoherence.mat')

Bicoherence = WIBicoherence;
Bico = NaN(15, 155, 155);
for Sess = 1:size(Bicoherence, 1)
    Bico(Sess, :, :) = Bicoherence.Bic{Sess,1};
end
Bico = squeeze(mean(Bico, 1));


Sig = Bicoherence.significance;

Significance = NaN(155, 155, size(Sig, 1));
for i = 1:15
    tmp = Sig{i, 1}; tmp(tmp > 0) = 1;
    Significance(:, :, i) = tmp;
end

Significance = mean(Significance, 3);


%%  plotting

fh = figure('Renderer', 'painters')

freqs = freq;
Bico = triu(Bico); Bico = flip(Bico);
Bico = tril(Bico); Bico = flip(Bico);
Bico(Bico == 0) = NaN;

contourf(freqs, freqs, Bico, 100,'linecolor','none')
hold on
Prop = Significance;
Prop(Prop<0.6) = 0;
Prop(Prop>0.6) = 1;
Prop(isnan(Bico)) = NaN;
contour(freqs, freqs, Prop, 1, 'LineColor', 'k')

axis xy
Reds=cbrewer('seq', 'Reds', 201); caxis([0 0.03]);
colormap(Reds);
cb = colorbar;
set(cb,'Position',[0.85 0.3 0.02 0.3]) 
xlabel('Frequency(Hz)'); ylabel('Frequency(Hz)')
ylim([0 75])
pbaspect([1 1 1])
set(gca, 'box', 'off', 'tickdir', 'out','LineWidth', 1.5)

set(gca,'ytick', [10 30 50 75])
set(gca,'xtick', [10 30 50 100 150])


FigureDirectory = 'D:\HiTGun\0 Main\Figures - refined'
FigureName = 'WiBicoherenceTopolineplot'
print(fh, fullfile(FigureDirectory, FigureName),'-dsvg','-r300')
print(fh, fullfile(FigureDirectory, FigureName),'-dpng','-r300')
savefig(fh, fullfile(FigureDirectory, FigureName))

close(fh)
