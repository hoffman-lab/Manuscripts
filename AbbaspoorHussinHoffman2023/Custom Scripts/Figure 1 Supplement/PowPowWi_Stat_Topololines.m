clear all, clc

load('D:\HiTGun\Wi TH Data\Analyses\Cross-Freq Power Correlation\PowPowPvalueTableWi_withStats_log.mat')

PowPowWi = PowPow(2:end, [1 4 5 6]);

load('D:\HiTGun\LU TH Data\Analyses\Cross-Freq Power Correlation\PowPow_uV_fLU.mat')

PowPow = [PowPow; PowPowWi];

%% Logscale plotting

fh = figure('Renderer', 'painters')


Search = 0;
for ev = 1:size(PowPow, 1)
    Search = Search + PowPow.SearchCorr{ev};
end
Search = Search/size(PowPow, 1);

Prop = 0;
for ev = 1:size(PowPow, 1)
    Prop = Prop + PowPow.Significance{ev};
end
Prop = Prop/size(PowPow, 1);

freqs = logspace(log10(1), log10(150), 400) %linspace(1, 150, 400); %logspace(log10(1), log10(200), nff);
contourf(freqs, freqs, Search, 101,'linecolor','none')
hold on

K = (1/16)*ones(4);
Prop = conv2(Prop,K,'same');

Prop(Prop<0.25) = 0;
Prop(Prop>=0.25) = 1;
contour(freqs, freqs, Prop, 3, 'LineColor', 'k')


axis xy
RdBu=cbrewer('div', 'RdBu', 101);
colormap(flip(RdBu)); caxis([-0.15 0.15])
cb = colorbar; set(cb,'Position',[0.87 0.3 0.02 0.3])
xlabel('Log Frequency'); ylabel('Log Frequency')
% title('Cross-freq power Correlation')
pbaspect([1 1 1])
set(gca, 'box', 'off', 'tickdir', 'out','LineWidth', 1.5)
set(gca,'yscale','log', 'xscale', 'log')

set(gca,'ytick', [5 10 30 50 100 150])
set(gca,'xtick', [5 10 30 50 100 150])

FigureDirectory = 'D:\HiTGun\0 Main\Figures - refined'
FigureName = 'WiLUPowPowCorrelation_topolineplot'
print(fh, fullfile(FigureDirectory, FigureName),'-dsvg','-r300')
print(fh, fullfile(FigureDirectory, FigureName),'-dpng','-r300')
savefig(fh, fullfile(FigureDirectory, FigureName))

close(fh)
