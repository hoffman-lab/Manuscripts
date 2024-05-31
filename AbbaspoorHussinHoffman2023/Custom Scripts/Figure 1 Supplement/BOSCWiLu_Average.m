clear all, clc
Directory = 'D:\HiTGun\Wi TH Data\Analyses\BOSC';
load(fullfile(Directory, 'BOSC_Wi_RestSearch'))

Directory = 'D:\HiTGun\LU TH Data\Analyses\BOSC';
load(fullfile(Directory, 'BOSC_LU_RestSearch'))

BOSC_Table = [WIBOSC; LUBOSC];

%%
Occupancy_search = cell2mat(BOSC_Table.Search)*100;
Occupancy_rest = cell2mat(BOSC_Table.Rest)*100;

BootSearch = bootstrp(10000, @mean, Occupancy_search);
BootRest = bootstrp(10000, @mean, Occupancy_rest);

CISearch = prctile(BootSearch,[5,95]);
CIRest   = prctile(BootRest,  [5,95]);


Freq2Test = [1:30 35 40 50 60]; Freq2Test=unique(Freq2Test);
pvalue = NaN(1, length(Freq2Test));
for ii = 1:length(Freq2Test)
    Idx = knnsearch(freqs',Freq2Test(ii));
    [pvalue(ii),h,stats] = signrank(Occupancy_search(:, Idx), Occupancy_rest(:, Idx));
end

[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvalue, 0.05);

Significant = NaN(size(freqs));
Significant(Freq2Test) = h;
Significant(Significant<1) = NaN;


fh = figure('Renderer', 'painters', 'Position', [50 50 800 800])


Significant = find(Significant == 1);
for ii = 1:length(Significant)
    fill([Significant(ii)-0.5 Significant(ii)-0.5 Significant(ii)+0.5 Significant(ii)+0.5], [4 16 16 4], [0.9 0.9 0.9], 'edgecolor', [0.9 0.9 0.9])
    hold on
end


hold on
Col = [255 8 32]/255;
plot(freqs,smooth(mean(Occupancy_search)), 'Color', [0 0 0], 'LineWidth', 2)
h1 = fill([freqs flip(freqs)], [smooth(CISearch(1, :)); smooth(flip(CISearch(2, :)))], Col, 'FaceAlpha', 0.7, 'EdgeColor', Col)

Col = [56 138 191]/255; %[30 105 188]  [140 192 220]
plot(freqs,smooth(mean(Occupancy_rest)), 'Color', [0 0 0], 'LineWidth', 2)
h2 = fill([freqs flip(freqs)], [smooth(CIRest(1, :)); smooth(flip(CIRest(2, :)))], Col, 'FaceAlpha', 0.7, 'EdgeColor', Col)
legend([h1 h2], {'Search', 'Rest'})
legend boxoff


set(gca, 'LineWidth', 3, 'Box', 'off', 'Tickdir', 'out', 'FontSize', 15,...
    'xtick', 0:10:70)
ylabel('Occupancy Rate (%)')
xlabel('Frequency (Hz)')
xlim([0 70])

pbaspect([1 1 1])


FigureDirectory = 'D:\HiTGun\Wi TH Data\Analyses'
FigureName = 'AverageBOSCWiLU'
print(fh, fullfile(FigureDirectory, FigureName),'-dsvg','-r300')
print(fh, fullfile(FigureDirectory, FigureName),'-dpng','-r300')
% close(fh)


