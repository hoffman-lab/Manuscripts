cd('D:\HiTGun\Data\fooof\Session Wise Mean')

Rest = readNPY('Rest.npy');
Search = readNPY('Search.npy');
load('freqs.mat')

%% Shuffling Test
% 3-10 (@2 Hz res); 10-40 (@5 Hz) and 40-90 (@10 Hz)

Freq2Test = [3:2:10, 10:5:40 50:10:90];

pvalue = NaN(1, length(Freq2Test));
for ii = 1:length(Freq2Test)
    Idx = knnsearch(freqs',Freq2Test(ii));
    [pvalue(ii),h,stats] = signrank(Search(:, Idx), Rest(:, Idx));
end

[h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvalue, 0.05,'dep');

%%

BootSearch = bootstrp(10000, @mean, Search);
BootRest = bootstrp(10000, @mean, Rest);

CISearch = prctile(BootSearch,[5,95]); 
CIRest   = prctile(BootRest,  [5,95]);


%% Plotting


fh = figure('Renderer', 'painters')

subplot(5, 1, 2:5)

for ii = [7:10,15:40,40:70]
%     fill([ii-0.5 ii-0.5 ii+0.5 ii+0.5], [-0.295 0.6 0.6 -0.295], [0.9 0.9 0.9], 'edgecolor', [0.9 0.9 0.9])
     plot([ii-0.5 ii+0.5], [0.3 0.3], 'Color',[0 0 0], 'LineWidth', 3)
  
hold on
end

hold on
Col = [95 194 192]/255;
h1 = plot(freqs,smooth(mean(BootSearch)), 'Color', Col, 'LineWidth', 2)
fill([freqs flip(freqs)], [smooth(CISearch(1, :)); smooth(flip(CISearch(2, :)))], Col, 'FaceAlpha', 0.7, 'EdgeColor', Col)

Col = [88 97 121]/255;
h2 = plot(freqs,smooth(mean(BootRest)), 'Color', Col, 'LineWidth', 2)
fill([freqs flip(freqs)], [smooth(CIRest(1, :)); smooth(flip(CIRest(2, :)))], Col, 'FaceAlpha', 0.7, 'EdgeColor', Col)

set(gca, 'box', 'off', 'tickdir', 'out','LineWidth', 1.5, 'xscale', 'log')
set(gca,'xtick',[5 10 15 20 30 50 100 150])
xtickangle(30)

xlabel('Log Frequency'); ylabel('Power')
set(gca, 'box', 'off', 'tickdir', 'out', 'LineWidth', 1.5, 'xscale', 'log')

legend([h1 h2], {'During Search', 'During Rest'})
legend boxoff


subplot(5, 1, 1)
area(freqs,smooth(mean(BootSearch))-smooth(mean(BootRest)), 'FaceColor', [0 0 0], 'LineWidth', 3)
ylabel('Search-Rest')
set(gca, 'box', 'off', 'tickdir', 'out', 'LineWidth', 1.5, 'xscale', 'log', 'xtick', [], 'Xcolor', [1 1 1])
ylim([-0.15 0.15])
% pbaspect([1 0.184 1])




%% Average FFT incet

load('D:\HiTGun\Data\LE_Pwelch_All_Sessions');
load('D:\HiTGun\Data\LU_PwelchLU_All_Sessions');

Pwelch = [Pwelch; PwelchLU];


%%

PowerSpectrumSearch = [];
for Sess = 1:size(Pwelch)
    PowerSpectrumSearch = [PowerSpectrumSearch; Pwelch.Search{Sess,1}];
end

Search = mean(PowerSpectrumSearch);


PowerSpectrumRest = [];
for Sess = 1:size(Pwelch)
    PowerSpectrumRest = [PowerSpectrumRest; Pwelch.Rest{Sess,1}];
end

PowerSpectrumDiff = [];
for Sess = 1:size(Pwelch)
    Diff = mean(pow2db(Pwelch.Search{Sess,1})) - mean(pow2db(Pwelch.Rest{Sess,1}));
    PowerSpectrumDiff = [PowerSpectrumDiff; Diff];
end

Rest = mean(PowerSpectrumRest);


axes('Position',[0.25 .52 .25 .2])
% Col = [95 194 192]/255;
Col = [88 97 121]/255;
h1 = plot(freqs,mean(pow2db(PowerSpectrumSearch)), 'Color', Col, 'LineWidth', 1)
hold on

% Col = [33 91 167]/255;
Col = [95 194 192]/255;
h2 = plot(freqs,mean(pow2db(PowerSpectrumRest)), 'Color', Col, 'LineWidth', 1)

set(gca, 'box', 'off', 'tickdir', 'out', 'LineWidth', 1, 'xscale', 'log')
set(gca,'xtick',[5 15 30 50 100 150])
xtickangle(90)

% pbaspect([1 1 1])

ylabel('10log10[power]')


%%

% FigureDirectory = 'D:\HiTGun\Figures\Grand'
% FigureName = 'AverageParametrizedFFT'
% print(fh, fullfile(FigureDirectory, FigureName),'-dsvg','-r300')
% print(fh, fullfile(FigureDirectory, FigureName),'-dpng','-r300')
% 
% close(fh)

