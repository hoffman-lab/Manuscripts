filedir = 'D:\HiTGun\Wi TH Data\Analyses\Power Spectral Density';
load(fullfile(filedir, 'PowerSpectralDensity_raw.mat'))

RestPSD   = cell2mat(PowerSpectralDensity.Rest);
SearchPSD = cell2mat(PowerSpectralDensity.Search);


filedir = 'D:\HiTGun\LU TH Data\Analyses\Power Spectral Density';
load(fullfile(filedir, 'PowerSpectralDensity_raw.mat'))


RestPSD   = [RestPSD; cell2mat(PowerSpectralDensity.Rest)];
SearchPSD = [SearchPSD; cell2mat(PowerSpectralDensity.Search)];



%%

RestPSD = zscore(pow2db(RestPSD), [], 2);
SearchPSD = zscore(pow2db(SearchPSD), [], 2);


BootSearch = bootstrp(10000, @mean, SearchPSD);
BootRest   = bootstrp(10000, @mean, RestPSD);

CISearch   = prctile(BootSearch,[5,95]);
CIRest     = prctile(BootRest,  [5,95]);

%%

fh = figure('position', [50 50 700 700], 'Renderer', 'painters')

hold on
Col = [255 8 32]/255;
h1 = plot(freqs,smooth(mean(BootSearch)), 'Color', [0 0 0], 'LineWidth', 2)
fill([freqs flip(freqs)], [smooth(CISearch(1, :)); smooth(flip(CISearch(2, :)))], Col, 'FaceAlpha', 0.7, 'EdgeColor', Col)

Col = [56 138 191]/255; %[30 105 188]  [140 192 220]
h2 = plot(freqs,smooth(mean(BootRest)), 'Color', [0 0 0], 'LineWidth', 2)
fill([freqs flip(freqs)], [smooth(CIRest(1, :)); smooth(flip(CIRest(2, :)))], Col, 'FaceAlpha', 0.7, 'EdgeColor', Col)

% set(gca, 'box', 'off', 'tickdir', 'out','LineWidth', 1.5, 'xscale', 'log')
xtickangle(30)

ylabel('Power (dB)')
xlabel('Frequency (Hz)')
set(gca, 'box', 'off', 'tickdir', 'out', 'LineWidth', 1.5)

%%
% subplot(6, 1, 6)
ax3 = axes('Position', [0.15 0.1 0.8 0.12])

area(freqs,smooth(mean(BootSearch))-smooth(mean(BootRest)), 'FaceColor', [0 0 0], 'LineWidth', 3)
set(gca, 'box', 'off', 'tickdir', 'out', 'LineWidth', 1.5, 'xscale', 'log')
set(gca,'xtick',[5 10 15 20 30 50 100 150])
xlabel('Log Frequency'); ylabel('Search-Rest')
ylim([-0.15 0.15])
% pbaspect([1 0.184 1])


