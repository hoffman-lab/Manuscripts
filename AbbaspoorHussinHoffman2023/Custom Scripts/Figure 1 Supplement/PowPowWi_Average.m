Rest = NaN(15, 400, 400);
for ii = 1:15
    Rest(ii, :, :) = PowPow.RestCorr{ii};
end
Rest = squeeze(mean(Rest, 1));

fh = figure('Renderer', 'painters', 'Position', [50 50 800 800])
Corr = Rest;
contourf(freqs, freqs, Corr, 101,'linecolor','none')
axis xy
RdBu=cbrewer('div', 'RdBu', 101);
colormap(flip(RdBu)); caxis([-0.15 0.15])
cb = colorbar; set(cb,'Position',[0.92 0.3 0.02 0.3])
xlabel('Log Frequency'); ylabel('Log Frequency')
% title('Cross-freq power Correlation')
pbaspect([1 1 1])
set(gca, 'box', 'off', 'tickdir', 'out','LineWidth', 1.5)
set(gca,'yscale','log', 'xscale', 'log')

set(gca,'ytick', [5 10 30 50 100 150])
set(gca,'xtick', [5 10 30 50 100 150])



%%
Search = NaN(6, 400, 400);
for ii = 10:15
    Search(ii-9, :, :) = PowPow.SearchCorr{ii};
end
Search = squeeze(mean(Search, 1));

fh = figure('Renderer', 'painters', 'Position', [50 50 800 800])
Corr = Search;
contourf(freqs, freqs, Corr, 101,'linecolor','none')
axis xy
RdBu=cbrewer('div', 'RdBu', 101);
colormap(flip(RdBu)); caxis([-0.15 0.15])
cb = colorbar; set(cb,'Position',[0.92 0.3 0.02 0.3])
xlabel('Log Frequency'); ylabel('Log Frequency')
% title('Cross-freq power Correlation')
pbaspect([1 1 1])
set(gca, 'box', 'off', 'tickdir', 'out','LineWidth', 1.5)
set(gca,'yscale','log', 'xscale', 'log')

set(gca,'ytick', [5 10 30 50 100 150])
set(gca,'xtick', [5 10 30 50 100 150])

