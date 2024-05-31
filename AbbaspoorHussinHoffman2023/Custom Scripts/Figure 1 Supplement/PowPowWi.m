%% Load Data
filedir = 'D:\HiTGun\Wi TH Data\Processed';
Sessions = struct2table(dir(filedir));
Sessions = Sessions(3:end, :);

PowPow = table;
for sess = 2:size(Sessions, 1)
    
    Session = Sessions.name{sess};
    PowPow.session(sess) = {Session};
    
%     %% Compute Power Spectrum - REST
%     ov       = 0.5;   % overlap
%     winlen   = 1024; % winlen
%     nff      = 400;
%     freqBins = logspace(log10(1), log10(150), nff);  %linspace(1, 150, nff); %
%     
%     PSD = cell(1, 2);
%     load(fullfile(filedir, Sessions.name{sess}, 'Sleep'))
%     [s,freqs,time, pd] = spectrogram(lfp,winlen,floor(ov*winlen), freqBins, 1000);
%     
%     PowPow.Rest(sess) = {pd};
%     PPCor = corrcoef(pd');
%     PowPow.RestCorr(sess) = {PPCor};
    
    
    %% Compute Power Spectrum - Search
    load(fullfile(filedir, Sessions.name{sess}, 'Treehouse'))
    to_remove = cellfun(@(x) isempty(x), TrialLFP);
    TrialLFP(to_remove) = []; TrialTime(to_remove) = [];
    
    PSD = cell(1, length(TrialLFP));
    for Trial = 1:length(TrialLFP)
        [s,freqs,time, pd] = spectrogram(TrialLFP{1, Trial},winlen,floor(ov*winlen), freqBins, 1000);
        PSD{Trial} = pd;
    end
    
    PSD = {cat(2, PSD{:})}; PSD = PSD{1};
    PowPow.Search(sess) = {PSD};
    PPCor = corrcoef(PSD');
    PowPow.SearchCorr(sess) = {PPCor};
    
    %% Statistical Significance
    n_permutes = 5000; pval = 0.05;
    
    [permmaps, mean_h0, std_h0, zval, cluster_thresh] = ...
        PowPowCor_Permutation_ClusterBased(PSD, n_permutes, pval);
    
    % now threshold real data...
    % first Z-score
    zmap = (PPCor-mean_h0) ./ std_h0;
    
    % threshold image at p-value, by setting subthreshold values to 0
    zmap(abs(zmap)<zval) = 0;
    
    
    % now find clusters in the real thresholded zmap
    % if they are "too small" set them to zero
    islands = bwconncomp(zmap);
    for i=1:islands.NumObjects
        % if real clusters are too small, remove them by setting to zero!
        if numel(islands.PixelIdxList{i}==i)<cluster_thresh
            zmap(islands.PixelIdxList{i})=0;
        end
    end
    
    %%
    fh = figure();
    fh.WindowState = 'maximized';
    
    subplot(121)
    contourf(freqs, freqs, PPCor, 100,'linecolor','none')
    axis xy
    RdBu=cbrewer('div', 'RdBu', 101);
    colormap(flip(RdBu)); caxis([-0.15 0.15])
    cb = colorbar;
    xlabel('Frequency'); ylabel('Frequency')
    pbaspect([1 1 1])
    set(gca, 'box', 'off', 'tickdir', 'out','LineWidth', 1.5)
    set(gca,'yscale','log', 'xscale', 'log')
    
    set(gca,'ytick', [5 10 30 50 100 150])
    set(gca,'xtick', [5 10 30 50 100 150])
    
    subplot(122)
    zmap(abs(zmap) ~= 0) = 1;
    contourf(freqs, freqs, zmap, 100,'linecolor','none')
    axis xy
    caxis([0 1])
    cb = colorbar;
    xlabel('Frequency'); ylabel('Frequency')
    pbaspect([1 1 1])
    set(gca, 'box', 'off', 'tickdir', 'out','LineWidth', 1.5)
    set(gca,'yscale','log', 'xscale', 'log')
    
    set(gca,'ytick', [5 10 30 50 100 150])
    set(gca,'xtick', [5 10 30 50 100 150])
    
    title('LE Cross-freq power Correlation')
    
    %%
    figName = [Session, '_Rest'];
    Directory = 'D:\HiTGun\Wi TH Data\Analyses\Cross-Freq Power Correlation'
    saveas(gcf, fullfile(Directory, [figName '.png']))
    close(gcf)
    
    PowPow.Significance(sess) = {zmap};
    
end

WiPSD = PowPow;
save(fullfile(Directory, 'PowPow_uV'), 'PowPow', 'freqs', '-v7.3')



%%

PPCor = NaN(14, 400, 400);
for ii = 2:15
    PPCor(ii-1, :, :) = PowPow.SearchCorr{ii,1};
end
PPCor = squeeze(mean(PPCor, 1));

fh = figure();
fh.WindowState = 'maximized';

subplot(121)
imagesc(freqs, freqs, PPCor)
axis xy
RdBu=cbrewer('div', 'RdBu', 101);
colormap(flip(RdBu)); caxis([-0.1 0.1])
cb = colorbar;
xlabel('Frequency'); ylabel('Frequency')
pbaspect([1 1 1])
set(gca, 'box', 'off', 'tickdir', 'out','LineWidth', 1.5)
set(gca,'yscale','log', 'xscale', 'log')

set(gca,'ytick', [5 10 30 50 100 150])
set(gca,'xtick', [5 10 30 50 100 150])


    