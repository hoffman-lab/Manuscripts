%% Load Data
filedir = 'D:\HiTGun\Wi TH Data\Processed';
Sessions = struct2table(dir(filedir));
Sessions = Sessions(3:end, :);

WIBicoherence = table;
for sess = 1:size(Sessions, 1)
    
    Session = Sessions.name{sess};
    
    %% Compute Power Spectrum - REST
    % load(fullfile(filedir, Sessions.name{sess}, 'Sleep'))
    % [bico,waxis] = bicoherence(lfp, 'nsamp', 1024);
    %
    % freq = waxis.*1000;
    % [c, f1] = min(abs(freq - 0));
    % [c, f2] = min(abs(freq - 150));
    % bico = bico(f1:f2, f1:f2);
    % freq = freq(f1:f2);
    %
    % Bicoherence.session(sess) = {Session};
    % Bicoherence.freq(sess) = {freq};
    % Bicoherence.BicREst(sess) = {bico};
    
    
    %% Compute Power Spectrum - Search
    load(fullfile(filedir, Sessions.name{sess}, 'Treehouse'))
    to_remove = cellfun(@(x) isempty(x), TrialLFP);
    TrialLFP(to_remove) = []; TrialTime(to_remove) = [];
    
    BICOHERENCE = cell(length(TrialLFP), 1);
    for Trial = 1:length(TrialLFP)
        [bic,waxis] = bicoherence(TrialLFP{1, Trial}, 'nsamp', 1024);
        BICOHERENCE{Trial} = bic;
    end
    
    bico = cat(3,BICOHERENCE{:});
    bico = mean(bico, 3);
    
    freq = waxis.*1000;
    [c, f1] = min(abs(freq - 0));
    [c, f2] = min(abs(freq - 150));
    bico = bico(f1:f2, f1:f2);
    freq = freq(f1:f2);
    
    %% Statistical Test
    
    if sess == 1
        load('D:\HiTGun\Bicoherence_montecarlo_cluster.mat')
    end
    
    % now threshold real data...
    % first Z-score
    zmap = (bico-mean_h0) ./ std_h0;
    
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
    
    %% Plot
    fh = figure();
    fh.WindowState = 'maximized';
        
    bico = triu(bico); bico = flip(bico);
    bico = tril(bico); bico = flip(bico);
    
    zmap = triu(zmap); zmap = flip(zmap);
    zmap = tril(zmap); zmap = flip(zmap);
    
    subplot(121)
    imagesc(freq, freq, bico)
    %         hold on
    %         contour(freq,freq,zmap,[1 1],'k','linewidth',5);
    axis xy; xlabel('Frequency'); ylabel('Frequency')
    Reds=cbrewer('seq', 'Reds', 201); caxis([0 0.1]);
    colormap(Reds);
    cb = colorbar; set(cb,'Position',[0.48 0.37 0.0083 0.1974])
    set(gca, 'box', 'off', 'FontSize', 15, 'tickdir', 'out', 'LineWidth', 1.5,...
        'xtick', [8 16 32 50 100 150], 'ytick', [8 16 32 50 75])
    ylim([0 75])
    pbaspect([1 1 1])
    
    
    subplot(122)
    imagesc(freq, freq, zmap)
    axis xy; xlabel('Frequency'); ylabel('Frequency')
    Reds=cbrewer('seq', 'Reds', 201); %caxis([0 0.03]);
    colormap(Reds);
    cb = colorbar; set(cb,'Position',[0.92 0.37 0.0083 0.1974])
    set(gca, 'box', 'off', 'FontSize', 15, 'tickdir', 'out', 'LineWidth', 1.5,...
        'xtick', [8 16 32 50 100 150], 'ytick', [8 16 32 50 75])
    ylim([0 75])
    pbaspect([1 1 1])
    
    %%
    figName = [Session]
    Directory = 'D:\HiTGun\Wi TH Data\Analyses\Bicoherence'
    saveas(gcf, fullfile(Directory, [figName '.png']))
    close(gcf)
    
    %%
    
    WIBicoherence.session(sess) = {Session};
    WIBicoherence.freq(sess) = {freq};
    WIBicoherence.Bic(sess) = {bico};
    WIBicoherence.significance(sess) = {zmap};
    
end

save(fullfile(Directory, 'WIBicoherence'), 'WIBicoherence', 'freq', '-v7.3')


%%

bico = 0;
for ev = 1:size(WIBicoherence, 1)
    if ev == 5
        continue
    end
    bico = bico + WIBicoherence.Bic{ev};
end
bico = bico/size(WIBicoherence, 1);



fh = figure();
fh.WindowState = 'maximized';

bico = WIBicoherence.BicSearch{1,1};
bico = bic;
freq = waxis.*1000;
[c, f1] = min(abs(freq - 0));
[c, f2] = min(abs(freq - 150));
bico = bico(f1:f2, f1:f2);
freq = freq(f1:f2);

bico = triu(bico); bico = flip(bico);
bico = tril(bico); bico = flip(bico);

bico(bico==0) = NaN;

contourf(freq, freq, bico, 100,'linecolor','none')

hold on
axis xy
Reds=cbrewer('seq', 'Reds', 201); %caxis([0 0.1]);
colormap(Reds);
cb = colorbar;
set(cb,'Position',[0.85 0.20 0.02 0.3])
xlabel('Frequency[Hz]'); ylabel('Frequency[Hz]')
ylim([0 75])
pbaspect([1 1 1])
set(gca, 'box', 'off', 'tickdir', 'out','LineWidth', 1.5)

set(gca,'ytick', [10 30 50 75])
set(gca,'xtick', [10 30 50 100 150])



