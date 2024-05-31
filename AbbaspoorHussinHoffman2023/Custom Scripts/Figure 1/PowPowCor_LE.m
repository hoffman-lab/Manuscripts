clear all, clc
addpath 'D:\Matlab Packages\fieldtrip-lite-20190224'


FilesDir = 'D:\HiTGun\Data\decimatedCSC_wholeSession';
TTLFilesDir = 'D:\HiTGun\Data\LE Change blindness Raw';

load('D:\HiTGun\Data\LETrialsTable');

info = struct2table(dir(FilesDir));
SessionTable = info.name;
SessionTable = SessionTable(3:end);

PowPow = table;
Row = 1;
for Sess = 1:size(SessionTable, 1)
    
    if Sess == 18 | Sess == 19; continue; end
    
    %% Load Data
    Session = SessionTable{Sess};
    
    %%% Find in-layer-channel
    % This line loads the LE_In_Layer_Channel table that contains the channel
    % of recording that SWRs were found on. This is used as a marker to
    % selecting in-layer-channel.
    load(fullfile('D:\HiTGun\Data', 'LE_In_Layer_Channel.mat'));
    Index = find(contains(string(LE_In_Layer_Channel.sid),Session));
    if isempty(Index); continue; end
    
    inchannel = LE_In_Layer_Channel.rip_cscNum(Index);
    channellist = [inchannel 13 14 15];
    channellist = [inchannel]; %unique(channellist);
    
    for ch = 1:length(channellist)
        
        channel = {['CSC', num2str(channellist(ch))]};
        
        %%% Load Channel LFP
        filedir = fullfile(FilesDir, Session, channel);
        load(filedir{1})
        
        
        filedir = fullfile(TTLFilesDir, Session);
        try
            [TimeStampsEV, ~, TTLs, ~, EventString Header] = ...
                Nlx2MatEV( fullfile(filedir, 'Events.nev'), [1 1 1 1 1], 1, ...
                1, 1 );
        catch
            continue
        end
        
        
        %% Find Event of Interest
        
        Q = floor(length(decimateCSC)/60000)-1;
        b = mod(length(decimateCSC), 60000);
        Trials = mat2cell( decimateCSC', 1, [ones(1, Q)*60000, 60000+b]);
        
        %         Trials{1}      = decimateCSC;
        tmp = table;
        tmp.trlfp = Trials';
        [Type{1:size(tmp,1)}] = deal('Task');
        tmp.type = Type';
        TaskTable = tmp;
        clear Type tmp
        
        %% Compute Power-Power Correlation
        Epochs = {'Task'};
        %         PowPowCor = cell(length(Epochs), 1);
        %         specgr = cell(length(Epochs), 1);
        PSD = cell(size(TaskTable, 1), 1);
        ov = 0.5;   % overlap
        winlen = 1024; % winlen
        nff = 400;
        freqBins = logspace(log10(1), log10(150), nff); %linspace(1, 150, nff); %logspace(log10(1), log10(200), nff);
        
        for ep = 1:size(TaskTable, 1)
            [s,freqs,time, pd] = spectrogram(TaskTable.trlfp{ep,1},winlen,floor(ov*winlen), freqBins, 1000);
            PSD{ep} = pd;
        end
        
        PSD = {cat(2, PSD{:})}; PSD = PSD{1};
        PowPow.session(Row) = {Session};
        PowPow.channel(Row) = channel;
        PowPow.Search(Row) = {PSD};
        PPCor = corrcoef(PSD');
        PowPowSearch = PPCor;
        PowPow.Corr(Row) = {PowPowSearch};
        
        
        %% Statistical Significance
        n_permutes = 5; pval = 0.05;
        
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
        imagesc(freqs, freqs, PPCor)
        axis xy
        load('vik.mat');
        colormap(jet); caxis([-0.3 0.3])
        cb = colorbar;
        xlabel('Frequency'); ylabel('Frequency')
        pbaspect([1 1 1])
        set(gca, 'box', 'off', 'FontSize', 15, 'tickdir', 'out', 'LineWidth', 1.5,...
            'xtick', [4 8 20 50 100 150], 'ytick', [4 8 20 50 100 150])
        
        subplot(122)
        zmap(abs(zmap) ~= 0) = 1;
        imagesc(freqs, freqs, zmap)
        axis xy
        caxis([0 1])
        cb = colorbar;
        xlabel('Frequency'); ylabel('Frequency')
        pbaspect([1 1 1])
        set(gca, 'box', 'off', 'FontSize', 15, 'tickdir', 'out', 'LineWidth', 1.5,...
            'xtick', [4 8 20 50 100 150], 'ytick', [4 8 20 50 100 150])
        
        title('LE Cross-freq power Correlation')
        
        %%
        figName = [Session '_' num2str(channellist(ch))]
        Figure_Output_Directory = 'D:\HiTGun\Data\PowPowCor Table Epoched'
        saveas(gcf, fullfile(Figure_Output_Directory, [figName '.png']))
        close(gcf)
        
        
        PowPow.Significance(Row) = {zmap};
        Row = Row + 1;
    end
end

LEPowPow = PowPow;
Directory = 'D:\HiTGun\Data\PowPowCor Table Epoched';
save(fullfile(Directory, 'PowPowPvalueTableLE_withStats_log2'), 'LEPowPow', 'freqs', '-v7.3')

