clear all, clc
addpath 'D:\Matlab Packages\fieldtrip-lite-20190224'


FilesDir = 'D:\HiTGun\Data\decimatedCSC_wholeSession';
TTLFilesDir = 'D:\HiTGun\Data\LE Change blindness Raw';

load('D:\HiTGun\Data\LETrialsTable');

info = struct2table(dir(FilesDir));
SessionTable = info.name;
SessionTable = SessionTable(3:end);

LEBicoherence = table;
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
    channellist = [inchannel]; %13 14 15
    channellist = unique(channellist);
    
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
        
        
        %% Compute Bicoherence
        Epochs = {'Task'};
        
        BICOHERENCE = cell(size(TaskTable, 1), 1);
        
        for ep = 1:size(TaskTable, 1)
            [bic,waxis] = bicoherence(TaskTable.trlfp{ep,1}, 'nsamp', 1024);
            BICOHERENCE{ep} = bic;
        end
        
        
        bico = cat(3,BICOHERENCE{:});
        bico = mean(bico, 3);
        
        freq = waxis.*1000;
        [c, f1] = min(abs(freq - 0));
        [c, f2] = min(abs(freq - 150));
        bico = bico(f1:f2, f1:f2);
        freq = freq(f1:f2);
        
        %% Statistical Significance
%         n_permutes = 5000; pval = 0.05; signal_length = 60000; Freq = [0 150];
%         
%         if Sess == 1 && ch == 1
%             [permmaps, mean_h0, std_h0, zval, cluster_thresh] = ...
%                 Bicoherence_Montecarlo_ClusterBased(n_permutes, Freq, pval, signal_length);
%             save('D:\HiTGun\Bicoherence_montecarlo_cluster.mat', 'permmaps', 'mean_h0', 'std_h0', 'zval', 'cluster_thresh')
%         end
        
        if Sess == 1 && ch == 1
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
        
%         zmap = zmap(f1:f2, f1:f2);

        %% Plot
        fh = figure();
        fh.WindowState = 'maximized';
        
        subplot(1, length(Epochs), 1)
        
        bico = triu(bico); bico = flip(bico);
        bico = tril(bico); bico = flip(bico);
        
        zmap = triu(zmap); zmap = flip(zmap);
        zmap = tril(zmap); zmap = flip(zmap);
        
        subplot(121)
        imagesc(freq, freq, bico)
        %         hold on
        %         contour(freq,freq,zmap,[1 1],'k','linewidth',5);
        axis xy; xlabel('Frequency'); ylabel('Frequency')
        cmap = brewermap([],'BuPu'); caxis([0 0.08]);
        load('devon.mat')
        colormap(flip(devon));
        cb = colorbar; set(cb,'Position',[0.92 0.37 0.0083 0.1974])
        set(gca, 'box', 'off', 'FontSize', 15, 'tickdir', 'out', 'LineWidth', 1.5,...
            'xtick', [8 16 32 50 100 150], 'ytick', [8 16 32 50 75])
        ylim([0 75])
        pbaspect([1 1 1])
        
        
        subplot(122)
        imagesc(freq, freq, zmap)
        axis xy; xlabel('Frequency'); ylabel('Frequency')
        cmap = brewermap([],'BuPu'); caxis([0 0.08]);
        load('devon.mat')
        colormap(flip(devon));
        cb = colorbar; set(cb,'Position',[0.92 0.37 0.0083 0.1974])
        set(gca, 'box', 'off', 'FontSize', 15, 'tickdir', 'out', 'LineWidth', 1.5,...
            'xtick', [8 16 32 50 100 150], 'ytick', [8 16 32 50 75])
        ylim([0 75])
        pbaspect([1 1 1])
        
        
        %%
        sgtitle ({'bicoherence', ['Session: ', Session, ' - Channel: ', channel{1}]})
        figName = [Session '_' num2str(channellist(ch))]
        
        Figure_Output_Directory = 'D:\HiTGun\Data\Bicoherence Table Epoched'
        saveas(gcf, fullfile(Figure_Output_Directory, [figName '.png']))
        close(gcf)
        
        %%
        LEBicoherence.session(Row) = {Session};
        LEBicoherence.channel(Row) = channel;
        LEBicoherence.freq(Row) = {freq};
        LEBicoherence.Bic(Row) = {bico};
        LEBicoherence.significance(Row) = {zmap};
        Row = Row + 1;
        
        clear bico waxis BICOHERENCE zmap
        
    end
end

Directory = 'D:\HiTGun\Data\Bicoherence Table Epoched';
save(fullfile(Directory, 'LEBicoherence_withStats_2'), 'LEBicoherence', '-v7.3')



