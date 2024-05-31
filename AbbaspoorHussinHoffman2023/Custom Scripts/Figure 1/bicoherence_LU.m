clear all, clc
addpath 'D:\Matlab Packages\fieldtrip-lite-20190224'

load(fullfile('D:\HiTGun\Data','LU_In_Layer_Channel.mat'))
load('D:\HiTGun\Data\LUTrialsTable');

s = {'LU_2014-11-27_11-04-38', ...
    'LU_2014-12-02_11-09-10', ...
    'LU_2014-12-03_11-26-53', ...
    'LU_2014-12-05_11-33-27', ...
    'LU_2014-12-06_11-27-34', ...
    'LU_2014-12-10_11-27-48', ...
    'LU_2014-12-16_11-36-02', ...
    'LU_2014-12-18_11-28-38', ...
    'LU_2014-12-23_11-19-55', ...
    'LU_2014-12-26_11-13-14', ...
    'LU_2014-12-31_11-17-20', ...
    'LU_2015-01-05_11-05-18', ...
    'LU_2015-01-08_10-28-27', ...
    'LU_2015-01-09_11-16-46', ...
    'LU_2015-01-12_12-18-20', ...
    'LU_2015-01-13_11-05-52', ...
    'LU_2015-01-14_11-09-15', ...
    'LU_2015-01-15_11-12-18'};


restIdxDir = 'D:\HiTGun\restIdx'
FilesDir = 'D:\HiTGun\Data\LU Change blindness decimated';

SessionTable = s';

LUBicoherence = table;
Row = 1;

for Sess = 1:size(SessionTable, 1)
    
    if Sess == 8 || Sess == 17 || Sess == 28; continue; end
    %% Load Data
    Session = SessionTable{Sess};
    inchannel = LU_In_Layer_Channel.rip_cscNum(find(contains(string(LU_In_Layer_Channel.sid), Session)));
    ChannelList = [inchannel 29 30 31];
    ChannelList = unique(ChannelList);
    ChannelList = [29];
    
    for ch = 1:length(ChannelList)
        
        channel = {['csc', num2str(ChannelList(ch))]};
        
        info = struct2table(dir(fullfile(FilesDir, Session)));
        id = find(contains(info.name, channel));
        if isempty(id); continue; end
        channel = info.name(id);
        channel = {channel{1}(1:end-4)};
        
        %%% Load Channel LFP
        try
            load(fullfile(FilesDir, Session, channel{1}))
        catch
            continue
        end
        
        %% Reduce Line Noise
        downSampCSC = removeLineNoise_SpectrumEstimation(downSampCSC', 1000, 'NH = 1, LF = 60, M = 1024');
        downSampCSC = downSampCSC';
        
        %% Find Event of Interest
        
        Q = floor(length(downSampCSC)/60000)-1;
        b = mod(length(downSampCSC), 60000);
        Trials = mat2cell( downSampCSC', 1, [ones(1, Q)*60000, 60000+b]);
        
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
        figName = [Session '_' num2str(ChannelList(ch))]
        
        Figure_Output_Directory = 'D:\HiTGun\Data\Bicoherence Table Epoched'
        saveas(gcf, fullfile(Figure_Output_Directory, [figName '.png']))
        close(gcf)

        %%
        LUBicoherence.session(Row) = {Session};
        LUBicoherence.channel(Row) = channel;
        LUBicoherence.freq(Row) = {freq};
        LUBicoherence.Bic(Row) = {bico};
        LUBicoherence.significance(Row) = {zmap};
        
        Row = Row + 1;
        
        clear bico waxis BICOHERENCE zmap
        
        
    end
end

Directory = 'D:\HiTGun\Data\Bicoherence Table Epoched';
save(fullfile(Directory, 'LUBicoherence_withStats_2'), 'LUBicoherence', '-v7.3')


