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

PowPow = table;
Row = 1;
for Sess = 1:size(SessionTable, 1)
    if Sess == 8 || Sess == 17; continue; end
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
        freqBins = logspace(log10(1), log10(150), nff); %linspace(1, 150, nff); %logspace(log10(1), log10(150), nff);
        
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
        figName = [Session '_' num2str(ChannelList(ch))]
        Figure_Output_Directory = 'D:\HiTGun\Data\PowPowCor Table Epoched'
        saveas(gcf, fullfile(Figure_Output_Directory, [figName '.png']))
        close(gcf)
        
        
        PowPow.Significance(Row) = {zmap};
        Row = Row + 1;
    end
end


LUPowPow = PowPow;
Directory = 'D:\HiTGun\Data\PowPowCor Table Epoched';
save(fullfile(Directory, 'PowPowPvalueTableLU_withStats_log'), 'LUPowPow', 'freqs', '-v7.3')
