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

HilEnvCorr = table;
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
        
        %% Compute Power-Envelope Correlation
        
        freqs = 1:1:150;
        BandWidth = 2;
        hil_env = NaN(length(freqs), length(downSampCSC));
        for fr = 1:length(freqs)
            Af1       = freqs(fr);
            Af2       = Af1 + BandWidth;
            
            [b a]     = butter(3, [Af1 Af2]/(1000/2) );
            
            flt_sig   = filtfilt(b, a, downSampCSC);
            hil_sig   = abs(hilbert(flt_sig));
            hil_sig   = zscore(hil_sig);
            hil_env(fr,:) = hil_sig;
        end
        
        HilEnvCorr.session(Row) = {Session};
        HilEnvCorr.channel(Row) = channel;
        env_corr = corrcoef(hil_env');
        EnvCorrSearch = env_corr;
        HilEnvCorr.Corr(Row) = {EnvCorrSearch};        
        
        %%
        fh = figure();
        fh.WindowState = 'maximized';
        
        subplot(121)
        imagesc(freqs, freqs, env_corr)
        axis xy
        RdBu=cbrewer('div', 'RdBu', 101);
        colormap(flip(RdBu)); caxis([-0.2 0.2])
        cb = colorbar;
        xlabel('Frequency'); ylabel('Frequency')
        pbaspect([1 1 1])
        set(gca, 'box', 'off', 'FontSize', 15, 'tickdir', 'out', 'LineWidth', 1.5,...
            'xtick', [4 8 20 50 100 150], 'ytick', [4 8 20 50 100 150])
        
        title('LE Cross-freq power Correlation')
        
        %%
        figName = [Session '_' num2str(ChannelList(ch))]
        Directory = 'D:\HiTGun\Data\HilEnvCorr';
        saveas(gcf, fullfile(Directory, [figName '.png']))
        close(gcf)
        Row = Row + 1;
    end
end


LUHilEnvCorr = HilEnvCorr;
save(fullfile(Directory, 'HilEnvCorr_LU'), 'LUHilEnvCorr', 'freqs', '-v7.3')

