clear all, clc
addpath 'D:\Matlab Packages\fieldtrip-lite-20190224'


FilesDir = 'D:\HiTGun\Data\decimatedCSC_wholeSession';
TTLFilesDir = 'D:\HiTGun\Data\LE Change blindness Raw';

load('D:\HiTGun\Data\LETrialsTable');

info = struct2table(dir(FilesDir));
SessionTable = info.name;
SessionTable = SessionTable(3:end);

HilEnvCorr = table;
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
        
        %% Compute Power-Envelope Correlation
        
        freqs = 1:1:150;
        BandWidth = 2;
        hil_env = NaN(length(freqs), length(decimateCSC));
        for fr = 1:length(freqs)
            Af1       = freqs(fr);
            Af2       = Af1 + BandWidth;
            
            [b a]     = butter(3, [Af1 Af2]/(1000/2) );
            
            
            flt_sig   = filtfilt(b, a, decimateCSC);
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
        figName = [Session '_' num2str(channellist(ch))]
        Directory = 'D:\HiTGun\Data\HilEnvCorr';
        saveas(gcf, fullfile(Directory, [figName '.png']))
        close(gcf)
        
        Row = Row + 1;
    end
end

LEHilEnvCorr = HilEnvCorr;
save(fullfile(Directory, 'HilEnvCorr_LE'), 'LEHilEnvCorr', 'freqs', '-v7.3')

