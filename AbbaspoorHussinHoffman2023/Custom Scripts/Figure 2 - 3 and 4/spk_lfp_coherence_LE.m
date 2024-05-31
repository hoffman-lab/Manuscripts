clear all, clc
addpath 'D:\Matlab Packages\fieldtrip-lite-20190224'

FilesDir = 'D:\HiTGun\Data\decimatedCSC_wholeSession';
load('D:\HiTGun\Data\LETrialsTable');
load(fullfile('D:\HiTGun\Data', 'Spike_Table_LE.mat'));
load('D:\HiTGun\Data\SWR_LEandLU'); clear table
Tetrod = readtable('D:\HiTGun\Data\LE_TTNumbers.xlsx');

SessionTable = unique(Spikes.session);


Row = 1;
for Sess = 1:size(SessionTable, 1)
    %% Load Data
    Session = SessionTable{Sess};
    
    findUnit = find(contains(Spikes.session, Session));
    if isempty(findUnit); continue; end
    TTID = Spikes.TTID(findUnit);
    
    %%% Find in-layer-channel
    % This line loads the LE_In_Layer_Channel table that contains the channel
    % of recording that SWRs were found on. This is used as a marker to
    % selecting in-layer-channel.
    load(fullfile('D:\HiTGun\Data', 'LE_In_Layer_Channel.mat'));
    Index = find(contains(string(LE_In_Layer_Channel.sid),Session));
    channel = LE_In_Layer_Channel.rip_cscNum(Index);
    TT = Tetrod.Tetrod(find(contains(Tetrod.Channel, num2str(channel)))); TT = str2num(TT{1}(3:end));
    
    if ~isempty(channel)
        for CHANNEL = 1:length(TTID)
            if TTID(CHANNEL, 1) ~= TT
                TTID(CHANNEL, 2) = channel;
            elseif TTID(CHANNEL, 1) ~= 5
                TTID(CHANNEL, 2) = 7;
            else
                TTID(CHANNEL, 2) = 35;
            end
        end
    else
        for CHANNEL = 1:length(TTID)
            TTID(CHANNEL, 2) = 35;
        end
    end
    
    %%
    sessUnits = Spikes.Spike_ts(findUnit);
    numUnits = size(sessUnits, 1);
    
    %% Spike-LFP synchrony
    
    for u = 2:numUnits+1
        
        channel = ['CSC', num2str(TTID(u-1, 2))];
        %%% Load Channel LFP
        try
            load(fullfile(FilesDir, Session, channel));
        catch
            continue
        end
        
        %%
        spkLFP = zeros(2, length(decimateCSC));
        spkLFP(1,:) = decimateCSC;
        Idx = interp1(TimeStamps(:,2), TimeStamps(:,1), [sessUnits{u-1, 1}.*1e6]);
        Idx = round(Idx); Idx(isnan(Idx)) = [];
        spkLFP(2, Idx) = 1;
        
        %% Surrogate Shuffling Test
        %%% Part 1: Compute LFP spectra
        surrogate_spk_LFP_Convol = SPCshuffleTest(decimateCSC);
        
        %%
        Spikes.Session(Row) = {Session};
        Spikes.UnitNum(Row) = u-1;
        
        fh = figure();
        fh.WindowState = 'maximized';
        
        %%
        acg_metrics = calc_ACG_metrics(sessUnits{u-1}, 32000);
        Spikes.TMI(Row) = acg_metrics.thetaModulationIndex;
        
        %%
        %%%%%%%% Prepare data
        Batches = 1:10000:length(spkLFP);
        
        trials{1} = spkLFP;
        times{1} = [0:0.001:(length(spkLFP)/1000)-0.001];
        tmpspkLFP = [];
        tmpspkLFP.label   = {'LFP', 'Spike'}';
        tmpspkLFP.time    = times;
        tmpspkLFP.trial   = trials;
        tmpspkLFP.fsample = 1000;
        
        %% Surrogate Shuffling Test
        %%% Part 2
        count = sum(spkLFP(2, :));
        nPermutation = 1000;
        
        nfr = length(surrogate_spk_LFP_Convol.freq);
        SurrogatePPC = NaN(nPermutation, nfr);
        
        for permi = 1:nPermutation
            rng('shuffle')
            id = randperm(size(surrogate_spk_LFP_Convol.fourierspctrm{1, 1}, 1), count); id = sort(id);
            tmp_Convol = surrogate_spk_LFP_Convol;
            tmp_Convol.fourierspctrm = {tmp_Convol.fourierspctrm{1, 1}(id, :, :)};
            tmp_Convol.time = {tmp_Convol.time{1, 1}(id, 1)};
            tmp_Convol.trial = {tmp_Convol.trial{1, 1}(id, 1)};
            [spk_LFP_statSts] = stLFP_synchrony(tmp_Convol, 'method', 'ppc0');
            SurrogatePPC(permi, :) = spk_LFP_statSts.ppc0;
        end
        
        SurrogateThresh = prctile(SurrogatePPC,100-(100*0.05), 1);
        
        
        %%
        %%%%%%%%%%%% STA
        cfg              = [];
        cfg.timwin       = [-0.25 0.25]; % take 500 ms
        cfg.spikechannel = tmpspkLFP.label{2}; % first unit
        cfg.channel      = tmpspkLFP.label(1); % first four chans
        cfg.latency      = [0.3 10];
        staPost          = ft_spiketriggeredaverage(cfg, tmpspkLFP);
        Spikes.STA(Row)    = {staPost.avg};
        
        % plot the sta
        subplot(2, 3, 1)
        plot(staPost.time, staPost.avg(:,:)', 'Color', [0 0 0])
        set(gca, 'box', 'off', 'FontSize', 15, 'tickdir', 'out', 'LineWidth', 1.5)
        xlabel('time (s)')
        title('Spike Triggered Average')
        xlim([-0.25 0.25])
        set(gca, 'box', 'off', 'FontSize', 15, 'tickdir', 'out', 'LineWidth', 1.5,...
            'xtick', -0.25:0.05:0.25, 'xticklabel', -250:50:250)
        pbaspect([1 1 1])
        
        %% LFP Spectra
        spk_LFP_Convol = spkspectrum(trials, times);
        Spikes.spkspectrum(Row) = spk_LFP_Convol;
        
        %% Spike - LFP Synchronization
        %%%%%%%%%%%% (ppc0)
        [spk_LFP_statSts] = stLFP_synchrony(spk_LFP_Convol, 'method', 'ppc0');
        Spikes.PPC(Row) = {spk_LFP_statSts.ppc0};
        Significant = spk_LFP_statSts.ppc0 > SurrogateThresh; Significant = double(Significant);
        Significant(Significant==0) = NaN;
        spk_LFP_statSts.ppc0(spk_LFP_statSts.ppc0 < 0) = 0;
        Spikes.significant(Row) = {Significant};
        
        % plot the results
        ax1 = subplot(2, 3, 2)
        plot_shaded(spk_LFP_statSts.freq(1:25), spk_LFP_statSts.ppc0(1:25)', 'Color', [72 35 60]/255, 'Alpha', 1)
        hold on
        plot(spk_LFP_statSts.freq(1:25), Significant(1:25).*max(spk_LFP_statSts.ppc0(1:25))+0.0005, 'LineStyle', 'none', 'Marker', '.', 'Color', [0 0 0])
        set(gca, 'box', 'off', 'FontSize', 15, 'tickdir', 'out', 'LineWidth', 1.5)
        xlabel('frequency')
        pbaspect([1 1 1])
        xlim([spk_LFP_statSts.freq(1) spk_LFP_statSts.freq(25)])
        title('PPC')
        
        ax2 = subplot(2, 3, 3)
        plot_shaded(spk_LFP_statSts.freq(26:end), spk_LFP_statSts.ppc0(26:end)', 'Color', [72 35 60]/255, 'Alpha', 1)
        hold on
        plot(spk_LFP_statSts.freq(26:end), Significant(26:end).*max(spk_LFP_statSts.ppc0(26:end))+0.0005, 'LineStyle', 'none', 'Marker', '.', 'Color', [0 0 0])
        set(gca, 'box', 'off', 'FontSize', 15, 'tickdir', 'out', 'LineWidth', 1.5, 'YTickLabel',[])
        pbaspect([1 1 1])
        xlim([spk_LFP_statSts.freq(26) spk_LFP_statSts.freq(end)])
        linkaxes([ax1, ax2], 'y');
        
        %%
        sub = {[7 8], [9 10]}
        for fr = 1:2
            subplot(2, 3, 3+fr)
            freqq = {[1:9], [10:19]}; id = NaN(1,3);
            for f = 1:length(freqq)
                [c, id(f)] = max(spk_LFP_statSts.ppc0(freqq{f}));
                id(f) = freqq{f}(id(f));
            end
            fourierspctrm = squeeze(spk_LFP_Convol.fourierspctrm{1});
            angles = angle(fourierspctrm(:, id(fr)));
            angles = mod(angles,2*pi);%covert to 0-2pi rather than -pi:pi
            
            % compute Rayleigh's p-value
            [pval, z] = circ_rtest(angles);
            
            subplot(2, 6, sub{fr}(1))
            [acount, acenters] = hist(angles, 20);
            hist(angles,20);
            h = findobj(gca,'Type','patch');
            h.FaceColor = [72 35 60]/255;
            h.EdgeColor = 'w';
            hold on
            plot([0:0.01:2*pi],cos(0:0.01:2*pi)*max(acount)/2+max(acount)/2,'color',[0 0 0], 'LineWidth', 3)
            xlabel('Phase angle'), ylabel('Count')
            set(gca,'xlim',[0 2*pi], 'xtick', [0:0.5*pi:2*pi], 'xticklabel', {'0','pi/2','pi', '3/2pi', '2pi'})
            title({['Rayleigh p = ' num2str(pval)], ['Frequency: ' num2str(spk_LFP_statSts.freq(id(fr)))]})
            pbaspect([1 1 1])
            
            % and as polar distribution
            subplot(2, 6, sub{fr}(2))
            h = polarhistogram(angles, 40), hold on
            h.FaceColor = [72 35 60]/255;
        end
        
        subplot(2, 6, [11 12])
        bar(acg_metrics.thetaModulationIndex, 'FaceColor', [72 35 60]/255)
        title({['Theta Modulation'], 'From ACG'})
        set(gca, 'box', 'off', 'tickdir', 'out', 'LineWidth', 1.5)
        pbaspect([0.4 1 1])
        
        %%%%%%%%%%%% Save figure
        nsample = sum(spkLFP(2, :));
        sgtitle ({'Spike-Field Synchrony', ['Session: ', Session], ['Channel:', channel], ['Unit: ', num2str(u-1)],['Number of samples: ', num2str(nsample)]})
        figName = [Session,'_','unit', num2str(u-1), channel];
        Figure_Output_Directory = 'D:\HiTGun\Figures\LE\Spike-Field Synchrony'
        saveas(gcf, fullfile(Figure_Output_Directory, [figName '.png']))
        %             savefig(gcf, fullfile(Figure_Output_Directory, figName))
        close(gcf)
        Row = Row + 1;
        
    end
    
end


save('LESpikesTablePPC', 'Spikes',  '-v7.3')



