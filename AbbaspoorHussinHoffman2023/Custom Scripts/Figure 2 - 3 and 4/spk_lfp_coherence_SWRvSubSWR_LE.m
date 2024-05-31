% clear all, clc
addpath 'D:\Matlab Packages\fieldtrip-lite-20190224'

FilesDir = 'D:\HiTGun\Data\decimatedCSC_wholeSession';
load('D:\HiTGun\Data\LETrialsTable');
load(fullfile('D:\HiTGun\', 'LESpikesTablePPC.mat'));
load('D:\HiTGun\Data\SWR_LEandLU'); clear table
% Tetrod = readtable('D:\HiTGun\Data\LE_TTNumbers.xlsx');
SessionTable = unique(Spikes.session);


% Row = 0;

for Sess = 2:size(SessionTable, 1)
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
    %     channelList = [2 6 13 15 26 35];
    %     channelList = [35];
    
    if ~isempty(channel)
        TT = Tetrod.Tetrod(find(contains(Tetrod.Channel, num2str(channel)))); TT = str2num(TT{1}(3:end));
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
        Row = findUnit(u-1);
        Spikes.Session(Row) = {Session};
        Spikes.Unit(Row) = u-1;
                
        channel = ['CSC', num2str(TTID(u-1, 2))];
        %%% Load Channel LFP
        try
            load(fullfile(FilesDir, Session, channel));
        catch
            continue
        end
        
        
        spkLFP = zeros(3, length(decimateCSC));
        spkLFP(1,:) = decimateCSC;
        Idx = interp1(TimeStamps(:,2), TimeStamps(:,1), [sessUnits{u-1, 1}.*1e6]);
        Idx = round(Idx); Idx(isnan(Idx)) = [];
        spkLFP(2, Idx) = 1;
        
        %% Find Event of Interest
        
        %%%%%%%%%%%% SWR
        id = find(contains(string(le.sid),Session));
        SWRmid = le.rip_mid_ts(id);
        SWRDur = le.rip_Dur(id);
        
        spkTrials = {};
        eventmid = NaN(1, length(SWRmid));
        for Tr = 1:length(SWRmid)
            eventmid(Tr) = interp1(TimeStamps(:,2), TimeStamps(:,1), [SWRmid(Tr)]); eventmid(Tr) = round(eventmid(Tr));
            half = ceil(SWRDur(Tr)/2);
            spkLFP(end, eventmid(Tr)-(half+600):eventmid(Tr)+(half+600)) = 1; % signal of interest in +/-250ms window around timestamp of interest
            % added 20ms before and 100 after actual SWR dur
        end
        
        
        %% Shuffle Test
        
        RippleWhere = find(spkLFP(3, :) == 1);
        tmpspkLFP = spkLFP; tmpspkLFP(2, RippleWhere) = 0;
        nsample = sum(tmpspkLFP(2, :));
        nsampleSWR = sum(spkLFP(2, :)) - nsample;

        if nsampleSWR > 20
            
            trials{1} = spkLFP([1 2], :);
            times{1} = [0:0.001:(length(spkLFP)/1000)-0.001];
            tmpspkLFP = [];
            tmpspkLFP.label   = {'LFP', 'Spike'}';
            tmpspkLFP.time    = times;
            tmpspkLFP.trial   = trials;
            tmpspkLFP.fsample = 1000;
            spk_LFP_Convol = spkspectrum(trials, times);
            
            
            NPermutation = 1000;
            ShuffleDiff = NaN([NPermutation, length(spk_LFP_Convol.freq)]);
            for Iteration = 1:NPermutation
                rng('shuffle')
                ShInd = randperm(sum(spkLFP(2, :)));
                NoSWR = ShInd(1:nsample);  SWR = ShInd(nsample+1:end);
                spk_LFP_Convol_SWR = spk_LFP_Convol;
                spk_LFP_Convol_SWR.fourierspctrm{1,1} = spk_LFP_Convol_SWR.fourierspctrm{1,1}(SWR, :, :);
                spk_LFP_Convol_SWR.time{1,1} = spk_LFP_Convol_SWR.time{1,1}(SWR, :);
                spk_LFP_Convol_SWR.trial{1,1} = spk_LFP_Convol_SWR.trial{1,1}(SWR, :);
                [spk_LFP_statSts] = stLFP_synchrony(spk_LFP_Convol_SWR, 'method', 'ppc0');
                stLFPsynSWR.SWR = spk_LFP_statSts.ppc0;
                
                
                spk_LFP_Convol_NoSWR = spk_LFP_Convol;
                spk_LFP_Convol_NoSWR.fourierspctrm{1,1} = spk_LFP_Convol_NoSWR.fourierspctrm{1,1}(NoSWR, :, :);
                spk_LFP_Convol_NoSWR.time{1,1} = spk_LFP_Convol_NoSWR.time{1,1}(NoSWR, :);
                spk_LFP_Convol_NoSWR.trial{1,1} = spk_LFP_Convol_NoSWR.trial{1,1}(NoSWR, :);
                [spk_LFP_statSts] = stLFP_synchrony(spk_LFP_Convol, 'method', 'ppc0');
                stLFPsynSWR.noSWR = spk_LFP_statSts.ppc0;
                ShuffleDiff(Iteration, :) = stLFPsynSWR.noSWR - stLFPsynSWR.SWR;
            end
            
            Spikes.ShuffleDiff(Row) = {ShuffleDiff};
            
            %         thresh_lo = prctile(max_val(:,1),    100*pval); % what is the
            %         thresh_hi = prctile(max_val(:,2),100-100*pval); % true p-value?
            
            %         SurrogateThresh = prctile(ShuffleDiff,100-(100*0.05), 1);
        else
            continue
        end
        
        
        %%
        CONT = 1;
        
        for ep = 1:2
            
            %%
            %%%%%%%% Prepare data
            tmpspkLFP = spkLFP;
            if ep == 1 % SWR
                fh = figure();
                fh.WindowState = 'maximized';
                
                RippleWhere = find(tmpspkLFP(3, :) == 0); % Find noSWR segments
                tmpspkLFP(2, RippleWhere) = 0;
                col = [87 136 108]/255;
            elseif ep == 2 % noSWR
                RippleWhere = find(tmpspkLFP(3, :) == 1); % Find SWR segments
                tmpspkLFP(2, RippleWhere) = 0;
                col = [14 15 25]/255;
            end
            nsample = sum(tmpspkLFP(2, :));
            if nsample < 20
                CONT = 0;
                continue
            end
            if CONT == 0
                CONT = 1;
                Row = Row + 1;
                close(gcf)
                continue
            end
            
            trials{1} = tmpspkLFP([1 2], :);
            times{1} = [0:0.001:(length(tmpspkLFP)/1000)-0.001];
            tmpspkLFP = [];
            tmpspkLFP.label   = {'LFP', 'Spike'}';
            tmpspkLFP.time    = times;
            tmpspkLFP.trial   = trials;
            tmpspkLFP.fsample = 1000;

            
            %%
            %%%%%%%%%%%% STA
            cfg              = [];
            cfg.timwin       = [-0.25 0.25]; % take 400 ms
            cfg.spikechannel = tmpspkLFP.label{2}; % first unit
            cfg.channel      = tmpspkLFP.label(1); % first four chans
            cfg.latency      = [0.3 10];
            staPost          = ft_spiketriggeredaverage(cfg, tmpspkLFP);
            if ep == 1; Spikes.STAswr(Row) = {staPost.avg}; end
            if ep == 2; Spikes.STAnoswr(Row) = {staPost.avg}; end
                
            % plot the sta
            if ep == 1; subplot(2, 7, [1 2]); end
            if ep == 2; subplot(2,7, [8 9]); end
            plot(staPost.time, staPost.avg(:,:)', 'Color', col, 'LineWidth', 2)
            set(gca, 'box', 'off', 'FontSize', 15, 'tickdir', 'out', 'LineWidth', 1.5)
            xlabel('time (s)')
            title('Spike Triggered Average')
            xlim([-0.25 0.25])
            set(gca, 'box', 'off', 'FontSize', 15, 'tickdir', 'out', 'LineWidth', 1.5,...
                'xtick', -0.25:0.05:0.25, 'xticklabel', -250:50:250)
            pbaspect([1 1 1])
            
            %% LFP Spectra
            spk_LFP_Convol = spkspectrum(trials, times);
            Spikes.freqRange(Row) = {spk_LFP_Convol.freq};
            
%             if ep == 1
%                 Spikes.spkspectrumSWR(Row) = spk_LFP_Convol;
%             elseif ep == 2
%                 Spikes.spkspectrumNoSWR(Row) = spk_LFP_Convol;
%             end
%             spk_LFP_Convol_SWR = spk_LFP_Convol;
%             spk_LFP_Convol_NoSWR = spk_LFP_Convol;
            
            %%
            %%%%%%%%%%%% Spike - LFP Synchronization (ppc0)
            [spk_LFP_statSts] = stLFP_synchrony(spk_LFP_Convol, 'method', 'ppc0');
            if ep == 1
                Spikes.PPCswr(Row) = {spk_LFP_statSts.ppc0};
            else
                Spikes.PPCnoswr(Row) = {spk_LFP_statSts.ppc0};
            end
            spk_LFP_statSts.ppc0(spk_LFP_statSts.ppc0 < 0) = 0;
            
            if ep == 1
                % plot the results
                ax11 = subplot(2, 7, 3)
                area(spk_LFP_statSts.freq(1:25), spk_LFP_statSts.ppc0(1:25)', 'Facecolor', col)
                set(gca, 'box', 'off', 'FontSize', 15, 'tickdir', 'out', 'LineWidth', 1.5)
                xlabel('frequency')
                pbaspect([1 1 1])
                xlim([spk_LFP_statSts.freq(1) spk_LFP_statSts.freq(25)])
                title(['PPC ', '(N: ', num2str(nsample), ')'])
                
                ax12 = subplot(2, 7, 4), hold on
                area(spk_LFP_statSts.freq(26:end), spk_LFP_statSts.ppc0(26:end)', 'facecolor', col)
                set(gca, 'box', 'off', 'FontSize', 15, 'tickdir', 'out', 'LineWidth', 1.5)
                pbaspect([1 1 1])
                xlim([spk_LFP_statSts.freq(26) spk_LFP_statSts.freq(end)])
                
            elseif ep == 2
                % plot the results
                ax21 = subplot(2, 7, 10)
                area(spk_LFP_statSts.freq(1:25), spk_LFP_statSts.ppc0(1:25)', 'Facecolor', col)
                set(gca, 'box', 'off', 'FontSize', 15, 'tickdir', 'out', 'LineWidth', 1.5)
                xlabel('frequency')
                pbaspect([1 1 1])
                xlim([spk_LFP_statSts.freq(1) spk_LFP_statSts.freq(25)])
                title(['PPC ', '(N: ', num2str(nsample), ')'])
                
                ax22 = subplot(2, 7, 11)
                area(spk_LFP_statSts.freq(26:end), spk_LFP_statSts.ppc0(26:end)', 'facecolor', col)
                set(gca, 'box', 'off', 'FontSize', 15, 'tickdir', 'out', 'LineWidth', 1.5)
                pbaspect([1 1 1])
                xlim([spk_LFP_statSts.freq(26) spk_LFP_statSts.freq(end)])
            end
            
            
            %%
            sub = {[5 6], [12 13]}
            if ep == 1
                [c, SWRMxPPCId] = max(spk_LFP_statSts.ppc0([3:10]));
            end
%             if ep == 1
%                 ppc_struct.swr(u) = c;
%             else
%                 ppc_struct.noswr(u) = c;
%             end
            
            fourierspctrm = squeeze(spk_LFP_Convol.fourierspctrm{1});
            angles = angle(fourierspctrm(:, SWRMxPPCId));
            angles = mod(angles,2*pi);%covert to 0-2pi rather than -pi:pi
            
            % compute Rayleigh's p-value
            [pval, z] = circ_rtest(angles);
            prefAngle = circ_mean(angles);
            itpc      = abs(mean(exp(1i*angles)));
            
            
            if ep == 1; subplot(2, 7, sub{1}(1)); end
            if ep == 2; subplot(2, 7, sub{2}(1)); end
            [acount, acenters] = hist(angles, 20);
            hist(angles,20);
            h = findobj(gca,'Type','patch');
            h.FaceColor = col;
            h.EdgeColor = 'w';
            hold on
            plot([0:0.01:2*pi],cos(0:0.01:2*pi)*max(acount)/2+max(acount)/2,'color',[0.5 0.5 0.5], 'LineWidth', 3)
            xlabel('Phase angle'), ylabel('Count')
            set(gca,'xlim',[0 2*pi], 'xtick', [0:0.5*pi:2*pi], 'xticklabel', {'0','pi/2','pi', '3/2pi', '2pi'})
            title({['Rayleigh p = ' num2str(pval)], ['Frequency: ' num2str(spk_LFP_statSts.freq(SWRMxPPCId))]})
            pbaspect([1 1 1])
            
            % and as polar distribution
            if ep == 1; subplot(2, 7, sub{1}(2)); end
            if ep == 2; subplot(2, 7, sub{2}(2)); end
            h = polarhistogram(angles, 40), hold on
            h.FaceColor = col;
            
            %%
            %%%%%%%%%%%% Save figure
            if ep == 2

                OrigDiff = Spikes.PPCnoswr{Row} - Spikes.PPCswr{Row};
                Spikes.PPCdiff(Row) = {OrigDiff};
%                 DiffSigvlue = NaN(size(OrigDiff)); DiffSiglog = zeros(size(OrigDiff));
%                 i = find(OrigDiff > SurrogateThresh);
%                 DiffSigvlue(i) = OrigDiff(i); DiffSiglog(i) = OrigDiff(i);
%                 Spikes.PPCdiffsig(Row) = {double(logical(DiffSiglog))};
                
                subplot(2, 7, 7)
                plot(spk_LFP_statSts.freq(1:25), OrigDiff(1:25)', 'Color', [0 0 0], 'LineWidth', 1)
%                 hold on
%                 plot(spk_LFP_statSts.freq(1:25), DiffSigvlue(1:25)', 'Color', [0 0 1], 'LineWidth', 1.5)
                set(gca, 'box', 'off', 'FontSize', 15, 'tickdir', 'out', 'LineWidth', 1.5)
                xlabel('frequency')
                pbaspect([1 1 1])
                xlim([spk_LFP_statSts.freq(1) spk_LFP_statSts.freq(25)])
                title({'SWR-noSWR', 'PPC Difference'})
                
                subplot(2, 7, 14)
                plot(spk_LFP_statSts.freq(26:end), OrigDiff(26:end)', 'Color', [0 0 0], 'LineWidth', 1)
%                 hold on
%                 plot(spk_LFP_statSts.freq(26:end), DiffSigvlue(26:end)', 'Color', [0 0 1], 'LineWidth', 1.5)
                set(gca, 'box', 'off', 'FontSize', 15, 'tickdir', 'out', 'LineWidth', 1.5)
                xlabel('frequency')
                pbaspect([1 1 1])
                xlim([spk_LFP_statSts.freq(26) spk_LFP_statSts.freq(end)])
                title({'SWR-noSWR', 'PPC Difference'})
                
                linkaxes([ax11, ax21], 'y');
                linkaxes([ax12, ax22], 'y');
                sgtitle ({['Session: ', Session], ['Channel:', channel], ['Unit:', num2str(u-1)]})
                figName = [Session,'_', channel,'_unit', num2str(u-1), '_SWR'];
                Figure_Output_Directory = 'D:\HiTGun\Figures\LE\Spike-Field Synchrony\SWR'
                saveas(gcf, fullfile(Figure_Output_Directory, [figName '.png']))
                %                     savefig(gcf, fullfile(Figure_Output_Directory, figName))
                close(gcf)

            end
        end
    end
end

save('LESpikesTablePPCSWR','Spikes', '-v7.3');

