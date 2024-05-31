clear all, clc
FilesDir = 'D:\HiTGun\Data\decimatedCSC_wholeSession';
TTLFilesDir = 'D:\HiTGun\Data\LE Change blindness Raw';

load('D:\HiTGun\Data\LETrialsTable');

info = struct2table(dir(FilesDir));
SessionTable = info.name;
SessionTable = SessionTable(3:end);

BOSC = table;
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
    channellist = inchannel;
    
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
        %%%%%%%%%%%% REST
        TaskStarts = find(TTLs == 23); TaskStarts = TaskStarts(end);
        TaskEnds = find(TTLs == 24); TaskEnds = TaskEnds(end);
        
        Trials_ev = NaN(2, 1);
        Trials_ev(1) = TimeStampsEV(TaskStarts);
        Trials_ev(2) = TimeStampsEV(TaskEnds);
        
        Trials = {};
        for Tr = 1:length(Trials_ev)
            rtime = interp1(TimeStamps(:,2), TimeStamps(:,1), [Trials_ev(Tr,1)]); rtime = round(rtime);
            if Tr == 1
                Trials{Tr}      = decimateCSC(1:rtime);
            elseif Tr ==2
                Trials{Tr}      = decimateCSC(rtime:end);
            end
        end
        
        tmp = table;
        tmp.trlfp = Trials';
        Type{1} ='PreRest';
        Type{2} = 'PostRest';
        tmp.type = Type';
        TaskTable = tmp;
        clear Type tmp
        
        Trials = {};
        rtime = NaN(1, 2);
        for Tr = 1:length(Trials_ev)
            rtime(Tr) = interp1(TimeStamps(:,2), TimeStamps(:,1), [Trials_ev(Tr,1)]); rtime = round(rtime);
        end
        Trials{1}      = decimateCSC(rtime(1):rtime(2));
        
        tmp = table;
        tmp.trlfp = Trials';
        [Type{1:size(tmp,1)}] = deal('Search');
        tmp.type = Type';
        TaskTable = [TaskTable; tmp];
        clear Type tmp
        
        %% Compute Power-Power Correlation
        Epochs = {'Search', 'PreRest', 'PostRest'};
        
        for ep = 1:length(Epochs)
            
            Index = find(contains(TaskTable.type,Epochs{ep}));
            signal = TaskTable.trlfp{Index,1}';
            
            % general setup
            cfg.BOSC.F             = 1:100;
            cfg.BOSC.wavenumber	   = 6;
            cfg.BOSC.fsample       = 1000;
            
            cfg.BOSC.threshold.duration	   = 3;         % vector of duration thresholds at each frequency (previously: ncyc)
            cfg.BOSC.threshold.percentile  = .95;       % percentile of background fit for power threshold
            ThetaRange = 6:9;
            GammaRange = 20:50;
            
            [TF,T,F]=BOSC_tf(signal, cfg.BOSC.F, cfg.BOSC.fsample, cfg.BOSC.wavenumber);
            [pv,meanpower]=BOSC_bgfit(F,TF);
            [powthresh,durthresh]=BOSC_thresholds(1000, cfg.BOSC.threshold.percentile, cfg.BOSC.threshold.duration, cfg.BOSC.F, meanpower);
            
            detected = zeros(size(TF));
            for f = 1:length(cfg.BOSC.F)
                detected(f,:) = BOSC_detect(TF(f,:),powthresh(f),durthresh(f),1000);
            end; clear f
            
            Occupancy = NaN(1, size(detected, 1));
            
            for fr = 1:size(detected, 1)
                Occupancy(fr) = numel(find(detected(fr, :) == 1))/size(detected, 2);
            end
            
            BOSC.session(Row) = {Session};
            BOSC.channel(Row) = channel;
            
            
            
            thresholded = sum(detected(ThetaRange,:))>=1;
            start = find( diff(thresholded) > 0 );
            stop = find( diff(thresholded)<0 );
            if length(stop) == length(start)-1; start = start(1:end-1); end
            if length(stop)-1 == length(start); stop = stop(2:end); end
            if start(1) > stop(1); stop(1) = []; start(end) = []; end
            ThetaBurstTimes = stop-start;
            
            
            thresholded = sum(detected(GammaRange,:))>=1;
            start = find( diff(thresholded) > 0 );
            stop = find( diff(thresholded)<0 );
            if length(stop) == length(start)-1; start = start(1:end-1); end
            if length(stop)-1 == length(start); stop = stop(2:end); end
            if start(1) > stop(1); stop(1) = []; start(end) = []; end
            GammaBurstTimes = stop-start;
            
            
            if ep == 1;
                BOSC.Search(Row) = {Occupancy};
                BOSC.SearchThetaDur(Row) = {ThetaBurstTimes};
                BOSC.SearchGammaDur(Row) = {GammaBurstTimes};
                
            elseif ep == 2;
                BOSC.PreRest(Row) = {Occupancy};
                BOSC.PreRestThetaDur(Row) = {ThetaBurstTimes};
                BOSC.PreRestGammaDur(Row) = {GammaBurstTimes};
                
            elseif ep == 3;
                BOSC.PostRest(Row) = {Occupancy};
                BOSC.PostRestThetaDur(Row) = {ThetaBurstTimes};
                BOSC.PostRestGammaDur(Row) = {GammaBurstTimes};
                
                Row = Row + 1;
            end
            
        end

    end
end

LEBOSC = BOSC;
Directory = 'D:\HiTGun\Data';
freqs = cfg.BOSC.F;
save(fullfile(Directory, 'BOSC_LE_RestSearch'), 'LEBOSC', 'freqs', '-v7.3')

