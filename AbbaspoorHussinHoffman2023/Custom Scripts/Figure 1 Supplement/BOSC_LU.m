clear all, clc

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

BOSC = table;
Row = 1;
for Sess = 1:size(SessionTable, 1)
    if Sess == 8 || Sess == 12 || Sess == 15 || Sess == 17 || Sess == 28; continue; end
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
        %%%%%%%%%%%% REST
        restIdxinfo = struct2table(dir(fullfile(restIdxDir, Session(1:13))));
        restIdxID = find(contains(restIdxinfo.name, 'blankScr'));
        restIdxinfo = restIdxinfo(restIdxID, :);

        Trials = {};
        for Rs = 1:size(restIdxinfo, 1)
            load(fullfile(restIdxinfo.folder{Rs}, restIdxinfo.name{Rs}))
            reststart = interp1(TimeStamps(:,2), TimeStamps(:,1), [restIdx.nlx_rest_start_ts]);
            restend   = interp1(TimeStamps(:,2), TimeStamps(:,1), [restIdx.nlx_rest_end_ts]);
            reststart = round(reststart); restend = round(restend);  
            Trials{Rs}      = downSampCSC(reststart:restend); % signal of interest in +/-250ms window around timestamp of interest
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
        for Rs = 1:size(restIdxinfo, 1)
            load(fullfile(restIdxinfo.folder{Rs}, restIdxinfo.name{Rs}))
            reststart = interp1(TimeStamps(:,2), TimeStamps(:,1), [restIdx.nlx_rest_start_ts]);
            restend   = interp1(TimeStamps(:,2), TimeStamps(:,1), [restIdx.nlx_rest_end_ts]);
            reststart = round(reststart); restend = round(restend);   
            if Rs == 1
                rtime(1) = restend;
            elseif Rs ==2
                rtime(2) = reststart;
            end  
        end
        
        Trials{1} = downSampCSC(rtime(1):rtime(2)); % signal of interest in +/-250ms window around timestamp of interest


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


LUBOSC = BOSC;
Directory = 'D:\HiTGun\Data';
freqs = cfg.BOSC.F;
save(fullfile(Directory, 'BOSC_LU_RestSearch'), 'LUBOSC', 'freqs', '-v7.3')
