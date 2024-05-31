clear all, clc
%% Load Data
filedir = 'D:\HiTGun\Wi TH Data\Processed';
Sessions = struct2table(dir(filedir));
Sessions = Sessions(3:end, :);

BOSC = table;
for sess = 1:size(Sessions, 1)-1
    
    Session = Sessions.name{sess};
    BOSC.session(sess) = {Session};
    
    %% Parameters
    
    cfg.BOSC.F             = 1:150;
    cfg.BOSC.wavenumber	   = 6;
    cfg.BOSC.fsample       = 1000;
    
    cfg.BOSC.threshold.duration	   = 3;         % vector of duration thresholds at each frequency (previously: ncyc)
    cfg.BOSC.threshold.percentile  = .95;       % percentile of background fit for power threshold
    ThetaRange = 6:9;
    GammaRange = 20:50;
    
    %% Detect Oscillatory Bouts - REST
    
    
    Occupancy = NaN(1, length(cfg.BOSC.F));
    
    load(fullfile(filedir, Sessions.name{sess}, 'Sleep'))
    
    [TF,T,F]=BOSC_tf(zscore(lfp), cfg.BOSC.F, cfg.BOSC.fsample, cfg.BOSC.wavenumber);
    [pv,meanpower]=BOSC_bgfit(F,TF);
    [powthresh,durthresh]=BOSC_thresholds(1000, cfg.BOSC.threshold.percentile, cfg.BOSC.threshold.duration, cfg.BOSC.F, meanpower);
    
    detected = zeros(size(TF));
    for f = 1:length(cfg.BOSC.F)
        detected(f,:) = BOSC_detect(TF(f,:),powthresh(f),durthresh(f),1000);
    end; clear f
    
    for fr = 1:size(detected, 1)
        Occupancy(fr) = numel(find(detected(fr, :) == 1))/size(detected, 2);
    end
    
    BOSC.Rest(sess) = {Occupancy};
    
    
    
    thresholded = sum(detected(ThetaRange,:))>=1;
    start = find( diff(thresholded) > 0 );
    stop = find( diff(thresholded)<0 );
    if length(stop) == length(start)-1; start = start(1:end-1); end
    if length(stop)-1 == length(start); stop = stop(2:end); end
    if start(1) > stop(1); stop(1) = []; start(end) = []; end
    BurstTimes = stop-start;
    BOSC.RestThetaDur(sess) = {BurstTimes};
    
    thresholded = sum(detected(GammaRange,:))>=1;
    start = find( diff(thresholded) > 0 );
    stop = find( diff(thresholded)<0 );
    if length(stop) == length(start)-1; start = start(1:end-1); end
    if length(stop)-1 == length(start); stop = stop(2:end); end
    if start(1) > stop(1); stop(1) = []; start(end) = []; end
    BurstTimes = stop-start;
    BOSC.RestGammaDur(sess) = {BurstTimes};
    
    
    
    %% Detect Oscillatory Bouts - SEARCH
    
    load(fullfile(filedir, Sessions.name{sess}, 'Treehouse'))
    to_remove = cellfun(@(x) isempty(x), TrialLFP);
    TrialLFP(to_remove) = []; TrialTime(to_remove) = [];
    
    
    Occupancy = NaN(length(TrialLFP), length(cfg.BOSC.F));
    startTheta=[]; stopTheta=[]; startGamma=[]; stopGamma=[];
    
    for TRIAL = 1:length(TrialLFP)
        
        [TF,T,F]=BOSC_tf(zscore(TrialLFP{1, TRIAL}), cfg.BOSC.F, cfg.BOSC.fsample, cfg.BOSC.wavenumber);
        [pv,meanpower]=BOSC_bgfit(F,TF);
        [powthresh,durthresh]=BOSC_thresholds(1000, cfg.BOSC.threshold.percentile, cfg.BOSC.threshold.duration, cfg.BOSC.F, meanpower);
        
        detected = zeros(size(TF));
        for f = 1:length(cfg.BOSC.F)
            detected(f,:) = BOSC_detect(TF(f,:),powthresh(f),durthresh(f),1000);
        end; clear f
        
        for fr = 1:size(detected, 1)
            Occupancy(TRIAL, fr) = numel(find(detected(fr, :) == 1))/size(detected, 2);
        end
        
        thresholded = sum(detected(ThetaRange,:))>=1;
        tmpStart = find( diff(thresholded) > 0 );
        tmpStop = find( diff(thresholded)<0 );
        if ~isempty(tmpStart) & ~isempty(tmpStop)
            if length(tmpStop) == length(tmpStart)-1; tmpStart = tmpStart(1:end-1); end
            if length(tmpStop)-1 == length(tmpStart); tmpStop = tmpStop(2:end); end
            if tmpStart(1) > tmpStop(1); tmpStop(1) = []; tmpStart(end) = []; end
            startTheta = [startTheta, tmpStart];
            stopTheta  = [stopTheta, tmpStop];
        end
        
        
        thresholded = sum(detected(GammaRange,:))>=1;
        tmpStart = find( diff(thresholded) > 0 );
        tmpStop = find( diff(thresholded)<0 );
        if ~isempty(tmpStart) & ~isempty(tmpStop)
            if length(tmpStop) == length(tmpStart)-1; tmpStart = tmpStart(1:end-1); end
            if length(tmpStop)-1 == length(tmpStart); tmpStop = tmpStop(2:end); end
            if tmpStart(1) > tmpStop(1); tmpStop(1) = []; tmpStart(end) = []; end
        end
        startGamma = [startGamma, tmpStart];
        stopGamma  = [stopGamma, tmpStop];
    end
    
    
    Occupancy = mean(Occupancy, 1);
    BOSC.Search(sess) = {Occupancy};
    
    BurstTimes = stopTheta-startTheta;
    BOSC.SearchThetaDur(sess) = {BurstTimes};
    BurstTimes = stopGamma-startGamma;
    BOSC.SearchGammaDur(sess) = {BurstTimes};
    
    %%
%     Occupancy_search = BOSC.Search{sess};
%     Col = [255 8 32]/255;
%     hold on
%     area(cfg.BOSC.F, Occupancy_search, 'EdgeColor', 'None',...
%         'FaceColor', [Col], 'FaceAlpha', 0.5)
%     
%     Occupancy_rest = BOSC.Rest{sess};
%     
%     Col = [56 138 191]/255;
%     area(cfg.BOSC.F, Occupancy_rest, 'EdgeColor', 'None',...
%         'FaceColor', [Col], 'FaceAlpha', 0.5)
%     
%     set(gca, 'LineWidth', 3, 'Box', 'off', 'Tickdir', 'out', 'FontSize', 15)
%     ylabel('Occupancy Rate (%)')
%     xlabel('Frequency (Hz)')
%     xlim([1 70])
%     
%     legend('Search', 'Rest')
%     
%     %%
%     figName = Session;
    Directory = 'D:\HiTGun\Wi TH Data\Analyses\BOSC'
%     saveas(gcf, fullfile(Directory, [figName '.png']))
%     close(gcf)
    
end

WIBOSC = BOSC;
freqs = cfg.BOSC.F;
save(fullfile(Directory, 'BOSC_Wi_RestSearch_forDurDist'), 'WIBOSC', 'freqs', '-v7.3')



%%



% general setup



if ep == 1;
    BOSC.Search(Row) = {Occupancy};
elseif ep == 2;
    BOSC.PreRest(Row) = {Occupancy};
elseif ep == 3;
    BOSC.PostRest(Row) = {Occupancy};
    Row = Row + 1;
end

end


Occupancy_search = BOSC.Search{Row-1};
Col = [255 8 32]/255;
hold on
area(cfg.BOSC.F, Occupancy_search, 'EdgeColor', 'None',...
    'FaceColor', [Col], 'FaceAlpha', 0.5)

Occupancy_rest = (BOSC.PreRest{Row-1}+BOSC.PostRest{Row-1})/2;

Col = [56 138 191]/255;
area(cfg.BOSC.F, Occupancy_rest, 'EdgeColor', 'None',...
    'FaceColor', [Col], 'FaceAlpha', 0.5)

set(gca, 'LineWidth', 3, 'Box', 'off', 'Tickdir', 'out', 'FontSize', 15)
ylabel('Occupancy Rate (%)')
xlabel('Frequency (Hz)')
xlim([1 70])

legend('Search', 'Rest')
