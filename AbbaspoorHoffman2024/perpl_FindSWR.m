%% Description
% This code detects sharp-wave ripple (SWR) events based on the methods described in Buzsaki
% Modified by Saman Abbaspoor - Feb 2020


% The SWR are detected based on threshold the preprocessed signal.
% Preprocessing:
% Signal is bandpass filtered in high frequency range [e.g. 100-250Hz]
% Then, signal is rectified,
% the rectified signal is lowpass filtered with a specific low cutoff frequency [e.g. 20Hz]
% the resultant signal is normalized

% The preprocessed signal is then thresholded with a low threshold criterion
% Then, the detected events that have a minimum interval less than a
% specific criterion [e.g. 30ms] are merged
% ripples that have a peak power less than the high threshold criterion are
% discarded
% Too short or too long detected events are also discarded

% If a noise channel (noise) is provided, the algorithm looks for
% overlapping detected events and if the event is also detected on the
% noise channel, that event is discarded.


% if MUA is provided a ovelapping threshold with MUA detected events is done.
% for more information, look at:
% Alyssa A. Carey, et al - Reward revaluation biases hippocampal replay content away from the preferred outcome
% also check the github https://github.com/vandermeerlab/vandermeerlab/blob/master/code-matlab/tasks/Replay_Analysis/getSWR.m


% Output
% The output is a matrix called ripples with 4 column which are:
% SWR start timestamp, SWR peak timestamp, SWR end timestamp
% and peak amplitude

% timestamps must be in seconds

function ripples = perpl_FindSWR(lfp, noise, mua, timestamps, varargin)
%% Input Parameters and Data
% Input Parameters
p = inputParser;
addParameter(p,'frequency',1000,@isnumeric)           % Sampling Frequency [e.g. 1000]
addParameter(p,'thresholds',[1 3],@isvector)         % Envelope Threshold [e.g. 2 5 ]
addParameter(p,'durations',[20 40 500],@isnumeric)    % inter-ripple interval, minimum ripple duration and maximum ripple duration (in ms)
addParameter(p,'overlap_thr',10,@isnumeric)           
% addParameter(p,'downsample_factor',30,@isnumeric)      %  downsample_factor = sr_recording/sr_lfp (e.g 30000/1000)    
addParameter(p,'freq_range',[100 180],@isnumeric)
parse(p,varargin{:})

% Assign Parameters
frequency                = p.Results.frequency;
lowThresholdFactor       = p.Results.thresholds(1);
highThresholdFactor      = p.Results.thresholds(2);
minInterRippleInterval   = p.Results.durations(1);
minRippleDuration        = p.Results.durations(2);
maxRippleDuration        = p.Results.durations(3);
overlap_thr              = p.Results.overlap_thr;
% downsample_factor        = p.Results.downsample_factor;
freq_range               = p.Results.freq_range;

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  SharpWave Ripple Detection  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Prepare Signal

% 1- Bandpass Filter in High Frequency Range [e.g. 100-250Hz]
[b a] = butter(3, freq_range/ (frequency/2) );
signal = filtfilt(b, a, lfp);

% 2- Rectify the filtered signal,
rectifiedSignal = signal.^2; %abs(signal);

% 3- Bandpass Filter signal in Low Frequency Range [e.g. 20Hz]
% 4- Normalize the Signal

[b a] = butter(3, [1 20]/(frequency/2) );
normalizedSquaredSignal = zscore( filtfilt(b, a, rectifiedSignal) );


%% Detect ripple periods by thresholding normalized squared signal

thresholded = normalizedSquaredSignal > lowThresholdFactor;
start = find( diff(thresholded) > 0 );
stop = find( diff(thresholded)<0 );

% Exclude last ripple if it is incomplete
if length(stop) == length(start)-1,
    start = start(1:end-1);
end

% Exclude first ripple if it is incomplete
if length(stop)-1 == length(start),
    stop = stop(2:end);
end

% Correct special case when both first and last ripples are incomplete
if start(1) > stop(1),
    stop(1) = [];
    start(end) = [];
end


firstPass = [start,stop];
if isempty(firstPass),
    disp('Detection by thresholding failed');
    return
else
    disp(['After detection by thresholding: ' num2str(length(firstPass)) ' events.']);
end


%% Merge ripples if inter-ripple period is too short

minInterRippleSamples = minInterRippleInterval/1000*frequency;

secondPass = [];
ripple = firstPass(1,:);

for i = 2:size(firstPass,1)
    
    if firstPass(i,1) - ripple(2) < minInterRippleSamples
        
        % Merge
        ripple = [ripple(1) firstPass(i,2)];
        
    else
        
        secondPass = [secondPass ; ripple];
        ripple = firstPass(i,:);
        
    end
    
end

secondPass = [secondPass ; ripple];
if isempty(secondPass),
    disp('Ripple merge failed');
    return
else
    disp(['After ripple merge: ' num2str(length(secondPass)) ' events.']);
end

%% Discard ripples with a peak power < highThresholdFactor

thirdPass = [];
peakNormalizedPower = [];
peakPosition = [];
for i = 1:size(secondPass,1)
    
    [maxValue,maxIndex] = max(normalizedSquaredSignal([secondPass(i,1):secondPass(i,2)]));
    if maxValue > highThresholdFactor,
        thirdPass = [thirdPass ; secondPass(i,:)];
        peakNormalizedPower = [peakNormalizedPower ; maxValue];
%         maxIndex = maxIndex + secondPass(i,1) - 1;
%         peakPosition = [peakPosition, maxIndex];
    end
    
end


if isempty(thirdPass)
    disp('Peak thresholding failed.');
    ripples = 0;
    return
else
    disp(['After peak thresholding: ' num2str(length(thirdPass)) ' events.']);
end


%%

if ~isempty(thirdPass)
    
%     Detect negative peak position for each ripple
    peakPosition = zeros(size(thirdPass,1),1);
    for i=1:size(thirdPass,1)
        [minValue,minIndex] = min(signal(thirdPass(i,1):thirdPass(i,2)));
        peakPosition(i) = minIndex + thirdPass(i,1) - 1;
    end
    
    
    ripples = [timestamps(thirdPass(:,1)) timestamps(peakPosition) ...
        timestamps(thirdPass(:,2)) peakNormalizedPower];
    
    
%     duration = (ripples(:,3)-ripples(:,1))/downsample_factor;
    duration = (ripples(:,3)-ripples(:,1));
    

    % Discard ripples that are way too long
    ripples(duration>maxRippleDuration/1000,:) = NaN;
    
    % Discard ripples that are too short
    ripples(duration<minRippleDuration/1000,:) = NaN;
    
    ripples = ripples((all((~isnan(ripples)),2)),:);
    disp(['After duration test: ' num2str(size(ripples,1)) ' events.']);
    
    
end


%% If a noise channel was provided, find ripple-like events and exclude them

bad = [];
if ~isempty(noise)
    
    % 1- Bandpass Filter in High Frequency Range [e.g. 100-250Hz]
    [b a] = butter(3, freq_range/ (frequency/2) );
    noise = filtfilt(b, a, double(noise));
    rectifiedNoise = abs(noise);    
    [b a] = butter(3, [1 20]/(frequency/2) );
    normalizedSquaredNoise = zscore( filtfilt(b, a, rectifiedNoise) );
    
    
    excluded = logical(zeros(size(ripples,1),1));
    
    % Exclude ripples when concomittent noise crosses high detection threshold
    for i = 1:size(ripples,1)
        
        tmp_ts = find( timestamps(:,1) >= ripples(i,1) & timestamps(:,1) <= ripples(i,3) );
        indices = [tmp_ts(1), tmp_ts(end)];

        if any(normalizedSquaredNoise(indices(1):indices(2))>highThresholdFactor) %highThresholdFactor %highThresholdFactor
            excluded(i) = 1;
        end
        
    end
    
    bad = ripples(excluded,:);
    ripples = ripples(~excluded,:);
    disp(['After ripple-band noise removal: ' num2str(size(ripples,1)) ' events.']);
end


%% MUA overlap Thresholding
keep_idx = [];

if ~isempty(mua)
    
    EVENT = mua_Detection(mua, timestamps);
    numEvent = length(EVENT);
    
    for iIV = 1:length(ripples)
        curr_tstart = ripples(iIV, 1) - overlap_thr;
        curr_tend   = ripples(iIV, 3) + overlap_thr;
        
        check1 = any(find(curr_tstart < EVENT(:, 1) & EVENT(:, 1) < curr_tend));
        check2 = any(find(curr_tstart < EVENT(:, 2) & EVENT(:, 2) < curr_tend));
        
        if check1 || check2
            keep_idx = cat(1,keep_idx,iIV);
        end
        
    end
    
    ripples = ripples(keep_idx, :);
    
    disp(['After MUA thresholding: ' num2str(size(ripples,1)) ' events.']);
end

end

