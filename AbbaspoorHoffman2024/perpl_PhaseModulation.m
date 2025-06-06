%% Phase Modulation
%%% Modified from %[PhaseLockingData] = bz_PhaseModulation(varargin)
%%% Saman Abbaspoor 03/12/2024 - Hoffman Lab - Vanderbilt

function [PhaseLockingData] = perpl_PhaseModulation(spikes, lfp, passband, varargin)
% USAGE
%[PhaseLockingData] = bz_PhaseModulation(varargin)
%
% INPUTS
% spikes        -spike time cellinfo struct
%
% lfp           -lfp struct with a single channel from bz_GetLFP()
%
% passband      -frequency range for phase modulation [lowHz highHz] form
%
% intervals     -(optional) may specify timespans over which to calculate
%               phase modulation.  Formats accepted: tstoolbox intervalSet
%               or a 2column matrix of [starts stops] in seconds
%
% samplingRate  -specifies lfp sampling frequency default=1250
%
% method        -method selection for how to generate phase,
%               possibilties are: 'hilbert' (default) or 'wavelet'
%
% powerThresh   -integer power threshold to use as cut off,
%               measured in standard deviations (default = 2)
%
% plotting      -logical if you want to plot, false if not, default=true
%
% saveMat       -logical to save cellinfo .mat file with results, default=false
%
%
% OUTPUTS
%
% phasedistros  - Spike distribution perecentages for each cell in each bin
%               specified by phasebins
%
% phasebins     - 180 bins spanning from 0 to 2pi
%
% phasestats    - ncellsx1 structure array with following (via
%                 CircularDistribution.m from FMAToolbox)
%                    phasestats.m        mean angle
%                    phasestats.mode     distribution mode
%                    phasestats.k        concentration
%                    phasestats.p        p-value for Rayleigh test
%                    phasestats.r        mean resultant length
%
%
% Calculates distribution of spikes over various phases from a specified
% cycle of an lfp vector.   Phase 0 means peak of lfp wave.
%
% Brendon Watson 2015
% edited by david tingley, 2017

%% defaults
p = inputParser;
addParameter(p,'intervals',[0 inf],@isvector)
addParameter(p,'samplingRate',1000,@isnumeric)
addParameter(p,'method','hilbert',@isstr)
addParameter(p,'Thresholded',[],@isvector)
addParameter(p,'remove_sample',[],@isvector)
addParameter(p,'plotting',false,@islogical)
addParameter(p,'numBins',180,@isnumeric)
addParameter(p,'powerThresh',2,@isnumeric)
addParameter(p,'saveMat',false,@islogical)

parse(p,varargin{:})


intervals = p.Results.intervals; % interval(s) over which to calculate
samplingRate = p.Results.samplingRate; % sampling rate of continuous signal (LFP)
method = p.Results.method;
Thresholded = p.Results.Thresholded;
remove_sample = p.Results.remove_sample;
plotting = p.Results.plotting;
numBins = p.Results.numBins;
powerThresh = p.Results.powerThresh;
saveMat = p.Results.saveMat;

%% Get phase for every time point in LFP
switch lower(method)
    case ('hilbert')
        [b a] = butter(3, passband/(samplingRate/2)); % order 3
        
        %         [b a] = butter(3,[passband(1)/(samplingRate/2) passband(2)/(samplingRate/2)],'bandpass'); % order 3
        %         [b a] = cheby2(4,20,passband/(samplingRate/2));
        filt = FiltFiltM(b,a,double(lfp.data(:,1)));
        %         power = fastrms(filt,ceil(samplingRate./passband(1)));  % approximate power is frequency band
        hilb = hilbert(filt);
        power = abs(hilb);
        lfpphase = mod(angle(hilb),2*pi);
        clear fil
    case ('wavelet')% Use Wavelet transform to calulate the signal phases
        [wave,f,t,coh,wphases,raw,coi,scale,priod,scalef]=getWavelet(double(lfp.data(:,1)),samplingRate,passband(1),passband(2),8,0);
        [~,mIdx]=max(wave);%get index max power for each timepiont
        pIdx=mIdx'+[0;size(f,2).*cumsum(ones(size(t,1)-1,1))];%converting to indices that will pick off single maxamp index from each of the freq-based phases at eacht timepoint
        lfpphase=wphases(pIdx);%get phase of max amplitude wave at each timepoint
        lfpphase = mod(lfpphase,2*pi);%covert to 0-2pi rather than -pi:pi
        power = rms(abs(wave))';
end

%% update intervals to remove sub-threshold power periods
if isempty(Thresholded)
    Thresholded = zscore(power)>powerThresh;
    Thresholded(remove_sample) = 0;
    Thresholded = find(Thresholded);
else
    Thresholded(remove_sample) = 0;
    Thresholded=find(Thresholded);
end

%% Get phases for each spike for each cell
h = [];
% cum_spkphases = [];
phasebins=[];
spkphases = cell(1,length(spikes.times));
for a = 1:length(spikes.times)
    s = spikes.times{a}(ismember(spikes.times{a}, Thresholded));
    
    if isempty(s)
        phasedistros(:,a) = zeros(numBins,1);
        phasestats.m(a) = nan;
        phasestats.r(a) = nan;
        phasestats.k(a) = nan;
        phasestats.p(a) = nan;
        phasestats.mode(a) = nan;
        phasestats.samplenum(a) = {nan};
        spkphases{a} = nan;
    else
        spkphases{a} = lfpphase(s);
        
        %% Gather binned counts and stats (incl Rayleigh Test)
        [phasedistros(:,a),phasebins,ps]=CircularDistribution(spkphases{a},'nBins',numBins);
        phasestats.m(a) = mod(ps.m,2*pi);
        phasestats.r(a) = ps.r;
        phasestats.k(a) = ps.k;
        phasestats.p(a) = ps.p;
        phasestats.mode(a) = ps.mode;
        phasestats.samplenum(a) = {s};
        
        %% plotting
        if plotting
            if ~exist('PhaseModulationFig','dir')
                mkdir('PhaseModulationFig');
            end
            h(end+1) = figure;
            hax = subplot(1,2,1);
            rose(spkphases{a})
            title(hax,['Cell #' num2str(a) '. Rayleigh p = ' num2str(phasestats.p(a)) '.'])
            
            hax = subplot(1,2,2);
            bar(phasebins*180/pi,phasedistros(:,a))
            xlim([0 360])
            set(hax,'XTick',[0 90 180 270 360])
            hold on;
            plot([0:360],cos(pi/180*[0:360])*0.05*max(phasedistros(:,a))+0.95*max(phasedistros(:,a)),'color',[.7 .7 .7])
            set(h(end),'name',['PhaseModPlotsForCell' num2str(a)]);
            print(fullfile('PhaseModulationFig',['PhaseModPlotsForCell' num2str(a)]),'-dpng','-r0');
        end
    end
end

%%
detectorName = 'bz_PhaseModulation';
channels = lfp.channels;
detectorParams = v2struct(intervals,samplingRate,method,plotting,numBins,...
    passband,powerThresh,channels);

PhaseLockingData = v2struct(phasedistros,phasebins,...
    phasestats,spkphases,...
    detectorName, detectorParams);
% try
% PhaseLockingData.region = spikes.region;
% catch
% PhaseLockingData.region = [];
% end
% PhaseLockingData.UID = spikes.UID;
% PhaseLockingData.sessionName = spikes.sessionName;

if saveMat
    save([lfp.Filename(1:end-4) '.PhaseLockingData.cellinfo.mat'],'PhaseLockingData');
end

end