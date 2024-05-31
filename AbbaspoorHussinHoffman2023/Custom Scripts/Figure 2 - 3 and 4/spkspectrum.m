function spk_LFP_Convol = spkspectrum(trialvec, timevec, varargin)

%%
p = inputParser;
addParameter(p,'frequency', 2:1:200,@isivector)           % Sampling Frequency [e.g. 1000]
addParameter(p,'LABELS',{'LFP'; 'Spike'},@iscell)         % Envelope Threshold [e.g. 2 5 ]
addParameter(p,'samplingFreq', 1000,@isnumeric)
parse(p,varargin{:})

% Assign Parameters
frequency                = p.Results.frequency;
LABELS                   = p.Results.LABELS;
samplingFreq             = p.Results.samplingFreq;

%%
tmpspkLFP = [];
tmpspkLFP.label   = LABELS;
tmpspkLFP.time    = timevec;
tmpspkLFP.trial   = trialvec;
tmpspkLFP.fsample = samplingFreq;

%%%%%%%%%%%% Spike - LFP Synchronization (ppc2)
cfg           = [];
cfg.method    = 'mtmconvol';
cfg.foi       = frequency;
cfg.t_ftimwin = 5./cfg.foi; % 5 cycles per frequency
cfg.taper     = 'hanning';
spk_LFP_Convol     = ft_spiketriggeredspectrum(cfg, tmpspkLFP);

end