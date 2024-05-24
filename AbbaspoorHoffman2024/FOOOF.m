%% FOOOF
%%% The code requires FOOOOF Matlab wrapper


BPTable = readtable('Y:\CoDy\The LaDy\Sheets\Data_Processed - Copy.xlsx');
% Sessions = find(BPTable.RSC_Concatenated == 1 & BPTable.Ripple_detected==0);
% BPTable.directory = cellfun(@(X) strrep(X,'Y:\', 'I:\'), BPTable.directory, 'UniformOutput', false);

%%
ov       = 0.5;   % overlap
winlen   = 1024; % winlen
nff      = 96*4;
freqBins = linspace(4, 100, nff); %logspace(log10(1), log10(200), nff);

% FOOOF settings
settings                = struct();  % Use defaults
settings.aperiodic_mode = 'knee';
settings.verbose        = 0;

f_range = [4, 100];

[b a] = butter(3, 350/15000);

%%

for Sess = 1:35
    if contains(BPTable.Animal_ID{Sess}, 'FN')
        nTotalChannels = 64;
        channels = 1:64;
        tmpdir = BPTable.directory{Sess};
        tmpdir = strrep(tmpdir,'Y:\CoDy\', 'H:\');
    elseif contains(BPTable.Animal_ID{Sess}, 'Wi')
        nTotalChannels = 128;
        channels = 1:2:128;
        tmpdir = BPTable.directory{Sess};
        tmpdir = strrep(tmpdir,'Y:\CoDy\', 'G:\');
    end
    
    
    display(['Processing Session:', tmpdir])
    
    TH_end = BPTable.TH_end{Sess};
    TH_end = str2num(TH_end); TH_end = TH_end(1);
    
    %% Load files
    data = perpl_LoadBinary(fullfile(tmpdir, 'aHPC_B_cnct.dat'),...
        'frequency', 30000,...
        'offset', 0,...
        'samples', 30*60*30000,...
        'nChannels', nTotalChannels,...
        'channels', channels,... %[BPTable.Ripple_Channel(Sess) 40]
        'downsample', 1,...
        'bitVolts', 0.195);
    
    THLFP = FiltFiltM(b,a,data', 2); THLFP = downsample(THLFP', 30)';
    %     THLFP = removeLineNoise_SpectrumEstimation(THLFP, 1000, 'NH = 3, LF = 60, M = 1024');
    
    try
        data = perpl_LoadBinary(fullfile(tmpdir, 'aHPC_B_cnct.dat'),...
            'frequency', 30000,...
            'offset', TH_end+30*60*30000,...
            'samples', 30*60*30000,...
            'nChannels', nTotalChannels,...
            'channels', channels,... %[BPTable.Ripple_Channel(Sess) 40]
            'downsample', 1,...
            'bitVolts', 0.195);
    catch
        data = perpl_LoadBinary(fullfile(tmpdir, 'aHPC_B_cnct.dat'),...
            'frequency', 30000,...
            'offset', TH_end,...
            'samples', 30*60*30000,...
            'nChannels', nTotalChannels,...
            'channels', channels,... %[BPTable.Ripple_Channel(Sess) 40]
            'downsample', 1,...
            'bitVolts', 0.195);
    end
    
    SleepLFP = FiltFiltM(b,a,data', 2); SleepLFP = downsample(SleepLFP', 30)';
    %     SleepLFP = removeLineNoise_SpectrumEstimation(SleepLFP, 1000, 'NH = 3, LF = 60, M = 1024');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    THFOOOF = table;
    for Ch = 1:size(THLFP, 1)
        [psd, freqs] = pwelch(THLFP(Ch, :),winlen,floor(ov*winlen), freqBins, 1000);
        % Transpose, to make inputs row vectors
        freqs = double(freqs'); psd = double(psd');
        % Run FOOOF, also returning the model
        fooof_results = fooof(freqs, psd, f_range, settings, true);
        psd = fooof_results.fooofed_spectrum - fooof_results.ap_fit;
        THFOOOF.freqs(Ch)             = {fooof_results.freqs};
        THFOOOF.power_spectrum(Ch)    = {fooof_results.power_spectrum};
        THFOOOF.fooofed_spectrum(Ch)  = {fooof_results.fooofed_spectrum};
        THFOOOF.ap_fit(Ch)            = {fooof_results.ap_fit};
        THFOOOF.AdjPower(Ch)          = {psd};
    end
    
    SleepFOOOF = table;
    for Ch = 1:size(SleepLFP, 1)
        [psd, freqs] = pwelch(SleepLFP(Ch, :),winlen,floor(ov*winlen), freqBins, 1000);
        % Transpose, to make inputs row vectors
        freqs = double(freqs'); psd = double(psd');
        % Run FOOOF, also returning the model
        fooof_results = fooof(freqs, psd, f_range, settings, true);
        psd = fooof_results.fooofed_spectrum - fooof_results.ap_fit;
        SleepFOOOF.freqs(Ch)             = {fooof_results.freqs};
        SleepFOOOF.power_spectrum(Ch)    = {fooof_results.power_spectrum};
        SleepFOOOF.fooofed_spectrum(Ch)  = {fooof_results.fooofed_spectrum};
        SleepFOOOF.ap_fit(Ch)            = {fooof_results.ap_fit};
        SleepFOOOF.AdjPower(Ch)          = {psd};
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    save(fullfile(tmpdir, 'FOOOF.mat'), 'THFOOOF', 'SleepFOOOF', '-v7.3');
    
end
