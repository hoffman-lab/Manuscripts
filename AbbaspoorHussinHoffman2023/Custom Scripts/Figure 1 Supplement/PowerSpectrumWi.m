%% Load Data
filedir = 'D:\HiTGun\Wi TH Data\Processed';
Sessions = struct2table(dir(filedir));
Sessions = Sessions(3:end, :);

PowerSpectralDensity = table;
for sess = 1:size(Sessions, 1) 
    
Session = Sessions.name{sess};
    
    
%% Compute Power Spectrum - REST

ov       = 0.5;   % overlap
winlen   = 1024; % winlen
nff      = 400;
freqBins = linspace(1, 150, nff); %logspace(log10(1), log10(200), nff);

PSD = NaN(2, nff);
% for ep = 1:2
%     load(fullfile(filedir, Sessions.name{sess}, ['rest', num2str(ep)]))
%     Rest = removeLineNoise_SpectrumEstimation(lfp(2, :), 1000, 'NH = 4, LF = 60, M = 1024');
%     [pd, freqs] = pwelch(Rest,winlen,floor(ov*winlen), freqBins, 1000);
%     PSD(ep, :) = pd;
% end

load(fullfile(filedir, Sessions.name{sess}, 'Sleep'))
% Rest = removeLineNoise_SpectrumEstimation(lfp(2, :), 1000, 'NH = 4, LF = 60, M = 1024');
[pd, freqs] = pwelch(lfp,winlen,floor(ov*winlen), freqBins, 1000);
PowerSpectralDensity.session(sess) = {Session};
PowerSpectralDensity.Rest(sess) = {pd};


%% Compute Power Spectrum - Search

load(fullfile(filedir, Sessions.name{sess}, 'Treehouse'))
to_remove = cellfun(@(x) isempty(x), TrialLFP);
TrialLFP(to_remove) = []; TrialTime(to_remove) = [];

PSD = NaN(length(TrialLFP), nff);
for Trial = 1:length(TrialLFP)
    [pd, freqs] = pwelch(TrialLFP{1, Trial},winlen,floor(ov*winlen), freqBins, 1000);
    PSD(Trial, :) = pd;
end

PSD = mean(PSD, 1); %PSD = zscore(pow2db(PSD));
PowerSpectralDensity.Search(sess) = {PSD};


%%

% fh = figure();
% fh.WindowState = 'maximized';
% 
% Col = [255 8 32]/255;
% plot(freqs, PowerSpectralDensity.Search{sess}, 'LineWidth', 2, 'Color', Col)
% hold on
% Col = [56 138 191]/255;
% plot(freqs, PowerSpectralDensity.Rest{sess}, 'LineWidth', 2, 'Color', Col)
% 
% axis xy
% xlabel('Frequency'); ylabel('Power (dB)')
% pbaspect([1 1 1])
% set(gca, 'box', 'off', 'FontSize', 15, 'tickdir', 'out', 'LineWidth', 1.5,...
%     'xtick', [4 8 20 50 100 150])
% 
% title(Session)

%%
% figName = [Session]
Directory = 'D:\HiTGun\Wi TH Data\Analyses\Power Spectral Density'
% saveas(gcf, fullfile(Directory, [figName '.png']))
% close(gcf)


end

WiPSD = PowerSpectralDensity;
save(fullfile(Directory, 'PowerSpectralDensity_raw'), 'PowerSpectralDensity', 'freqs', '-v7.3')

