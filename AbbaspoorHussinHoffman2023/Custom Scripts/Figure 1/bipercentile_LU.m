% clear all, clc

load(fullfile('D:\HiTGun\Data','LU_In_Layer_Channel.mat'))
load('D:\HiTGun\Data\LUTrialsTable');



s = {'LU_2014-11-27_11-04-38', ...
    'LU_2014-11-29_13-14-04', ...
    'LU_2014-12-02_11-09-10', ...
    'LU_2014-12-03_11-26-53', ...
    'LU_2014-12-04_11-09-08', ...
    'LU_2014-12-05_11-33-27', ...
    'LU_2014-12-06_11-27-34', ...
    'LU_2014-12-10_11-27-48', ...
    'LU_2014-12-16_11-36-02', ...
    'LU_2014-12-18_11-28-38', ...
    'LU_2014-12-23_11-19-55', ...
    'LU_2014-12-24_11-25-36', ...
    'LU_2014-12-26_11-13-14', ...
    'LU_2014-12-29_11-26-47', ...
    'LU_2014-12-31_11-17-20', ...
    'LU_2015-01-05_11-05-18', ...
    'LU_2015-01-06_11-23-22', ...
    'LU_2015-01-07_11-01-23', ...
    'LU_2015-01-08_10-28-27', ...
    'LU_2015-01-09_11-16-46', ...
    'LU_2015-01-12_12-18-20', ...
    'LU_2015-01-13_11-05-52', ...
    'LU_2015-01-14_11-09-15', ...
    'LU_2015-01-15_11-12-18'};

restIdxDir = 'D:\HiTGun\restIdx'
FilesDir = 'D:\HiTGun\Data\LU Change blindness decimated';

info = struct2table(dir(FilesDir));
SessionTable = info.name;
SessionTable = SessionTable(3:end);

ChannelList = [29 30 31];

biper_si = cell([size(SessionTable, 1) 1]);
biper_Rest   = cell([size(SessionTable, 1) 1]);

for Sess = 1:size(SessionTable, 1)
    
    if Sess == 8 || Sess == 17 || Sess == 28; continue; end
    
    %% Load Data
    Session = SessionTable{Sess};
    
    ch = 1;
    channel = {['csc', num2str(ChannelList(ch))]};
    
    info = struct2table(dir(fullfile(FilesDir, Session)));
    id = find(contains(info.name, channel));
    if isempty(id); continue; end
    channel = info.name(id);
    channel = {channel{1}(1:end-4)};
    
    %%% Load Channel LFP
    load(fullfile(FilesDir, Session, channel{1}))
    
    %% Reduce Line Noise
    downSampCSC = removeLineNoise_SpectrumEstimation(downSampCSC', 1000, 'NH = 1, LF = 60, M = 1024');
    downSampCSC = downSampCSC';
    
    
    %% Find Event of Interest
    
    Trials{1} = downSampCSC;
    tmp = table;
    tmp.trlfp = Trials';
    [Type{1:size(tmp,1)}] = deal('Task');
    tmp.type = Type';
    TaskTable = tmp;
    clear Type tmp
    
    %% Remove trials shorter than 1 seconds
    
    SIZE = cellfun(@(x) length(x), TaskTable.trlfp);
    IDX = find(SIZE < 2048);
    TaskTable(IDX, :) = [];

    
    %% Compute Power-Power Correlation
    Epochs = {'Task'};
    
    Index = find(contains(TaskTable.type,Epochs{1}));
    
    ov = 0.5;   % overlap
    winlen = 1024; % winlen
    nff = 150;
    freqBins = linspace(1, 150, nff); %logspace(log10(1), log10(200), nff);
    
    
    PSD = cell(length(Index), 1);
    for Tr = 1:length(Index)
        [s,freq,time, pd] = spectrogram(TaskTable.trlfp{Index(Tr)},winlen,floor(ov*winlen), freqBins, 1000);
        PSD{Tr} = pd;
    end
    
    trialBins{1} = cellfun(@(x) size(x, 2), PSD);
    PSD = {cat(2, PSD{:})}; PSD = PSD{1}; %specgr{ep} = PSD;
    
    [c id1] = min(abs(freq - 20));
    [c id2] = min(abs(freq - 30));
    sd = PSD(id1:id2, :); sd = mean(sd, 1);
    [~, id] = sort(sd, 'ascend');
    SortedPSD = PSD(:, id);
    SortedPSD = SortedPSD./median(SortedPSD, 2);
    
    if Sess == 1;
        save('SortedIDSess1LE', 'id')
    end
    
    
    MSortedPSD = NaN(length(freq), 50);
    batch = floor(linspace(1, size(SortedPSD, 2), 51));
    for b = 1:50
        MSortedPSD(:, b) = mean(SortedPSD(:, batch(b):batch(b+1)), 2)';
    end
    biper_si{Sess} = MSortedPSD;
end
