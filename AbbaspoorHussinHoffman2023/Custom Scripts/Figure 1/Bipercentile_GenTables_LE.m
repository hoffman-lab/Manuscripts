% clear all, clc
addpath 'D:\Matlab Packages\fieldtrip-lite-20190224'


FilesDir = 'D:\HiTGun\Data\decimatedCSC_wholeSession';
TTLFilesDir = 'D:\HiTGun\Data\LE Change blindness Raw';

load('D:\HiTGun\Data\LETrialsTable');

info = struct2table(dir(FilesDir));
SessionTable = info.name;
SessionTable = SessionTable(3:end);

Row = 1;

biper_si_LE = cell([size(SessionTable, 1) 1]);

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
    channellist = [15];   %inchannel 13 14
    channellist = unique(channellist);
    
    for ch = 1%:length(channellist)
        
        channel = {['CSC', num2str(channellist(ch))]};
        
        %%% Load Channel LFP
        filedir = fullfile(FilesDir, Session, channel);
        load(filedir{1})
        
        %% Find Event of Interest
        Trials{1} = decimateCSC;
        tmp = table;
        tmp.trlfp = Trials';
        [Type{1:size(tmp,1)}] = deal('Task');
        tmp.type = Type';
        TaskTable = tmp;
        clear Type tmp
        
        %% Compute FFT
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
        biper_si_LE{Sess} = MSortedPSD;
    end
end
