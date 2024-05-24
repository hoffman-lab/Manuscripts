%% Pairwise Phase Consistency
%%% This script computed PPC0 independent of FieldTrip.
%%% It first TF decompose the signal using Morlet Wavelets and then rus PPC
%%% on each frequency.

%%% Saman Abbaspoor 03/12/2024 - Hoffman Lab - Vanderbilt



function perpl_PairwisePhaseConsistency
clear all, clc
load('H:\The LaDy\Data\CellTable.mat');

SessionTable = unique(CellTable.date);

PPCTable = CellTable(:, [1:5 8]);

BPTable = readtable('Y:\CoDy\The LaDy\Sheets\Data_Processed - Copy.xlsx');

for Sess = 1:size(SessionTable, 1)
    %% Load Unit Data
    Session = SessionTable{Sess};
    
    findUnit = find(contains(CellTable.date, Session));
    
    sessUnits = CellTable.timestamps(findUnit);
    numUnits = size(sessUnits, 1);
    
    Indices = find(contains(string(datestr(BPTable.Date)), string(datestr(Session))));
    
    %% Load LFP Data
    if contains(CellTable.AnimalID(findUnit(1)), 'WI')
        nTotalChannels = 128;
        Selected_Channels = [BPTable.Ripple_Chan_TH(Indices)];
        tmpdir = BPTable.directory{Indices};
        tmpdir = strrep(tmpdir,'Y:\CoDy','G:');
    elseif contains(CellTable.AnimalID(findUnit(1)), 'FN')
        nTotalChannels = 64;
        Selected_Channels = [BPTable.Ripple_Chan_TH(Indices)];
        tmpdir = BPTable.directory{Indices};
        tmpdir = strrep(tmpdir,'Y:\CoDy','H:');
    end
    
    data = perpl_LoadBinary(fullfile(tmpdir,'aHPC_B_cnct.dat'),...
        'frequency', 30000,...
        'offset', 0,...
        'samples', inf,...
        'nChannels', nTotalChannels,...
        'channels', [Selected_Channels],...
        'downsample', 1,...
        'bitVolts', 0.195);
    
    [b a] = butter(3, 350/30000);
    RipChan = FiltFiltM(b,a,data', 2); RipChan = downsample(RipChan', 30)';
    
    load(fullfile(tmpdir,'timestamps.mat')); timestamps = downsample(timestamps', 30)';
    
    lfp.data = RipChan';
    lfp.timestamps = timestamps;
    lfp.samplingRate = 1000;
    [wavespec] = MorletWaveSpec(lfp, 'frange', [4 200],...
        'nfreqs', 98,...
        'ncyc', 7);
    
    spksample = cellfun(@(X) interp1(timestamps, 1:numel(timestamps), X, 'nearest'),...
        sessUnits, 'UniformOutput', false);
    
    tic
    PPCstruct = [];
    for c = 1:numUnits
        spkPhases = angle(wavespec.data(spksample{c, 1}, :));
        ppc0 = zeros([1 98]);
        parfor i = 1:98
            ppc0(i) = ppc(spkPhases(:, i));
        end
        PPCstruct.Val{c} = ppc0;
        PPCstruct.freq{c} = wavespec.freqs;
    end
    toc
    
    tmp = [];
    [tmp{1:numUnits}]            = deal(tmpdir);
    PPCTable.directory(findUnit) = tmp;
    
    PPCTable.PPCval(findUnit)    = PPCstruct.Val;
    PPCTable.PPCfreq(findUnit)   = PPCstruct.freq;
    
end

save(fullfile('H:\The LaDy\Data\FN', 'PPCTable'), 'PPCTable',  '-v7.3')

end
