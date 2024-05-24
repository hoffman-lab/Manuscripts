function perpl_CreateSurrogateActivity()

clear all, clc
load('H:\The LaDy\Data\CellTable_Complete_8.mat');
SessionTable = unique(CellTable.date);
BPTable = readtable('Y:\CoDy\The LaDy\Sheets\Data_Processed - Copy.xlsx');

%%
SurrogateActivityTable = CellTable(:, [1:5 8 42]);

for Sess = 1:size(SessionTable, 1)

    
    ['Session: ', num2str(Sess), '\', num2str(size(SessionTable, 1))]
    Session = SessionTable{Sess}
    
    findUnit = find(contains(CellTable.date, Session));
    
    sessUnits = CellTable.timestamps(findUnit);
    numUnits = size(sessUnits, 1);
    
    Indices = find(contains(string(datestr(BPTable.Date)), string(datestr(Session))));
    
    %%
    if contains(CellTable.AnimalID(findUnit(1)), 'WI')
        nTotalChannels = 128;
        Selected_Channels = [BPTable.Ripple_Chan_TH(Indices)];
        riptmpdir = BPTable.directory{Indices};
        riptmpdir = strrep(riptmpdir,'Y:\CoDy','G:');
        
        load(fullfile(riptmpdir, 'timestamps.mat'))
    elseif contains(CellTable.AnimalID(findUnit(1)), 'FN')
        nTotalChannels = 64;
        Selected_Channels = [BPTable.Ripple_Chan_TH(Indices)];
        riptmpdir = BPTable.directory{Indices};
        riptmpdir = strrep(riptmpdir,'Y:\CoDy','H:');
        
        load(fullfile(riptmpdir, 'timestamps.mat'))
    end
    
    %% Ripple
    %%% Create Binary Spike Matrix
    UnitActivity = CellTable.timestamps(findUnit); UnitActivity = UnitActivity';
    [spikemat] = bz_SpktToSpkmat(UnitActivity, 'win', [timestamps(1) timestamps(end)],...
        'dt', 0.001, 'binsize', 0.001,...
        'bintype', 'boxcar',...
        'units', 'counts');
    
    
    %%%%% Create Ripple-aligned Spike Matrix
%     ts = timestamps(datasample(1:numel(timestamps), 1000, 'Replace', false));
    ts = timestamps(randsample(numel(timestamps), 1000));
    ripple_sn = interp1(spikemat.timestamps, 1:length(spikemat.timestamps), ts, 'nearest');
    ripple_sn(ripple_sn - 1000 < 0) = ripple_sn(ripple_sn - 1000 < 0) + 1000;
    ripple_sn(ripple_sn + 1000 > length(spikemat.data)) = ripple_sn(ripple_sn + 1000 > length(spikemat.data)) - 1000;

    PopVec = NaN(numUnits, size(ripple_sn, 1), 2001);
    for Event = 1:size(ripple_sn, 1)
        PopVec(:, Event, :) = [spikemat.data(ripple_sn(Event)-1000:ripple_sn(Event)+1000,:)]';
    end
    PopVec = logical(PopVec);
    for c = 1:numUnits
        SurrogateActivityTable.ripspkmat(findUnit(c)) = {squeeze(PopVec(c, :, :))};
    end
end

save(fullfile('H:\The LaDy\Data', 'SurrogateActivityTable'), 'SurrogateActivityTable',  '-v7.3')


end

