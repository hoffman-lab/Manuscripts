function perpl_CreateRippleSurrogate()

clear all, clc
load('H:\The LaDy\Data\CellTable_Complete_8.mat');
SessionTable = unique(CellTable.date);
BPTable = readtable('Y:\CoDy\The LaDy\Sheets\Data_Processed - Copy.xlsx');

%%
load(fullfile('H:\The LaDy\Data', 'SurrogateRippleSpkTable'))

for Sess = 1:size(SessionTable, 1)
    
    %     if Sess == 29
    %         continue
    %     end
    
    ['Session: ', num2str(Sess), '\', num2str(size(SessionTable, 1))]
    Session = SessionTable{Sess}
    
    findUnit = find(contains(CellTable.date, Session));
    
    sessUnits = CellTable.timestamps(findUnit);
    numUnits = size(sessUnits, 1);
    
    Indices = find(contains(string(datestr(BPTable.Date)), string(datestr(Session))));
    
    %%
    if contains(CellTable.AnimalID(findUnit(1)), 'WI')
        %         nTotalChannels = 128;
        %         Selected_Channels = [BPTable.Ripple_Chan_TH(Indices)];
        riptmpdir = BPTable.directory{Indices};
        riptmpdir = strrep(riptmpdir,'Y:\CoDy','G:');
        riptablename = 'Ripple_Table.mat';
        load(fullfile(riptmpdir, riptablename))
        load(fullfile(riptmpdir, 'timestamps.mat'));
        
    elseif contains(CellTable.AnimalID(findUnit(1)), 'FN')
        %         nTotalChannels = 64;
        %         Selected_Channels = [BPTable.Ripple_Chan_TH(Indices)];
        riptmpdir = BPTable.directory{Indices};
        riptmpdir = strrep(riptmpdir,'Y:\CoDy','H:');
        riptablename = 'Ripple_Table.mat';
        load(fullfile(riptmpdir, riptablename))
        load(fullfile(riptmpdir, 'timestamps.mat'))
    end
    
    Confirmed = Ripple_Table.ConfirmedRipples;
    Ripple_Table = Ripple_Table(logical(Confirmed), :);
    
    
    %% Ripple
    %%% Create Binary Spike Matrix
    UnitActivity = CellTable.timestamps(findUnit); UnitActivity = UnitActivity';
    [spikemat] = bz_SpktToSpkmat(UnitActivity, 'win', [timestamps(1) timestamps(end)],...
        'dt', 0.001, 'binsize', 0.001,...
        'bintype', 'boxcar',...
        'units', 'counts');
    
    %%%%% Create Ripple-aligned Spike Matrix
    ripple_sn = interp1(spikemat.timestamps, 1:length(spikemat.timestamps), Ripple_Table.peak_ts, 'nearest');
    ripple_sn = select_random_ts(ripple_sn, length(spikemat.timestamps));
    
    for U = 1:numUnits
        PSTH = NaN(numel(ripple_sn)-2, 2001);
        fs = min(find(spikemat.data(:, U) > 0));
        ls = max(find(spikemat.data(:, U) > 0));
        for Event = 1:numel(ripple_sn)-2
            if ripple_sn(Event) > fs & ripple_sn(Event) < ls
                PSTH(Event, :) = [spikemat.data(ripple_sn(Event)-1000:ripple_sn(Event)+1000,U)]';
            end
        end
        to_keep = find(~isnan(sum(PSTH')));
        SurrogateRippleSpkTable.ripspkmat(findUnit(U)) = {logical(PSTH(to_keep, :))};
    end
    
    
    save(fullfile('H:\The LaDy\Data', 'SurrogateRippleSpkTable'), 'SurrogateRippleSpkTable',  '-v7.3')
    
end

end


function rand_sn = select_random_ts(event_sn, npoints)
NUMEVENTS = numel(event_sn);
tmpts = ones(npoints, 1);

for i = 1:size(event_sn, 1)
    try
        tmpts([event_sn(i)-2000:event_sn(i)+2000]) = 0;
    catch
    end
end


rand_sn = [];
for i = 1:size(event_sn, 1)
    try
        tmpid = [event_sn(i)-60000:event_sn(i)-10000 event_sn(i)+10000:event_sn(i)+60000];
        y = randsample(find(tmpts(tmpid) == 1),1);
        rand_sn = [rand_sn, tmpid(y)];
    catch
    end
end

end
