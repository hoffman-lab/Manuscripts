%% RippleSpikingContent
%%% This code calculates the spiking probability, firing rate, and FR
%%% ratio of units during ripples.

%%% Saman Abbaspoor 03/12/2024 - Hoffman Lab - Vanderbilt


function perpl_RippleSpikingContent()

clear all, clc
load('H:\The LaDy\Data\CellTable_Complete_6.mat');
SessionTable = unique(CellTable.date);
BPTable = readtable('Y:\CoDy\The LaDy\Sheets\Data_Processed - Copy.xlsx');

%%
RippleSpkTable = CellTable(:, [1:5 8 42]);
SpindleSpkTable = CellTable(:, [1:5 8 42]);
SOSpkTable = CellTable(:, [1:5 8 42]);
DeltaSpkTable = CellTable(:, [1:5 8 42]);

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
    
    for U = 1:numUnits
        PSTH = NaN(size(ripple_sn, 1)-2, 2001);
        fs = min(find(spikemat.data(:, U) > 0));
        ls = max(find(spikemat.data(:, U) > 0));
        for Event = 1:size(ripple_sn, 1)-2
            if ripple_sn(Event) > fs & ripple_sn(Event) < ls
                PSTH(Event, :) = [spikemat.data(ripple_sn(Event)-1000:ripple_sn(Event)+1000,U)]';
            end
        end
        to_keep = find(~isnan(sum(PSTH')));
        RippleSpkTable.ripspkmat(findUnit(U)) = {logical(PSTH(to_keep, :))};
    end
    
    %%%% Create Ripple-aligned Spike Matrix
    % $$$$$$$$ check Within_Ripple_Duration and ripple durations from ripple table
    ripple_sn = interp1(spikemat.timestamps, 1:length(spikemat.timestamps), Ripple_Table.peak_ts, 'nearest');
    ripple_start = interp1(spikemat.timestamps, 1:length(spikemat.timestamps), Ripple_Table.start_ts, 'nearest');
    ripple_stop  = interp1(spikemat.timestamps, 1:length(spikemat.timestamps), Ripple_Table.stop_ts, 'nearest');
    

    for U = 1:numUnits
        Within_Ripple_Spikes = [];
        Within_Ripple_Duration = [];
        fs = min(find(spikemat.data(:, U) > 0));
        ls = max(find(spikemat.data(:, U) > 0));
        for Event = 1:size(ripple_start, 1)-2
            if ripple_sn(Event) > fs & ripple_sn(Event) < ls
                tmp = (spikemat.data(ripple_start(Event):ripple_stop(Event),U));
                Within_Ripple_Spikes = [Within_Ripple_Spikes, sum(tmp)'];
                Within_Ripple_Duration = [Within_Ripple_Duration, length(tmp)];
            end
        end
        Within_Ripple_FR = sum(Within_Ripple_Spikes') ./ (sum(Within_Ripple_Duration)/1000);
        RippleSpkTable.Ripple_FR(findUnit(U))     = Within_Ripple_FR;
        RippleSpkTable.Ripple_Spknum(findUnit(U)) = sum(Within_Ripple_Spikes);
        RippleSpkTable.Ripple_Ratio(findUnit(U))  = Within_Ripple_FR ./ CellTable.firingRate(findUnit(U));

        Ripple_Participation = sum(Within_Ripple_Spikes>0)/numel(Within_Ripple_Spikes);
        rip_active_num       = sum(Within_Ripple_Spikes>0);
        
        RippleSpkTable.Ripple_num(findUnit(U))              = size(Ripple_Table, 1);
        RippleSpkTable.Ripple_active_num(findUnit(U))       = rip_active_num;
        RippleSpkTable.Ripple_Participation_Probability(findUnit(U)) = Ripple_Participation;
        RippleSpkTable.Within_Ripple_Spikes(findUnit(U))    = {Within_Ripple_Spikes};
        RippleSpkTable.Within_Ripple_Duration(findUnit(U))  = {Within_Ripple_Duration};
    end
    
    
    save(fullfile('H:\The LaDy\Data', 'RippleSpkTable'), 'RippleSpkTable',  '-v7.3')

end

end

