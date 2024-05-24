%% Spike Cross Correlogram
%%% This code computes CCG between a pair of units.
%%% It also performes a shuffling statistical test for significance and
%%% performs multiple comparison correction

%%% Saman Abbaspoor 03/12/2024 - Hoffman Lab - Vanderbilt

clear all, clc
load('H:\The LaDy\Data\CellTable_Complete_8.mat')
SessionTable = unique(CellTable.date);
BPTable = readtable('Y:\CoDy\The LaDy\Sheets\Data_Processed - Copy.xlsx');

CellTypeOrder = [2010 31 2020 15 6 2030 0 3 18 10];

CellCCGsTable = CellTable(:, [1:5 8 42]);

for Sess = 1:size(SessionTable, 1)
    ['Session: ', num2str(Sess), '\', num2str(size(SessionTable, 1))]
    Session = SessionTable{Sess}
    
    findUnit = find(contains(CellTable.date, Session));
    
    sessUnits = CellTable.timestamps(findUnit);
    numUnits = size(sessUnits, 1);
    
    Indices = find(contains(string(datestr(BPTable.Date)), string(datestr(Session))));
    % Compute original CCGs
    binSize = 0.001; duration= 0.1;
    spkTimes = CellTable.timestamps(findUnit);
    [ccg,t] = CCG(spkTimes,[],'binSize', binSize, 'duration', duration, 'norm', 'count', 'Fs',1/30000);
    
    % Create surrogate dataset and compute CCG per surrogate data
    permutations = 5000;   
    binSize = 0.001; duration= 0.04;
    surccg = NaN([duration*1000+1 size(ccg, [2 3]) permutations]);

    parfor Perm = 1:permutations
        ISIs = cellfun(@(X) diff(X), spkTimes, 'UniformOutput', false);
        ShuffledISIs = cellfun(@(X) X(randperm(length(X))), ISIs, 'UniformOutput', false);
        surrogateTimes = cellfun(@(X, Y) cumsum([X(1); Y]), spkTimes, ShuffledISIs, 'UniformOutput', false);
        [surccg(:, :, :, Perm), ~] = CCG(surrogateTimes,[],'binSize', binSize, 'duration', duration,...
            'norm', 'count', 'Fs',1/30000);
    end
    
    % Compute Pvalues
    Pvals = (sum(surccg - ccg(51-20:51+20, :, :) > 0, 4)+1) ./ (permutations+1);

    % Multiple Comparison Correction per pairwise CCG
    h = NaN(size(Pvals));
    for ref = 1:length(findUnit)
        for tar = 1:length(findUnit)
            [h(:, ref, tar), crit_p, adj_ci_cvrg, adj_Pvals]=fdr_bh(Pvals(:, ref, tar), 0.01 ,'pdep');
%             [c_pvalues, c_alpha, h, extra] = fwer_holmbonf(Pvals(:, ref, tar), 0.05, true)
        end
    end
        
    
    for C = 1:length(findUnit)
        CellCCGsTable.CCG(findUnit(C))  = {squeeze(ccg(:, C, :))};
        CellCCGsTable.Pvals(findUnit(C))  = {squeeze(Pvals(:, C, :))};
        CellCCGsTable.h(findUnit(C))  = {squeeze(h(:, C, :))};
        CellCCGsTable.tbins(findUnit(C))  = {t};
    end
    
end

save(fullfile('H:\The LaDy\Data', 'CellCCGsTable'), 'CellCCGsTable',  '-v7.3')
