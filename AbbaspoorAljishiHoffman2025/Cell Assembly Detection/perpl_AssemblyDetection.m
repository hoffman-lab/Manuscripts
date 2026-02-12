%% Assembly During Sleep
%%% Ref: V Lopes-dos-Santos 2013; Detecting cell assemblies in large neuronal populations
%%% To run the code you need Cell Assembly Detection toolbox from:
%%% https://github.com/tortlab/Cell-Assembly-Detection
%%% This code runs PCA on the normalized spike matrix, and comoutes the
%%% number of significant components/cell assemblies using shuffling or
%%% analytical methods. Next, it runs ICA on the PCA-reduced space to demix
%%% the signal and extract cell assemblies.

%%% Saman Abbaspoor, 3/12/2024 - Hoffman Lab - Vanderbilt

%%
function perpl_AssemblyDetection()

clear all, clc
load('H:\The LaDy\Data\CellTable_Complete_8.mat');
SessionTable = unique(CellTable.date);
BPTable = readtable('Y:\CoDy\The LaDy\Sheets\Data_Processed - Copy.xlsx');

%%
Assembly_Table = table;

for Sess = 1:size(SessionTable, 1)
    
    ['Session: ', num2str(Sess), '\', num2str(size(SessionTable, 1))]
    Session = SessionTable{Sess}
    
    findUnit = find(contains(CellTable.date, Session));
    
    sessUnits = CellTable.timestamps(findUnit);
    numUnits = size(sessUnits, 1);
    
    Indices = find(contains(string(datestr(BPTable.Date)), string(datestr(Session))));
    
    %%
    if contains(CellTable.AnimalID(findUnit(1)), 'WI')
        tmpdir = BPTable.directory{Indices};
        tmpdir = strrep(tmpdir,'Y:\CoDy','G:');
        load(fullfile(tmpdir, 'timestamps.mat'))
        
        %         riptablename = 'Ripple_Table.mat';
        %         load(fullfile(tmpdir, riptablename))
        
    elseif contains(CellTable.AnimalID(findUnit(1)), 'FN')
        tmpdir = BPTable.directory{Indices};
        tmpdir = strrep(tmpdir,'Y:\CoDy','H:');
        load(fullfile(tmpdir, 'timestamps.mat'))
        
        %         riptablename = 'Ripple_Table.mat';
        %         load(fullfile(tmpdir, riptablename))
    end
    
    %     Confirmed = Ripple_Table.ConfirmedRipples;
    %     Ripple_Table = Ripple_Table(logical(Confirmed), :);
    %
    %% Create the Activity Spike Matrix
    TH_end = CellTable.TH_end; TH_end = cell2mat(TH_end);
    TH_end_sample = TH_end(findUnit, 1); TH_end_sample = num2cell(TH_end_sample)';
    TH_end = TH_end(findUnit, 2); TH_end = num2cell(TH_end)';
    
    UnitActivity = CellTable.timestamps(findUnit); UnitActivity = UnitActivity';
    
    UnitActivity = cellfun(@(X, Y) X(X>Y), UnitActivity, TH_end, 'UniformOutput', false);
    EndTime = max(cell2mat(UnitActivity'));
    
    fprintf('Creating Spike Matrix ... \n')
    dt = 0.001; binsize = 0.09; bintype = 'gaussian'; units = 'counts';
    [spikemat] = bz_SpktToSpkmat(UnitActivity, 'win', [timestamps(TH_end_sample{1}) EndTime],... %timestamps(end)
        'dt', dt, 'binsize', binsize,...
        'bintype', bintype,...
        'units', units);
    
    Activitymatrix = spikemat.data';
    
    % Activitymatrix = smoothdata(Activitymatrix, 'gaussian', 10);
    % Activitymatrix = zscore(Activitymatrix, [], 2);
    
    %%
    fprintf('Computing PCs ... \n')
    
    opts = [];
    opts.threshold.permutations_percentile = 95;
    % opts.threshold.number_of_permutations = 20;
    opts.threshold.method = 'MarcenkoPastur';
    opts.Patterns.method  = 'PCA';
    AssemblyTemplates_PCA = assembly_patterns(Activitymatrix,opts);
    
    %% Apply ICA to the PCA subspace
    fprintf('Apply ICA ... \n')
    
    Zproj = AssemblyTemplates_PCA'*Activitymatrix;
    
    number_of_iterations = 1000; %;
    AssemblyTemplates_ICA=...
        fast_ica(Zproj,size(AssemblyTemplates_PCA, 2),number_of_iterations);
    
    
    AssemblyTemplates = AssemblyTemplates_PCA*AssemblyTemplates_ICA;
    
    %% Sign Correction
    for i = 1:size(AssemblyTemplates, 2)
        [~, I] = max(abs(AssemblyTemplates(:,i)));
        if AssemblyTemplates(I,i) < 0
            AssemblyTemplates(:,i) = -1.*AssemblyTemplates(:,i);
        end
    end
    
    %% Detect significant assembly members
    SigWeights = zeros(size(AssemblyTemplates));
    
    for Ass = 1:size(AssemblyTemplates, 2)
        tmpass = abs(AssemblyTemplates(:,Ass));
        MEAN = mean(tmpass); STD = std(tmpass); level = MEAN + 2*STD;
        ixx = abs(AssemblyTemplates(:,Ass)) > level;
        SigWeights(ixx, Ass) = 1;
    end
    
    
    %%
    Assembly_Table.date(Sess)          = CellTable.date(findUnit(1));
    Assembly_Table.AnimalID(Sess)      = CellTable.AnimalID(findUnit(1));
    Assembly_Table.BrainRegion(Sess)   = CellTable.BrainRegion(findUnit(1));
    
    Assembly_Table.AssemblyTemplates(Sess) = {AssemblyTemplates};
    Assembly_Table.SigWeights(Sess)        = {SigWeights};
    
    Assembly_Table.putativeCellType(Sess) = {CellTable.putativeCellType(findUnit)};
    Assembly_Table.SupDeep(Sess)         = {CellTable.SupDeep(findUnit)};
    
    Assembly_Table.Ref(Sess)   = {'Sleep'};
    Assembly_Table.dt(Sess)      = dt;
    Assembly_Table.binsize(Sess) = binsize;
    Assembly_Table.bintype(Sess) = {bintype};
    Assembly_Table.units(Sess)   = {units};
    
    %%
    save(fullfile('H:\The LaDy\Data', 'Assembly_Table_RefSleep_90ms_gauss'), 'Assembly_Table',  '-v7.3')
end

end

