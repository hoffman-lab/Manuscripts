clear all, clc

load(fullfile('D:\HiTGun\Data','LU_In_Layer_Channel.mat'))
load('D:\HiTGun\Data\LUTrialsTable');

SessionTable = {'LU_2014-11-27_11-04-38', ...
    'LU_2014-12-02_11-09-10', ...
    'LU_2014-12-03_11-26-53', ...
    'LU_2014-12-05_11-33-27', ...
    'LU_2014-12-06_11-27-34', ...
    'LU_2014-12-10_11-27-48', ...
    'LU_2014-12-16_11-36-02', ...
    'LU_2014-12-18_11-28-38', ...
    'LU_2014-12-23_11-19-55', ...
    'LU_2014-12-26_11-13-14', ...
    'LU_2014-12-31_11-17-20', ...
    'LU_2015-01-05_11-05-18', ...
    'LU_2015-01-08_10-28-27', ...
    'LU_2015-01-09_11-16-46', ...
    'LU_2015-01-12_12-18-20', ...
    'LU_2015-01-13_11-05-52', ...
    'LU_2015-01-14_11-09-15', ...
    'LU_2015-01-15_11-12-18'}';

restIdxDir = 'D:\HiTGun\restIdx'
FilesDir = 'D:\HiTGun\Data\LU Change blindness decimated';

% ChannelList = [6 10 17 21 29 30 31];

PwelchLU = table;
Row = 1;
for Sess = 1:size(SessionTable, 1)
    if Sess == 8 || Sess == 17 || Sess == 28; continue; end
    %% Load Data
    Session = SessionTable{Sess};
    inchannel = LU_In_Layer_Channel.rip_cscNum(find(contains(string(LU_In_Layer_Channel.sid), Session)));
    ChannelList = [29];
    ch = 1
    
    channel = {['csc', num2str(ChannelList(ch))]};
    
    info = struct2table(dir(fullfile(FilesDir, Session)));
    id = find(contains(info.name, channel));
    if isempty(id); continue; end
    channel = info.name(id);
    channel = {channel{1}(1:end-4)};
    
    %%% Load Channel LFP
    try
        load(fullfile(FilesDir, Session, channel{1}))
    catch
        continue
    end
    
    %% Reduce Line Noise
    downSampCSC = removeLineNoise_SpectrumEstimation(downSampCSC', 1000, 'NH = 1, LF = 60, M = 1024');
    downSampCSC = downSampCSC';
    
    %% Find Event of Interest
    
    
    %%%%%%%%%%%%Search
    id = find( contains(LUTrialsTable.session,Session(1:13)) &...
        contains(LUTrialsTable.type, 'Search')&...
        LUTrialsTable.performance == 1 );
    
    tmpTable = LUTrialsTable(id, :);
    
    Trials = {};
    for Tr = 1:size(tmpTable, 1)
        trialstart = interp1(TimeStamps(:,2), TimeStamps(:,1), [tmpTable.trialstart(Tr)]);
        trialend   = interp1(TimeStamps(:,2), TimeStamps(:,1), [tmpTable.trialend(Tr)]);
        trialstart = round(trialstart); trialend = round(trialend);
        Trials{Tr} = downSampCSC(trialstart:trialend); % signal of interest in +/-250ms window around timestamp of interest
    end
    
    tmp = table;
    tmp.trlfp = Trials';
    [Type{1:size(tmp,1)}] = deal('Search');
    tmp.type = Type';
    TaskTable = tmp;
    clear Type tmp
    
    
    %%%%%%%%%%%% REST
    %%%%%%%%%%%% REST
    restIdxinfo = struct2table(dir(fullfile(restIdxDir, Session(1:13))));
    restIdxID = find(contains(restIdxinfo.name, 'blankScr'));
    restIdxinfo = restIdxinfo(restIdxID, :);
    
    Trials = {};
    for Rs = 1:size(restIdxinfo, 1)
        load(fullfile(restIdxinfo.folder{Rs}, restIdxinfo.name{Rs}))
        reststart = interp1(TimeStamps(:,2), TimeStamps(:,1), [restIdx.nlx_rest_start_ts]);
        restend   = interp1(TimeStamps(:,2), TimeStamps(:,1), [restIdx.nlx_rest_end_ts]);
        reststart = round(reststart); restend = round(restend);
        Trials{Rs}      = downSampCSC(reststart:restend); % signal of interest in +/-250ms window around timestamp of interest
    end
    
    tmp = table;
    tmp.trlfp = Trials';
    [Type{1:size(tmp,1)}] = deal('Rest');
    tmp.type = Type';
    TaskTable = [TaskTable; tmp];
    clear Type tmp
    
    %% Remove trials shorter than 1 seconds
    
    SIZE = cellfun(@(x) length(x), TaskTable.trlfp);
    IDX = find(SIZE < 2048);
    TaskTable(IDX, :) = [];
    
    %% Compute Power-Power Correlation
    Epochs = {'Search', 'Rest'};
    
    for ep = 1:length(Epochs)
        
        Index = find(contains(TaskTable.type,Epochs{ep}));
        
        ov = 0.5;   % overlap
        winlen = 1024; % winlen
        %     nff = max(256,2^nextpow2(winlen));
        nff = 150;
        freqBins = logspace(log10(1), log10(150), nff); %linspace(1, 100, nff);
        
        
        PSD = cell(length(Index), 1);
        PWELCH = NaN(length(Index), nff);
        
        for Tr = 1:length(Index)
            if ep == 1
                % if it's the search epoch, remove the first 500ms
                [pxx,f] = pwelch(TaskTable.trlfp{Index(Tr)}(500:end),winlen,floor(ov*winlen), freqBins, 1000);
            elseif ep == 2 & Tr == 1 & length(TaskTable.trlfp{Index(Tr)}) > 600000
                % Take 10min as rest
                [pxx,f] = pwelch(TaskTable.trlfp{Index(Tr)}(end-600000:end),winlen,floor(ov*winlen), freqBins, 1000);
            elseif ep == 2 & Tr == 2 & length(TaskTable.trlfp{Index(Tr)}) > 600000
                [pxx,f] = pwelch(TaskTable.trlfp{Index(Tr)}(1:600000),winlen,floor(ov*winlen), freqBins, 1000);
            else
                [pxx,f] = pwelch(TaskTable.trlfp{Index(Tr)},winlen,floor(ov*winlen), freqBins, 1000);
            end
            
            PWELCH(Tr, :) = pxx;
            
        end
        
        PwelchLU.session(Row) = {Session};
        PwelchLU.channel(Row) = channel;
        if ep == 1; PwelchLU.Search(Row) = {PWELCH}; end
        if ep == 2; PwelchLU.Rest(Row) = {PWELCH}; Row = Row + 1; end
        
    end
end
save('D:\HiTGun\Data\LU_PwelchLU_All_Sessions', 'PwelchLU')

%%

Pwelch_corrected = Pwelch;

Search = cellfun(@(x) mean(pow2db(x)), Pwelch_corrected.Search, 'UniformOutput', false);
Pwelch_corrected.Search = Search;
Rest = cellfun(@(x) mean(pow2db(x)), Pwelch_corrected.Rest, 'UniformOutput', false);
Pwelch_corrected.Rest = Rest;
Diff = cellfun(@(x, y) x-y, Pwelch_corrected.Rest, Pwelch_corrected.Search, 'UniformOutput', false);
Pwelch_corrected.Diff = Diff;


Search = cat(1, Pwelch_corrected.Search{:});
Rest = cat(1, Pwelch_corrected.Rest{:});
Diff = cat(1, Pwelch_corrected.Diff{:});

%% Plot
fh = figure();
% fh.WindowState = 'maximized';
ax_main = axes('Position', [0.1 0.2 0.8 0.55]);
ax_top = axes('Position', [0.1 0.8 0.8 0.1]);

%%%%%%%%%%%%
axes(ax_main);

Col = [150 29 78]/255;
h1 = plot_distribution(f,Search, 'Color', Col, 'LineWidth', 3)
hold on

Col = [224 164 88]/255;
h2 = plot_distribution(f,Rest, 'Color', Col, 'LineWidth', 3)

set(gca, 'box', 'off', 'tickdir', 'out', 'FontSize', 15, 'LineWidth', 1.5, 'xscale', 'log')
set(gca,'xtick',ceil(logspace(log10(f(1)),log10(f(end)),10)))

pbaspect([1 1 1])

xlabel('Log Frequency'); ylabel('10log10[power]')
legend([h1 h2], {'During Search', 'During Rest'})
legend boxoff

%%%%%%%%%%%%
axes(ax_top);
plot_distribution(f,Diff, 'Color', [0 0 0], 'LineWidth', 3)
ylabel('Rest-Search')
set(gca, 'box', 'off', 'tickdir', 'out', 'FontSize', 15,'LineWidth', 1.5, 'xscale', 'log')
ax_top.XTickLabel = [];
pbaspect([1 0.184 1])

sgtitle('Average Power Spectra Across All Trials')





