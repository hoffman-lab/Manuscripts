clear all, clc
FilesDir = 'D:\HiTGun\Data\decimatedCSC_wholeSession';
TTLFilesDir = 'D:\HiTGun\Data\LE Change blindness Raw';

load('D:\HiTGun\Data\LETrialsTable');
load(fullfile('D:\HiTGun\Data', 'LE_In_Layer_Channel.mat'));

info = struct2table(dir('D:\HiTGun\Data\LE Change blindness Raw'));
SessionTable = info.name;
SessionTable = SessionTable(3:end);

Pwelch = table;
Row = 1;
for Sess = 1:size(SessionTable, 1)
    
    %% Load Data
    Session = SessionTable{Sess};
    
    Index = find(contains(string(LE_In_Layer_Channel.sid),Session));
    inchannel = LE_In_Layer_Channel.rip_cscNum(Index);
    
    channellist = [inchannel 13 14 15];    
    channellist = unique(channellist);
    
    for ch = 1:length(channellist)
        
        channel = {['CSC', num2str(channellist(ch))]};
        
        %%% Load Channel LFP
        filedir = fullfile(FilesDir, Session, channel);
        load(filedir{1})
        
        
        filedir = fullfile(TTLFilesDir, Session);
        try
        [TimeStampsEV, ~, TTLs, ~, EventString Header] = ...
            Nlx2MatEV( fullfile(filedir, 'Events.nev'), [1 1 1 1 1], 1, ...
            1, 1 );
        catch
            continue
        end
                
        %% Find Event of Interest
        
        %%%%%%%%%%%%Search
        id = find( contains(LETrialsTable.session,Session) &...
            contains(LETrialsTable.type, 'Search')&...
            LETrialsTable.performance == 1 );
        
        tmpTable = LETrialsTable(id, :);
        
        Trials = {};
        for Tr = 1:size(tmpTable, 1)
            trialstart = interp1(TimeStamps(:,2), TimeStamps(:,1), [tmpTable.trialstart(Tr)]); trialstart = round(trialstart);
            trialend = interp1(TimeStamps(:,2), TimeStamps(:,1), [tmpTable.trialend(Tr)]); trialend = round(trialend);
            Trials{Tr}      = decimateCSC(trialstart:trialend); % signal of interest in +/-250ms window around timestamp of interest
        end
        
        tmp = table;
        tmp.trlfp = Trials';
        [Type{1:size(tmp,1)}] = deal('Search');
        tmp.type = Type';
        TaskTable = tmp;
        clear Type tmp
        
        %%%%%%%%%%%% REST
        %%%%%%%%%%%% REST
        TaskStarts = find(TTLs == 23); TaskStarts = TaskStarts(end);
        TaskEnds = find(TTLs == 24); TaskEnds = TaskEnds(end);
        
        Trials_ev = NaN(2, 1);
        Trials_ev(1) = TimeStampsEV(TaskStarts);
        Trials_ev(2) = TimeStampsEV(TaskEnds);
        
        
        Trials = {};
        for Tr = 1:length(Trials_ev)
            rtime = interp1(TimeStamps(:,2), TimeStamps(:,1), [Trials_ev(Tr,1)]); rtime = round(rtime);
            if Tr == 1
                Trials{Tr}      = decimateCSC(1:rtime);
            elseif Tr ==2
                Trials{Tr}      = decimateCSC(rtime:end);
            end
        end
        
        tmp = table;
        tmp.trlfp = Trials';
        [Type{1:size(tmp,1)}] = deal('Rest');
        tmp.type = Type';
        TaskTable = [TaskTable; tmp];
        clear Type tmp
        
        
        %% Remove trials shorter than 1 seconds
        
        SIZE = cellfun(@(x) length(x), TaskTable.trlfp);
        IDX = find(SIZE < 1024);
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
            
            Pwelch.session(Row) = {Session};
            Pwelch.channel(Row) = channel;
            if ep == 1; Pwelch.Search(Row) = {PWELCH}; end
            if ep == 2; Pwelch.Rest(Row) = {PWELCH}; Row = Row + 1; end
            
        end
    end
end

save('D:\HiTGun\Data\LE_Pwelch_All_Sessions_All_Channels', 'Pwelch')

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

% save('LELU_Pwelch_All_Sessions', 'Search', 'Rest')





%% Plot
fh = figure();
% fh.WindowState = 'maximized';
ax_main = axes('Position', [0.1 0.2 0.8 0.55]);
ax_top = axes('Position', [0.1 0.8 0.8 0.1]);
% incet = axes('Position', [0.33 0.25 0.2 0.2]);

%%%%%%%%%%%%
axes(ax_main);

Col = [150 29 78]/255;
h1 = plot_distribution(f,pow2db(Search), 'Color', Col, 'LineWidth', 3)
hold on

Col = [224 164 88]/255;
h2 = plot_distribution(f,pow2db(Rest), 'Color', Col, 'LineWidth', 3)

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


%%

figName = [Session]
Figure_Output_Directory = 'D:\HiTGun\Figures\LE'
savefig(fullfile(Figure_Output_Directory, figName))
saveas(gcf, fullfile(Figure_Output_Directory, [figName '.png']))
