%% This is modified from CellExplorer script
%%% Saman Abbaspoor - 03/12/2024 - Vanderbilt


function deepSuperficialfromRipple = perpl_classify_DeepSuperficial(basepath,ripple,varargin)
% SWR ripples classification of Deep-superficial channelss
% Defines the deep superficial boundary defined by the reversal of the average sharp wave.
%
% The algorith assigned both distance to the reversal point and labels.
% The assigned label are (ordered by depth):
%   'Cortical'    : Assigned to channels belonging to a spikegroup with the channelTag Cortical
%   'Deep'        : Assigned to channels above the reversal
%   'Superficial' : Assigned to channels below the reversal
%   ''            : Assigned to channels belonging to a spikegroup with the channelTag Bad
%
% Distance to the boundary in µm for electrode groups where the reversal is detected.
%
% INPUT
% ripple: ripple table where each row contains the info on a single ripple event
% with a column RipLFP which contains the ripple-aligned LFP (+-500ms) across all
% channels

% Part of CellExplorer: https://CellExplorer.org
% By Peter Petersen
% petersen.peter@gmail.com
% Last edited: 20-02-2021

p = inputParser;
addParameter(p,'ripple_channels',[1:64],@isvector);
addParameter(p,'channels_to_exclude',[],@isvector);
addParameter(p,'conv_length',2,@isnumerical);
addParameter(p,'VerticalSpacing',90,@isnumerical);
addParameter(p,'saveMat',true,@islogical); % Defines if a mat file is created
addParameter(p,'saveFig',false,@islogical); % Defines if the summary figure is saved to basepath


% Parsing inputs
parse(p,varargin{:})
ripple_channels{1} = p.Results.ripple_channels;
channels_to_exclude{1} = p.Results.channels_to_exclude;
conv_length = p.Results.conv_length;
VerticalSpacing = p.Results.VerticalSpacing;
saveMat = p.Results.saveMat;
saveFig = p.Results.saveFig;


ripple_average = [];
ripple_power = [];
ripple_amplitude = [];
ripple_time_axis = [-150:150];
nChannels = numel(ripple_channels{1});
deepSuperficial_ChClass3 = repmat({''},1,nChannels);
deepSuperficial_ChClass = repmat({''},1,nChannels);
deepSuperficial_ChDistance3 = nan(1,nChannels);
deepSuperficial_ChDistance = nan(1,nChannels);

jj = 1;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Main algorithm  (two versions below)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
ConfirmedRipples = ripple.ConfirmedRipples;
idx = find(ConfirmedRipples == 1);
LFP = ripple.RipLFP(idx, :);
lfp = zeros(numel(idx), size(ripple.RipLFP{1,1}, 1),  size(ripple.RipLFP{1,1}, 2));
for event = 1:numel(idx)
    lfp(event, :, :) = LFP{event, 1};
end
if size(lfp, 1) == 1
    ripple_average{jj} = squeeze(lfp(:, ripple_channels{1}, 100:400))';
else
    ripple_average{jj} = squeeze(mean(lfp(:, ripple_channels{1}, 100:400)))';
end
ripple_power{jj} = sum(ripple_average{jj}.^2)./max(sum(ripple_average{jj}.^2));
ripple_amplitude{jj} = mean(ripple_average{jj})/max(abs(mean(ripple_average{jj})));

SWR_slope = [];

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Second algorithm (newest)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Seems to do much better across mouse and rat recordings where the
% sharp wave can have very different shape and amplitude.

[~,ia,~] = intersect(ripple_channels{jj}, channels_to_exclude{1});
[ripple_channels2,ia2] = setdiff(ripple_channels{jj}, channels_to_exclude{1},'legacy');
ia2 = sort(ia2);
SWR_average = nanconv(ripple_average{jj},ones(40,1)/40,'edge');
SWR_average = SWR_average-SWR_average(100,:);
SWR_diff{jj} = sum(SWR_average(100:150,:)); %135
SWR_diff2 = SWR_diff{jj};
SWR_diff2(ia) = [];
SWR_diff2 = nanconv(SWR_diff2,[ones(1,conv_length)]/conv_length,'edge');

SWR_average2 = nanconv(ripple_average{jj},ones(20,1)/20,'edge');
SWR_amplitude{jj} = sum(abs(ripple_average{jj}(100:201,:)-SWR_average2(100:201,:)));
SWR_amplitude2 = SWR_amplitude{jj};
SWR_amplitude2(ia) = [];

coefs = polyfit([1:length(SWR_amplitude2)], SWR_amplitude2, 1);
SWR_slope = coefs(1);

SWR_diff{jj} = SWR_diff{jj}./max(abs(SWR_diff{jj}));
SWR_amplitude{jj} = (SWR_amplitude{jj}-min(SWR_amplitude{jj}))./max(abs(SWR_amplitude{jj}-min(SWR_amplitude{jj})));

if any(diff(SWR_diff2<0)==1) %&& ~any(jj == electrodeGroupsToExclude)
    indx = find(diff(SWR_diff2<0)==1);indx = indx(end);
    threshold = interp1(SWR_diff2([indx,indx+1]),[ia2(indx),ia2(indx+1)],0);
    deepSuperficial_ChClass(ripple_channels{jj}([1:threshold])) = repmat({'Deep'},length([1:threshold]),1);
    deepSuperficial_ChClass(ripple_channels{jj}([ceil(threshold):size(SWR_diff{jj},2)])) = repmat({'Superficial'},length([ceil(threshold):size(SWR_diff{jj},2)]),1);
    deepSuperficial_ChDistance(ripple_channels{jj}) = (1:size(SWR_diff{jj},2))-threshold;
elseif any(diff(SWR_diff2<0)==-1) %&& ~any(jj == electrodeGroupsToExclude)
    indx = find(diff(SWR_diff2<0)==-1);indx = indx(end);
    threshold = interp1(SWR_diff2([indx,indx+1]),[ia2(indx),ia2(indx+1)],0);
    deepSuperficial_ChClass(ripple_channels{jj}([1:threshold])) = repmat({'Superficial'},length([1:threshold]),1);
    deepSuperficial_ChClass(ripple_channels{jj}([ceil(threshold):size(SWR_diff{jj},2)])) = repmat({'Deep'},length([ceil(threshold):size(SWR_diff{jj},2)]),1);
    deepSuperficial_ChDistance(ripple_channels{jj}) = (1:size(SWR_diff{jj},2))-threshold;
else
    if SWR_slope > 0
        deepSuperficial_ChClass(ripple_channels{jj}) = repmat({'Deep'},length(ripple_channels{jj}),1); % Deep
        if SWR_diff2(end)*5<max(SWR_diff2) && SWR_diff2(end) > 0
            deepSuperficial_ChDistance(ripple_channels{jj}) = (1:size(SWR_diff{jj},2))-length(ripple_channels{jj})-1;
        end
    else
        deepSuperficial_ChClass(ripple_channels{jj}) = repmat({'Superficial'},length(ripple_channels{jj}),1); % Superficial
        if SWR_diff2(1)*5>min(SWR_diff2)  && SWR_diff2(1) < 0
            deepSuperficial_ChDistance(ripple_channels{jj}) = (1:size(SWR_diff{jj},2))+1;
        end
    end
end


clear signal


deepSuperficial_ChDistance3 = deepSuperficial_ChDistance3 * VerticalSpacing;
deepSuperficial_ChDistance = deepSuperficial_ChDistance * VerticalSpacing;

% Saving the result to basename.deepSuperficialfromRipple.channelinfo.mat
deepSuperficialfromRipple.channel = [1:length(deepSuperficial_ChDistance)]'; % 1-indexed
deepSuperficialfromRipple.channelClass = deepSuperficial_ChClass';
deepSuperficialfromRipple.channelDistance = deepSuperficial_ChDistance';
deepSuperficialfromRipple.ripple_power = ripple_power;
deepSuperficialfromRipple.ripple_amplitude = ripple_amplitude;
deepSuperficialfromRipple.ripple_average = ripple_average;
deepSuperficialfromRipple.ripple_time_axis = ripple_time_axis;
deepSuperficialfromRipple.ripple_channels = ripple_channels; %  index 1 for channels
deepSuperficialfromRipple.SWR_diff = SWR_diff;
deepSuperficialfromRipple.SWR_amplitude = SWR_amplitude;
% if isfield(ripples,'detectorinfo')
%     deepSuperficialfromRipple.detectorinfo = ripples.detectorinfo;
% end
deepSuperficialfromRipple.processinginfo.function = 'classification_deepSuperficial';
deepSuperficialfromRipple.processinginfo.date = now;
deepSuperficialfromRipple.processinginfo.params.verticalSpacing = VerticalSpacing;
% deepSuperficialfromRipple.processinginfo.params.electrodeGroupsToExclude = electrodeGroupsToExclude;
if saveMat
    save(fullfile(basepath, 'deepSuperficialfromRipple.channelinfo.mat'),'deepSuperficialfromRipple');
end

% Plotting the average ripple with sharp wave across all electrode groups


figure
subplot(1, 2, 1)
plot(SWR_diff{jj},-[1:size(SWR_diff{jj},2)],'-k','linewidth',2), hold on, grid on
plot(SWR_amplitude{jj},-[1:size(SWR_amplitude{jj},2)],'m','linewidth',1)

for jjj = 1:size(SWR_diff{jj},2)
    % Plotting depth (µm)
    if strcmp(deepSuperficial_ChClass(ripple_channels{jj}(jjj)),'Superficial')
        plot(SWR_diff{jj}(jjj),-(jjj),'or','linewidth',2)
    elseif strcmp(deepSuperficial_ChClass(ripple_channels{jj}(jjj)),'Deep')
        plot(SWR_diff{jj}(jjj),-(jjj),'ob','linewidth',2)
    elseif strcmp(deepSuperficial_ChClass(ripple_channels{jj}(jjj)),'Cortical')
        plot(SWR_diff{jj}(jjj),-(jjj),'og','linewidth',2)
    else
        plot(SWR_diff{jj}(jjj),-(jjj),'ok')
    end
end

xlabel('SWR diff')
ylabel('Channels')
plot([0 0], [1 -numel(ripple_channels{1})], '--', 'Color', [0 0 0])
ylim([-numel(ripple_channels{1}) -1])
set(gca, 'ytick', -numel(ripple_channels{1}):1:1, 'yticklabel', numel(ripple_channels{1}):-1:1, 'box', 'off', 'tickdir', 'out', 'FontSize', 15)

ripple_average = deepSuperficialfromRipple.ripple_average;
jjj = 1:size(SWR_diff{jj},2);
ripple_average = cellfun(@(X) zscore(X, [], 1)-(jjj-1)*3, ripple_average, 'UniformOutput', false);
Max = max(ripple_average{1}, [], 'all'); MIN = min(ripple_average{1}, [], 'all');

subplot(1, 2, 2)

for jjj = 1:size(SWR_diff{jj},2)
    % Plotting depth (µm)
    if strcmp(deepSuperficial_ChClass(ripple_channels{jj}(jjj)),'Superficial')
        color = [1 0 0];
    elseif strcmp(deepSuperficial_ChClass(ripple_channels{jj}(jjj)),'Deep')
        color = [0 0 1];
    elseif strcmp(deepSuperficial_ChClass(ripple_channels{jj}(jjj)),'Cortical')
        color = [0 1 0];
    else
        color = [0 0 0];
    end
    
    plot(ripple_time_axis,ripple_average{jj}(:,jjj), 'Color', color)
    hold on
end


ylim([MIN Max])

jjj = 1:size(SWR_diff{jj},2);
Points = ripple_average{1}(1,:);
set(gca, 'ytick', flip(Points), 'yticklabel', flip(jjj), 'box', 'off', 'tickdir', 'out', 'FontSize', 15)


% Saving figure
if saveFig
    saveas(gcf,'deepSuperficial_classification_fromRipples.png');
end
