%% Saman Abbaspoor - 03/12/2024 - Hoffman Lab

function [zmap] = ...
        FOOOF_Permutation_ClusterBased(Sample1, Sample2, n_permutes, pval)
   
    
diffmap = squeeze(mean(Sample1)) - squeeze(mean(Sample2));
    

num_frex = size(Sample1, 3);
num_chan = size(Sample1, 2);
ntrials  = size(Sample1, 1);
% convert p-value to Z value
zval = abs(norminv(pval));

% initialize null hypothesis maps
permmaps = zeros(n_permutes,num_chan,num_frex);

% for convenience, tf power maps are concatenated
%   in this matrix, trials 1:ntrials are from channel "1" 
%   and trials ntrials+1:end are from channel "2"
catfile = cat(1,Sample1,Sample2);


% generate maps under the null hypothesis
for permi = 1:n_permutes
    
    % randomize trials, which also randomly assigns trials to channels
    randorder = randperm(size(catfile,1));
    temp_catfile = catfile(randorder,:,:);
    
    % compute the "difference" map
    % what is the difference under the null hypothesis?
    permmaps(permi,:,:) = squeeze( mean(temp_catfile(1:ntrials, :,:),1) - mean(temp_catfile(ntrials+1:end, :,:),1) );
end


%% show non-corrected thresholded maps

% compute mean and standard deviation maps
mean_h0 = squeeze(mean(permmaps));
std_h0  = squeeze(std(permmaps));

% now threshold real data...
% first Z-score
zmap = (diffmap-mean_h0) ./ std_h0;

% threshold image at p-value, by setting subthreshold values to 0
zmap(abs(zmap)<zval) = 0;


%% corrections for multiple comparisons

% initialize matrices for cluster-based correction
max_cluster_sizes = zeros(1,n_permutes);
% ... and for maximum-pixel based correction
max_val = zeros(n_permutes,2); % "2" for min/max

% loop through permutations
for permi = 1:n_permutes
    
    % take each permutation map, and transform to Z
    threshimg = squeeze(permmaps(permi,:,:));
    threshimg = (threshimg-mean_h0)./std_h0;
    
    % threshold image at p-value
    threshimg(abs(threshimg)<zval) = 0;
    
    
    % find clusters (need image processing toolbox for this!)
    islands = bwconncomp(threshimg);
    if numel(islands.PixelIdxList)>0
        
        % count sizes of clusters
        tempclustsizes = cellfun(@length,islands.PixelIdxList);
        
        % store size of biggest cluster
        max_cluster_sizes(permi) = max(tempclustsizes);
    end
    
    
    % get extreme values (smallest and largest)
    temp = sort( reshape(permmaps(permi,:,:),1,[] ));
    max_val(permi,:) = [ min(temp) max(temp) ];
    
end

%% show histograph of maximum cluster sizes
% find cluster threshold (need image processing toolbox for this!)
% based on p-value and null hypothesis distribution
cluster_thresh = prctile(max_cluster_sizes,100-(100*pval));

%% plots with multiple comparisons corrections

% now find clusters in the real thresholded zmap
% if they are "too small" set them to zero
islands = bwconncomp(zmap);
for i=1:islands.NumObjects
    % if real clusters are too small, remove them by setting to zero!
    if numel(islands.PixelIdxList{i}==i)<cluster_thresh
        zmap(islands.PixelIdxList{i})=0;
    end
end  
    
end

