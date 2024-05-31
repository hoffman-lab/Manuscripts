function [permmaps, mean_h0, std_h0, zval, cluster_thresh] = ...
    Bicoherence_Montecarlo_ClusterBased(n_permutes, Freq, pval, signal_length)

% convert p-value to Z value
zval = abs(norminv(pval));

permmaps = bicoher_montecarlo2(signal_length, Freq, n_permutes, pval);

% compute mean and standard deviation maps
mean_h0 = squeeze(mean(permmaps, 3));
std_h0  = squeeze(std(permmaps, [], 3));


% initialize matrices for cluster-based correction
max_cluster_sizes = zeros(1,n_permutes);
% ... and for maximum-pixel based correction
max_val = zeros(n_permutes,2); % "2" for min/max

% loop through permutations
for permi = 1:n_permutes
    
    % take each permutation map, and transform to Z
    threshimg = squeeze(permmaps(:,:, permi));
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
    temp = sort( reshape(permmaps(:,:, permi), 1,[] ));
    max_val(permi,:) = [ min(temp) max(temp) ];
end

cluster_thresh = prctile(max_cluster_sizes,100-(100*pval));

end
