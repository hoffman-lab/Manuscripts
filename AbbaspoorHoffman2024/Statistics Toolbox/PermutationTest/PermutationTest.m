function y1 = PermutationTest(tf1,n1, tf2,n2, time, freq, n_permutes, pval1, pval2)

%% statistics via permutation testing
voxel_pval = pval1;
mcc_cluster_pval = pval2;
n_permutes = n_permutes;

permmaps = zeros(n_permutes,length(time),length(freq));


difference = mean(tf1, 3) - mean(tf2, 3);
difference = (difference - mean(difference))./std(difference)


tnum   = squeeze(mean(tf1, 3) - mean(tf2, 3));
tdenom = sqrt( (std(tf1,0,3).^2)./n1 +  (std(tf2,0,3).^2)./n2 );
real_t = tnum./tdenom;

% tf power maps are concatenated
tf3d = cat(3, tf1, tf2);

% generate maps under the null hypothesis
for permi = 1:n_permutes
    
    % randomize trials, which also randomly assigns trials to channels
    randorder = randperm(size(tf3d,3));
    temp_tf3d = tf3d(:,:,randorder);
    
    % compute the "difference" map
    tnum   = squeeze( mean(temp_tf3d(:,:,1:size(tf1, 3)),3) - mean(temp_tf3d(:,:,size(tf1, 3)+1:end),3) );
    tdenom = sqrt((std(temp_tf3d(:,:,1:size(tf1, 3)),0,3).^2)./n1 + (std(temp_tf3d(:,:,size(tf1, 3)+1:end),0,3).^2)./n2);
    tmap   = tnum./tdenom;
    
    permmaps(permi,:,:) = tmap;
    
    tmap(abs(tmap)<tinv(1-voxel_pval,n1-1))=0;
        
    clustinfo = bwconncomp(tmap);
    
    
    if numel(clustinfo.PixelIdxList)>0
        
        tempclustsizes = zeros(1, size(clustinfo.PixelIdxList, 2));
        
        for ii = 1:size(clustinfo.PixelIdxList, 2)
            
            tempclustsizes(ii) = sum(tmap(clustinfo.PixelIdxList{ii}));
            
        end
        
        max_clust_info(permi) = max(tempclustsizes);
        
    end
    
end

%% show non-corrected thresholded maps
% compute mean and standard deviation maps
zmap = (real_t-squeeze(mean(permmaps,1)))./squeeze(std(permmaps));

% apply cluster-level corrected threshold
zmapthresh = zmap;
% uncorrected pixel-level threshold
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=0;
% find islands and remove those smaller than cluster size threshold
clustinfo = bwconncomp(zmapthresh);
clust_info = cellfun(@numel,clustinfo.PixelIdxList);
clust_threshold = prctile(max_clust_info,100-mcc_cluster_pval*100);

% identify clusters to remove
whichclusters2remove = find(clust_info<clust_threshold);

% remove clusters
for i=1:length(whichclusters2remove)
    zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
end


%% plotting
figure('Renderer', 'painters', 'Position', [10 40 2400 500])

% hold on

h = subplot(1,3, 1)
p = get(h, 'position');
p(1) = p(1) - 0.1;
p(3) = p(3) + 0.07;
p(4) = p(4) - 0.03;
set(h, 'position', p);

imagesc(time - 1, freq, difference');
axis 'xy'
title('Difference maps')
line([0 0],[ylim],'Color','k','LineWidth',2,'LineStyle','-');
line([xlim],[10 10],'Color','k','LineWidth',2,'LineStyle','-.');
line([xlim],[20 20],'Color','k','LineWidth',2,'LineStyle','-.');

colorbar
cmap = twowaycol(256, 'RdYlBu10', 'pchip')
colormap(gca, cmap);

h = subplot(1,3, 2)
p = get(h, 'position');
p(1) = p(1) - 0.07;
p(3) = p(3) + 0.07;
p(4) = p(4) - 0.03;
set(h, 'position', p);

imagesc(time - 1, freq, zmapthresh');
axis 'xy'
title('Cluster-based-corrected thresholded maps')
line([0 0],[ylim],'Color','k','LineWidth',2,'LineStyle','-');
line([xlim],[10 10],'Color','k','LineWidth',2,'LineStyle','-.');
line([xlim],[20 20],'Color','k','LineWidth',2,'LineStyle','-.');

colorbar
cmap = twowaycol(256, 'RdYlBu10', 'pchip')
colormap(gca, cmap);

end