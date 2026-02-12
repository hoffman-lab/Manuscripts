%##################################################################
% Weighted clustering coefficient.
% Source: Barrat et al, The architecture of complex weighted networks
%
% INPUTS: weighted adjacency matrix, nxn
% OUTPUTs: vector of node weighted clustering coefficients, nx1
%
% Other routines used: degrees.m, kneighbors.m
% GB: last updated, Sep 30 2012
%##################################################################

function wC=weightedClustCoeff(adj)

[deg,~,~]=degrees(adj);
n=size(adj,1); % number of nodes
wC=zeros(n,1); % initialize weighted clust coeff

for i=1:n % across all nodes
    neigh=kneighbors(adj,i,1);
    if length(neigh)<2; continue; end
    
    s=0;
    for ii=1:length(neigh)
        for jj=1:length(neigh)
          
            if adj(neigh(ii),neigh(jj))>0; s=s+(adj(i,neigh(ii))+adj(i,neigh(jj)))/2; end
        
        end
    end
   
    wC(i)=s/(deg(i)*(length(neigh)-1));
end

end

%##################################################################
% Compute the total degree, in-degree and out-degree of a graph based on the adjacency matrix; 
% Note: Returns weighted degrees, if the input matrix is weighted
%
% INPUTS: adjacency matrix, nxn
% OUTPUTS: degree (1xn), in-degree (1xn) and out-degree (1xn) sequences
%
% Other routines used: isDirected.m
% GB: last updated, Sep 26, 2012
%##################################################################

function [deg,indeg,outdeg]=degrees(adj)

indeg = sum(adj);
outdeg = sum(adj');

if isDirected(adj)
  deg = indeg + outdeg; % total degree

else   % undirected graph: indeg=outdeg
  deg = indeg + diag(adj)';  % add self-loops twice, if any

end

end


function neighbors = kneighbors(A,ind,k)
%KNEIGHBORS Finds the adjacency list of nodes at distance k from 'ind'
%
% @input A, NxN adjacency matrix
% @input ind, a scalar node # of the search index
% @input k, a scalar for distance to search
% @output neighbors, a adjacency list of nodes reachable from 'ind' in 'k' hops 
% INPUTS: adjacency matrix (nxn), start node index, k - number of links
% OUTPUTS: vector of k-neighbors indices
%
% Updated: For readability.

% IB: last updated, 3/23/14

adjk = A;
for i=1:k-1 
    adjk = adjk*A; 
end;

neighbors = find(adjk(ind,:)>0);

end



%##################################################################
% Checks whether the graph is directed, using the matrix transpose function
%
% INPUTS: adjacency matrix, nxn
% OUTPUTS: boolean variable, 0 or 1
%
% Note: one-liner alternative: S=not(issymmetric(adj));
% GB: last updated, Sep 23, 2012
%##################################################################

function S=isDirected(adj)

S = true;
if adj==transpose(adj); S = false; end

end





% % ALTERNATIVE  =========================================================
% wadj=adj;
% adj=adj>0;
% 
% [wdeg,~,~]=degrees(wadj);
% [deg,~,~]=degrees(adj);
% n=size(adj,1); % number of nodes
% wC=zeros(n,1);
% 
% for i=1:n
%     if deg(i)<2; continue; end
%     
%     s=0;
%     for ii=1:n
%         for jj=1:n
%             s=s+adj(i,ii)*adj(i,jj)*adj(ii,jj)*(wadj(i,ii)+wadj(i,jj))/2;
%         end
%     end
% 
%     wC(i)=s/(wdeg(i)*(deg(i)-1));
% end