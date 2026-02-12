function [stateSeq_sorted, stateSeq_raw, estTR, estE, binEdges, TR0, E0] = ...
         perpl_MultistateHMM(x, nStates, nBins)
% HMM_ACTIVATION_MULTISTATE_KMEANS
% Use k-means to initialize HMM transition (TR0) and emission (E0) matrices,
% then train a discrete-emission HMM (MATLAB 2021) on activation signal x.
%
% x       : activation vector (T x 1 or 1 x T)
% nStates : number of HMM states (default: 3)
% nBins   : number of bins for discretization (default: 20)
%
% OUTPUT:
%   stateSeq_sorted : T x 1, states 1..nStates ordered low→high activation
%   stateSeq_raw    : T x 1, raw Viterbi state labels (before reordering)
%   estTR           : nStates x nStates, learned transition matrix
%   estE            : nStates x nBins, learned emission matrix
%   binEdges        : 1 x (nBins+1), edges used for discretization
%   TR0, E0         : initial transition/emission matrices from k-means

    % ----------------------------
    % Defaults
    % ----------------------------
    if nargin < 2 || isempty(nStates)
        nStates = 3;
    end
    if nargin < 3 || isempty(nBins)
        nBins = 20;
    end

    % ----------------------------
    % Ensure column vector & valid
    % ----------------------------
    x = x(:);
    T = numel(x);

    if min(x) == max(x)
        error('Activation vector is constant — cannot cluster.');
    end

    % ----------------------------
    % Discretize x into bins (for HMM emissions)
    % ----------------------------
    binEdges = linspace(min(x), max(x), nBins+1);
    seq = discretize(x, binEdges);       % 1..nBins
    seq(isnan(seq)) = nBins;             % handle edge case
    seq = seq(:)';                       % row vector for hmmtrain

    % ----------------------------
    % K-means on continuous x
    % ----------------------------
    % Initial clustering on activation values
    [idx, C] = kmeans(x, nStates);
    
%     , ...
%                       'Replicates', 10, ...
%                       'MaxIter',    200, ...
%                       'Display',   'off'

    % Sort clusters by mean activation (low→high)
    [~, orderC] = sort(C, 'ascend');    % orderC(k) = old label of kth lowest cluster

    newIdx = zeros(size(idx));
    for k = 1:nStates
        newIdx(idx == orderC(k)) = k;   % relabel clusters to 1..nStates ordered
    end
    idx = newIdx;                       % idx(t) ∈ {1..nStates}, ordered

    % ----------------------------
    % Initialize TR0 from k-means sequence
    % ----------------------------
    TR0 = zeros(nStates);
    for t = 1:(T-1)
        i = idx(t);
        j = idx(t+1);
        TR0(i,j) = TR0(i,j) + 1;
    end

    % Add pseudocounts to avoid zeros
    TR0 = TR0 + 1;
    TR0 = TR0 ./ sum(TR0, 2);          % row-normalize

    % ----------------------------
    % Initialize E0 from k-means & discretized seq
    % ----------------------------
    E0 = zeros(nStates, nBins);
    edges_sym = 0.5 : 1 : (nBins + 0.5);  % integer-centered bin edges for histcounts

    for k = 1:nStates
        mask = (idx == k);
        if any(mask)
            % Symbol counts for this state over discrete seq
            counts = histcounts(seq(mask), edges_sym);
            % Add pseudocounts for stability
            counts = counts + 1;
            E0(k,:) = counts / sum(counts);
        else
            % If somehow empty, use uniform
            E0(k,:) = ones(1, nBins) / nBins;
        end
    end

    % ----------------------------
    % Train HMM using built-in hmmtrain
    % ----------------------------
    [estTR, estE] = hmmtrain(seq, TR0, E0, ...
                             'Maxiterations', 2000, ...
                             'Tolerance',     1e-6, ...
                             'Verbose',       false);

    % ----------------------------
    % Viterbi decode
    % ----------------------------
    stateSeq_raw = hmmviterbi(seq, estTR, estE);
    stateSeq_raw = stateSeq_raw(:);   % T x 1

    % ----------------------------
    % Sort final states by mean activation
    % ----------------------------
    meanActivation = nan(nStates,1);
    for k = 1:nStates
        if any(stateSeq_raw == k)
            meanActivation(k) = mean(x(stateSeq_raw == k));
        end
    end

    [~, order] = sort(meanActivation, 'ascend');  % low → high

    stateSeq_sorted = zeros(T,1);
    for k = 1:nStates
        stateSeq_sorted(stateSeq_raw == order(k)) = k;
    end
end
