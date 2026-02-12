function [C, P, nPairs] = cross_transition_matrix(refStates, tgtStates, nStates, varargin)
% CROSS_TRANSITION_MATRIX
%   Flexible state–state relationship between ref and tgt assemblies.
%
%   [C, P, nPairs] = cross_transition_matrix(refStates, tgtStates, nStates)
%       Default: lag = 1 (ref(t) -> tgt(t+1)),
%       row-normalized conditional probabilities P(j|i).
%
%   [C, P, nPairs] = cross_transition_matrix(..., 'Lag', tau)
%       Uses ref(t) -> tgt(t+tau). tau must be an integer (can be 0).
%
%   Inputs:
%       refStates : vector of reference HMM states (1..nStates)
%       tgtStates : vector of target   HMM states (1..nStates)
%       nStates   : scalar, number of possible states
%
%   Outputs:
%       C       : nStates×nStates count matrix
%       P       : nStates×nStates probability matrix
%       nPairs  : total number of valid (i,j) pairs used

    % ---------- Parse optional inputs ----------
    p = inputParser;
    addParameter(p, 'Lag', 1, @(x) isnumeric(x) && isscalar(x));
    parse(p, varargin{:});

    lag  = round(p.Results.Lag);       % enforce integer

    if lag < 0
        error('Lag must be >= 0 for this implementation.');
    end

    refStates = refStates(:);
    tgtStates = tgtStates(:);
    
    maxLen = min(numel(refStates), numel(tgtStates) - lag);
    idxRef = (1:maxLen).';
    idxTgt = idxRef + lag;

    sRef = refStates(idxRef);
    sTgt = tgtStates(idxTgt);

    % ---------- Count matrix ----------
    C = zeros(nStates);
    for k = 1:numel(sRef)
        i = sRef(k);
        j = sTgt(k);
        if i>=1 && i<=nStates && j>=1 && j<=nStates
            C(i,j) = C(i,j) + 1;
        end
    end

    nPairs = sum(C,'all');

    % ---------- Probability matrix ----------
    P = zeros(size(C));

    % Row-normalized P(j | i)
    rowSums = sum(C,2);
    for i = 1:nStates
        if rowSums(i) > 0
            P(i,:) = C(i,:) ./ rowSums(i);
        end
    end
end
