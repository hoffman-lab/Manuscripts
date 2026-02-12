function results = perpl_selectSVMHyperparams(X, y, varargin)
% -------------------------------------------------------------------------
% PURPOSE
%   Grid-search hyperparameter tuning for a *linear* SVM classifier
%   by selecting the BoxConstraint value that maximizes K-fold cross-validation
%   accuracy.
%
% INPUTS
%   X : numeric matrix, size (nObservations x nFeatures)
%       Feature matrix. Each row is an observation/sample; each column is a feature.
%
%   y : vector/categorical/cellstr/string array, length nObservations
%       Class labels corresponding to rows of X.
%
% NAME-VALUE OPTIONAL INPUTS
%   'BoxConstraints' : numeric vector (default: logspace(-3,3,10))
%       Candidate BoxConstraint values to evaluate.
%
%   'KFold' : positive integer (default: 10)
%       Number of cross-validation folds.
%
%   'Standardize' : logical scalar (default: true)
%       If true, standardize predictors internally in fitcsvm.
%
%   'FeatureMask' : [] or logical vector or numeric indices (default: [])
%       Features to treat as unstable.
%       - If logical vector (1 x nFeatures), true entries are masked.
%       - If numeric indices, those columns are masked.
%
%   'MaskMode' : 'zero' or 'remove' (default: 'zero')
%       What to do with masked features:
%       - 'zero'   : set masked feature columns to 0
%       - 'remove' : remove masked feature columns entirely
%
%   'Verbose' : logical scalar (default: true)
%       If true, prints best accuracy and BoxConstraint.
%
% OUTPUTS
%   results : struct with fields
%       .bestBoxConstraint   : scalar, best BoxConstraint found
%       .bestAccuracy        : scalar, best cross-validated accuracy
%       .boxConstraints      : vector, grid searched
%       .accuracies          : vector, accuracy per grid value (same order)
%       .cvLosses            : vector, loss per grid value (same order)
%       .bestModel           : SVM model trained on full (processed) X with bestBoxConstraint
%       .processedX          : feature matrix actually used (after masking/removal)
%       .maskInfo            : struct describing applied masking
%
% NOTES
%   - Accuracy is computed as: 1 - kfoldLoss(crossval(model,'KFold',K)).
%   - This function does not do nested CV; the selected hyperparameter is chosen
%     using CV on the full dataset and then a final model is fit on the full data.
% -------------------------------------------------------------------------

% ----------------------------- Parse inputs ------------------------------
p = inputParser;
p.FunctionName = mfilename;

addRequired(p, 'X', @(z) isnumeric(z) && ismatrix(z) && ~isempty(z));
addRequired(p, 'y', @(z) ~isempty(z) && numel(z) == size(X,1));

addParameter(p, 'BoxConstraints', logspace(-3, 3, 10), ...
    @(z) isnumeric(z) && isvector(z) && all(z > 0));

addParameter(p, 'KFold', 10, ...
    @(z) isnumeric(z) && isscalar(z) && z >= 2 && mod(z,1) == 0);

addParameter(p, 'Standardize', true, ...
    @(z) islogical(z) && isscalar(z));

addParameter(p, 'FeatureMask', [], ...
    @(z) isempty(z) || islogical(z) || (isnumeric(z) && isvector(z)));

addParameter(p, 'MaskMode', 'zero', ...
    @(s) ischar(s) || isstring(s));

addParameter(p, 'Verbose', true, ...
    @(z) islogical(z) && isscalar(z));

parse(p, X, y, varargin{:});

boxConstraints = p.Results.BoxConstraints(:)';
KFold          = p.Results.KFold;
doStandardize  = p.Results.Standardize;
featureMask    = p.Results.FeatureMask;
maskMode       = lower(string(p.Results.MaskMode));
verbose        = p.Results.Verbose;

% ----------------------------- Preprocess X ------------------------------
Xproc = X;

maskInfo = struct('enabled', false, 'mode', char(maskMode), ...
                  'maskType', '', 'maskedCols', []);

if ~isempty(featureMask)
    maskInfo.enabled = true;

    if islogical(featureMask)
        if numel(featureMask) ~= size(Xproc,2)
            error('FeatureMask (logical) must have length nFeatures = size(X,2).');
        end
        maskedCols = find(featureMask);
        maskInfo.maskType = 'logical';
    else
        maskedCols = unique(featureMask(:))';
        if any(maskedCols < 1) || any(maskedCols > size(Xproc,2))
            error('FeatureMask indices must be within [1, size(X,2)].');
        end
        maskInfo.maskType = 'indices';
    end

    maskInfo.maskedCols = maskedCols;

    switch maskMode
        case "zero"
            Xproc(:, maskedCols) = 0;

        case "remove"
            keepCols = true(1, size(Xproc,2));
            keepCols(maskedCols) = false;
            Xproc = Xproc(:, keepCols);

        otherwise
            error("MaskMode must be 'zero' or 'remove'.");
    end
end

% ----------------------------- Grid search -------------------------------
bestAccuracy = -Inf;
bestBox      = NaN;

accuracies = nan(1, numel(boxConstraints));
cvLosses   = nan(1, numel(boxConstraints));

for i = 1:numel(boxConstraints)
    box = boxConstraints(i);

    % Train linear SVM
    Mdl = fitcsvm( ...
        Xproc, y, ...
        'KernelFunction', 'linear', ...
        'BoxConstraint', box, ...
        'Standardize', doStandardize);

    % K-fold CV
    CVMdl = crossval(Mdl, 'KFold', KFold);
    loss  = kfoldLoss(CVMdl);
    acc   = 1 - loss;

    cvLosses(i)   = loss;
    accuracies(i) = acc;

    if acc > bestAccuracy
        bestAccuracy = acc;
        bestBox      = box;
    end
end

% ---------------------- Fit final model on full data ---------------------
bestModel = fitcsvm( ...
    Xproc, y, ...
    'KernelFunction', 'linear', ...
    'BoxConstraint', bestBox, ...
    'Standardize', doStandardize);

% ----------------------------- Package output ----------------------------
results = struct();
results.bestBoxConstraint = bestBox;
results.bestAccuracy      = bestAccuracy;
results.boxConstraints    = boxConstraints;
results.accuracies        = accuracies;
results.cvLosses          = cvLosses;
results.bestModel         = bestModel;
results.processedX        = Xproc;
results.maskInfo          = maskInfo;

% ----------------------------- Verbose print -----------------------------
if verbose
    fprintf('Best CV Accuracy: %.4f\n', bestAccuracy);
    fprintf('Best BoxConstraint: %.6g\n', bestBox);
end
end
