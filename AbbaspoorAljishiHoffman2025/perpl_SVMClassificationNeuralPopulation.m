function LinearSVMConditionificationTable = perpl_SVMConditionificationNeuralPopulation(FR, Condition, FeatureMask, varargin)
% perpl_SVMConditionificationNeuralPopulation
% -------------------------------------------------------------------------
% PURPOSE
%   Train a linear SVM to Conditionify observations using a neural-population design
%   matrix, tune the BoxConstraint hyperparameter using Bayesian optimization,
%   evaluate cross-validated accuracy, and run multiple permutation-based tests
%   that probe different properties of the neural population code:
%
%     (1) Label permutation test (tests whether accuracy exceeds chance given labels)
%     (2) Activity heterogeneity permutation (shuffles neuron identities within trials)
%     (3) Within-subpopulation heterogeneity permutation (shuffles within +/− weight groups)
%     (4) Noise-correlation disruption permutation (shuffles trials within Condition per neuron)
%     (5) Cofiring structure analysis (+ permutation test):
%         compares within-group vs across-group correlations (groups defined by SVM weight sign)
%
% INPUTS
%   FR : numeric matrix, size (nTrials x nNeurons)
%       Design matrix used for Conditionification. Rows are trials/observations and
%       columns are neurons/features (e.g., firing rates in a chosen epoch).
%
%   Condition : vector/categorical/cellstr/string array, length nTrials
%       Condition labels for each trial (must match the number of rows in FR).
%
%   FeatureMask : logical or numeric mask for unstable neurons/features
%       If logical: length nNeurons, true entries are set to 0 in FR.
%       If numeric indices: those columns are set to 0 in FR.
%
% NAME-VALUE OPTIONAL INPUTS
%   'KFold' : positive integer (default: 10)
%       Number of cross-validation folds.
%
%   'BoxConstraintRange' : 1x2 numeric vector (default: [1e-10, 1])
%       Range used for Bayesian optimization of BoxConstraint (log-transformed).
%
%   'MaxObjectiveEvaluations' : positive integer (default: 100)
%       Maximum objective evaluations for Bayesian optimization.
%
%   'AcquisitionFunctionName' : char/string (default: 'expected-improvement-plus')
%       Acquisition function used by Bayesian optimization.
%
%   'NumPermutations' : positive integer (default: 5000)
%       Number of permutations for permutation tests (label/activity/subpop/noise/corr).
%
%   'NumBootstraps' : positive integer (default: 10000)
%       Number of bootstrap resamples for confidence intervals on permutation distributions.
%
%   'Alpha' : scalar in (0,1) (default: 0.05)
%       Significance level for bootstrap confidence intervals.
%
%   'UseParallel' : logical scalar (default: true)
%       If true, uses parfor for permutation loops when possible.
%
%   'ConditionValues' : [] or 1x2 cell array / string array / numeric array (default: [])
%       If provided, these two values define the two conditions for the
%       noise-correlation disruption permutation. If empty, the function
%       attempts to infer two unique Conditiones from Condition.
%
% OUTPUTS
%   LinearSVMConditionificationTable : struct
%       A struct with fields mirroring the original table variables:
%         .SVMMdl                    trained SVM model (best BoxConstraint)
%         .CVMdl                     cross-validated model object
%         .cvLoss                    cross-validation loss
%         .cvAccuracy                cross-validation accuracy
%         .ConditionPermAccuracy         permutation accuracies (label shuffle)
%         .ConditionPermAccuracyStats    [CI_low, mean, CI_high] for label shuffle
%         .ConditionPermPvalue           p-value for label shuffle
%         .weights                   SVM weights (Beta)
%         .normWeights               L2-normalized SVM weights
%         .bias                      SVM bias term
%         .ActPermAccuracy           permutation accuracies (neuron identity shuffle)
%         .ActPermAccuracyStats      CI stats for neuron identity shuffle
%         .ActPermPvalue             p-value for neuron identity shuffle
%         .MatchPermAccuracy         accuracies (shuffle within +/− weight subpops)
%         .MatchPermAccuracyStats    CI stats for within-subpop shuffle
%         .MatchPermPvalue           p-value for within-subpop shuffle
%         .NoisePermAccuracy         accuracies (destroy noise correlations)
%         .NoisePermAccuracyStats    CI stats for noise-corr shuffle
%         .NoisePermPvalue           p-value for noise-corr shuffle
%         .Cofiring                  [meanWithinCorr, meanAcrossCorr, difference]
%         .CofiringPerm              permutation distribution of corr differences
%         .CofiringStats             CI stats for corr-difference permutations
%         .CofiringPvalue            p-value for corr-difference permutation test
%
% -------------------------------------------------------------------------

% ----------------------------- Parse inputs ------------------------------
p = inputParser;
p.FunctionName = mfilename;

addRequired(p, 'FR', @(z) isnumeric(z) && ismatrix(z) && ~isempty(z));
addRequired(p, 'Condition', @(z) ~isempty(z) && numel(z) == size(FR,1));
addRequired(p, 'FeatureMask', @(z) isempty(z) || islogical(z) || (isnumeric(z) && isvector(z)));

addParameter(p, 'KFold', 10, @(z) isnumeric(z) && isscalar(z) && z>=2 && mod(z,1)==0);
addParameter(p, 'BoxConstraintRange', [1e-10, 1], @(z) isnumeric(z) && numel(z)==2 && all(z>0));
addParameter(p, 'MaxObjectiveEvaluations', 100, @(z) isnumeric(z) && isscalar(z) && z>=1 && mod(z,1)==0);
addParameter(p, 'AcquisitionFunctionName', 'expected-improvement-plus', @(z) ischar(z) || isstring(z));
addParameter(p, 'NumPermutations', 5000, @(z) isnumeric(z) && isscalar(z) && z>=1 && mod(z,1)==0);
addParameter(p, 'NumBootstraps', 10000, @(z) isnumeric(z) && isscalar(z) && z>=1 && mod(z,1)==0);
addParameter(p, 'Alpha', 0.05, @(z) isnumeric(z) && isscalar(z) && z>0 && z<1);
addParameter(p, 'UseParallel', true, @(z) islogical(z) && isscalar(z));
addParameter(p, 'ConditionValues', [], @(z) isempty(z) || numel(z)==2);

parse(p, FR, Condition, FeatureMask, varargin{:});

KFold                  = p.Results.KFold;
BoxConstraintRange     = p.Results.BoxConstraintRange;
MaxObjectiveEvaluations= p.Results.MaxObjectiveEvaluations;
AcquisitionFunctionName= string(p.Results.AcquisitionFunctionName);
numPermutations        = p.Results.NumPermutations;
numBootstraps          = p.Results.NumBootstraps;
alpha                  = p.Results.Alpha;
useParallel            = p.Results.UseParallel;
ConditionValues        = p.Results.ConditionValues;

% ----------------------------- Prepare data ------------------------------
% Z-score predictors (matches your code)
FR = zscore(FR);

if ~isempty(FeatureMask)
    if islogical(FeatureMask)
        if numel(FeatureMask) ~= size(FR,2)
            error('FeatureMask (logical) must have length = size(FR,2).');
        end
        FR(:, FeatureMask) = 0;
    else
        FR(:, FeatureMask) = 0;
    end
end

% ---------------------- Train model + hyperparameter opt -----------------
% Define the range for BoxConstraint (log scale)
boxConstraintVar = optimizableVariable('BoxConstraint', BoxConstraintRange, 'Transform', 'log');

Mdl = fitcsvm(FR, Condition, ...
    'Standardize', true, ...
    'KernelFunction', 'linear', ...
    'OptimizeHyperparameters', boxConstraintVar, ...
    'HyperparameterOptimizationOptions', struct( ...
        'AcquisitionFunctionName', char(AcquisitionFunctionName), ...
        'MaxObjectiveEvaluations', MaxObjectiveEvaluations));

% ---------------------- Cross-validation performance ---------------------
CVMdl = crossval(Mdl, 'KFold', KFold);
cvLoss = kfoldLoss(CVMdl);
cvAccuracy = 1 - cvLoss;

disp(['Cross-Validation Accuracy: ', num2str(cvAccuracy)]);

% ---------------------- Initialize output container ----------------------
LinearSVMConditionificationTable = struct();
LinearSVMConditionificationTable.SVMMdl     = Mdl;
LinearSVMConditionificationTable.CVMdl      = CVMdl;
LinearSVMConditionificationTable.cvLoss     = cvLoss;
LinearSVMConditionificationTable.cvAccuracy = cvAccuracy;

% ---------------------- Permutation test: labels -------------------------
permAccuracy = NaN(numPermutations, 1);

if useParallel
    parfor i = 1:numPermutations
        permutedCondition = Condition(randperm(length(Condition)));

        permMdl = fitcsvm(FR, permutedCondition, ...
            'Standardize', true, ...
            'KernelFunction', 'linear', ...
            'BoxConstraint', Mdl.ModelParameters.BoxConstraint);

        permCVMdl = crossval(permMdl, 'KFold', KFold);
        permLoss  = kfoldLoss(permCVMdl);
        permAccuracy(i) = 1 - permLoss;
    end
else
    for i = 1:numPermutations
        permutedCondition = Condition(randperm(length(Condition)));

        permMdl = fitcsvm(FR, permutedCondition, ...
            'Standardize', true, ...
            'KernelFunction', 'linear', ...
            'BoxConstraint', Mdl.ModelParameters.BoxConstraint);

        permCVMdl = crossval(permMdl, 'KFold', KFold);
        permLoss  = kfoldLoss(permCVMdl);
        permAccuracy(i) = 1 - permLoss;
    end
end

pValue = mean(permAccuracy >= cvAccuracy);
disp(['Original Cross-Validation Accuracy: ', num2str(cvAccuracy)]);
disp(['P-Value from Permutation Test: ', num2str(pValue)]);

ci = bootci(numBootstraps, {@mean, permAccuracy}, 'alpha', alpha);

LinearSVMConditionificationTable.ConditionPermAccuracy      = permAccuracy;
LinearSVMConditionificationTable.ConditionPermAccuracyStats = [ci(1) mean(permAccuracy) ci(2)];
LinearSVMConditionificationTable.ConditionPermPvalue        = pValue;

% ---------------------- Extract weights (coefficients) -------------------
weights = Mdl.Beta;
bias    = Mdl.Bias;
normWeights = weights / norm(weights, 2);

LinearSVMConditionificationTable.weights     = weights;
LinearSVMConditionificationTable.normWeights = normWeights;
LinearSVMConditionificationTable.bias        = bias;

% ---------------------- Permutation: heterogeneity across neurons --------
permAccuracy = NaN(numPermutations, 1);

if useParallel
    parfor i = 1:numPermutations
        permFR = arrayfun(@(row) FR(row, randperm(size(FR, 2))), 1:size(FR, 1), 'UniformOutput', false);
        permFR = cell2mat(permFR');

        permMdl = fitcsvm(permFR, Condition, ...
            'Standardize', true, ...
            'KernelFunction', 'linear', ...
            'BoxConstraint', Mdl.ModelParameters.BoxConstraint);

        permCVMdl = crossval(permMdl, 'KFold', KFold);
        permLoss  = kfoldLoss(permCVMdl);
        permAccuracy(i) = 1 - permLoss;
    end
else
    for i = 1:numPermutations
        permFR = arrayfun(@(row) FR(row, randperm(size(FR, 2))), 1:size(FR, 1), 'UniformOutput', false);
        permFR = cell2mat(permFR');

        permMdl = fitcsvm(permFR, Condition, ...
            'Standardize', true, ...
            'KernelFunction', 'linear', ...
            'BoxConstraint', Mdl.ModelParameters.BoxConstraint);

        permCVMdl = crossval(permMdl, 'KFold', KFold);
        permLoss  = kfoldLoss(permCVMdl);
        permAccuracy(i) = 1 - permLoss;
    end
end

pValue = mean(permAccuracy >= cvAccuracy);
ci = bootci(numBootstraps, {@mean, permAccuracy}, 'alpha', alpha);

LinearSVMConditionificationTable.ActPermAccuracy      = permAccuracy;
LinearSVMConditionificationTable.ActPermAccuracyStats = [ci(1) mean(permAccuracy) ci(2)];
LinearSVMConditionificationTable.ActPermPvalue        = pValue;

% ---------------------- Permutation: within functional subpops -----------
positiveWeights = find(normWeights > 0);
negativeWeights = find(normWeights < 0);

permAccuracy = NaN(numPermutations, 1);

if useParallel
    parfor i = 1:numPermutations
        permFR = FR;
        for j = 1:size(FR, 1)
            if ~isempty(positiveWeights)
                permFR(j, positiveWeights) = FR(j, positiveWeights(randperm(length(positiveWeights))));
            end
            if ~isempty(negativeWeights)
                permFR(j, negativeWeights) = FR(j, negativeWeights(randperm(length(negativeWeights))));
            end
        end

        permMdl = fitcsvm(permFR, Condition, ...
            'Standardize', true, ...
            'KernelFunction', 'linear', ...
            'BoxConstraint', Mdl.ModelParameters.BoxConstraint);

        permCVMdl = crossval(permMdl, 'KFold', KFold);
        permLoss  = kfoldLoss(permCVMdl);
        permAccuracy(i) = 1 - permLoss;
    end
else
    for i = 1:numPermutations
        permFR = FR;
        for j = 1:size(FR, 1)
            if ~isempty(positiveWeights)
                permFR(j, positiveWeights) = FR(j, positiveWeights(randperm(length(positiveWeights))));
            end
            if ~isempty(negativeWeights)
                permFR(j, negativeWeights) = FR(j, negativeWeights(randperm(length(negativeWeights))));
            end
        end

        permMdl = fitcsvm(permFR, Condition, ...
            'Standardize', true, ...
            'KernelFunction', 'linear', ...
            'BoxConstraint', Mdl.ModelParameters.BoxConstraint);

        permCVMdl = crossval(permMdl, 'KFold', KFold);
        permLoss  = kfoldLoss(permCVMdl);
        permAccuracy(i) = 1 - permLoss;
    end
end

pValue = mean(permAccuracy >= cvAccuracy);
ci = bootci(numBootstraps, {@mean, permAccuracy}, 'alpha', alpha);

LinearSVMConditionificationTable.MatchPermAccuracy      = permAccuracy;
LinearSVMConditionificationTable.MatchPermAccuracyStats = [ci(1) mean(permAccuracy) ci(2)];
LinearSVMConditionificationTable.MatchPermPvalue        = pValue;

% ---------------------- Permutation: destroy noise correlations ----------
permAccuracy = NaN(numPermutations, 1);

% Determine the two conditions used for within-Condition shuffling
if isempty(ConditionValues)
    u = unique(Condition);
    if numel(u) ~= 2
        error('Noise-correlation permutation requires exactly 2 unique Conditiones. Provide ''ConditionValues'' if needed.');
    end
    ConditionValues = u(:)';
end

if useParallel
    parfor i = 1:numPermutations
        permFR = FR;

        for cond = 1:2
            condIndices = find(contains(string(Condition), string(ConditionValues(cond))));
            for neuron = 1:size(FR, 2)
                permFR(condIndices, neuron) = FR(condIndices(randperm(length(condIndices))), neuron);
            end
        end

        permMdl = fitcsvm(permFR, Condition, ...
            'Standardize', true, ...
            'KernelFunction', 'linear', ...
            'BoxConstraint', Mdl.ModelParameters.BoxConstraint);

        permCVMdl = crossval(permMdl, 'KFold', KFold);
        permLoss  = kfoldLoss(permCVMdl);
        permAccuracy(i) = 1 - permLoss;
    end
else
    for i = 1:numPermutations
        permFR = FR;

        for cond = 1:2
            condIndices = find(contains(string(Condition), string(ConditionValues(cond))));
            for neuron = 1:size(FR, 2)
                permFR(condIndices, neuron) = FR(condIndices(randperm(length(condIndices))), neuron);
            end
        end

        permMdl = fitcsvm(permFR, Condition, ...
            'Standardize', true, ...
            'KernelFunction', 'linear', ...
            'BoxConstraint', Mdl.ModelParameters.BoxConstraint);

        permCVMdl = crossval(permMdl, 'KFold', KFold);
        permLoss  = kfoldLoss(permCVMdl);
        permAccuracy(i) = 1 - permLoss;
    end
end

pValue = mean(permAccuracy >= cvAccuracy);
ci = bootci(numBootstraps, {@mean, permAccuracy}, 'alpha', alpha);

LinearSVMConditionificationTable.NoisePermAccuracy      = permAccuracy;
LinearSVMConditionificationTable.NoisePermAccuracyStats = [ci(1) mean(permAccuracy) ci(2)];
LinearSVMConditionificationTable.NoisePermPvalue        = pValue;

% ---------------------- Cofiring comparison + permutation ----------------
positiveWeights = find(normWeights > 0);
negativeWeights = find(normWeights < 0);

corrMatrix = corr(FR);

withinPosCorr = corrMatrix(positiveWeights, positiveWeights);
withinNegCorr = corrMatrix(negativeWeights, negativeWeights);
acrossCorr    = corrMatrix(positiveWeights, negativeWeights);

WithinPosCorr  = withinPosCorr(triu(true(size(withinPosCorr)), 1));
WithinNegCorr  = withinNegCorr(triu(true(size(withinNegCorr)), 1));
meanWithinCorr = mean([WithinPosCorr; WithinNegCorr], 'omitnan');
meanAcrossCorr = mean(acrossCorr(:), 'omitnan');

origCorrDifference = meanWithinCorr - meanAcrossCorr;

LinearSVMConditionificationTable.Cofiring = [meanWithinCorr, meanAcrossCorr, origCorrDifference];

permCorrDifferences = zeros(numPermutations, 1);

for i = 1:numPermutations
    permutedSigns   = sign(normWeights);
    permutedSigns   = permutedSigns(randperm(length(permutedSigns)));
    permutedWeights = abs(normWeights) .* permutedSigns;

    permPosWeights = find(permutedWeights > 0);
    permNegWeights = find(permutedWeights < 0);

    permWithinPosCorr = corrMatrix(permPosWeights, permPosWeights);
    permWithinNegCorr = corrMatrix(permNegWeights, permNegWeights);
    permAcrossCorr    = corrMatrix(permPosWeights, permNegWeights);

    permWithinPosCorr = permWithinPosCorr(triu(true(size(permWithinPosCorr)), 1));
    permWithinNegCorr = permWithinNegCorr(triu(true(size(permWithinNegCorr)), 1));

    permMeanWithinCorr = mean([permWithinPosCorr; permWithinNegCorr], 'omitnan');
    permMeanAcrossCorr  = mean(permAcrossCorr(:), 'omitnan');

    permCorrDifferences(i) = permMeanWithinCorr - permMeanAcrossCorr;
end

pValue = mean(permCorrDifferences >= origCorrDifference);

disp(['Original Correlation Difference: ', num2str(origCorrDifference)]);
disp(['P-Value from Permutation Test: ', num2str(pValue)]);

% Bootstrap CI for permutation distribution of correlation differences
ci = bootci(1000, {@mean, permCorrDifferences}, 'alpha', alpha);

LinearSVMConditionificationTable.CofiringPerm   = permCorrDifferences;
LinearSVMConditionificationTable.CofiringStats  = [ci(1) mean(permCorrDifferences) ci(2)];
LinearSVMConditionificationTable.CofiringPvalue = pValue;

end
