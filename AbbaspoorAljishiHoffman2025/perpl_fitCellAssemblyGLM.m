function CellAssemblyGLMModelTable = perpl_fitCellAssemblyGLM( ...
    Activation, Memage, SeqEl, TrialPeriod, Outcome, LinearVelocity, AngularVelocity, varargin)
% perpl_fitCellAssemblyGLM
% -------------------------------------------------------------------------
% PURPOSE
%   Fit a separate GLM for each cell assembly (each column in Activation)
%   using predictors related to memory age, sequence element, trial period,
%   behavioral outcome, and kinematics (linear/angular velocity).
%
%   The dependent variable for each GLM is:
%       SignificantActivation = Activation(:, CellNum)
%
%   For each CellNum:
%     1) Create a data table with predictors + SignificantActivation
%     2) Fit a GLM using fitglm() with a specified model formula
%     3) Store the fitted model object in a cell array
%
% INPUTS
%   Activation : numeric matrix, size (nObservations x nCells)
%       Dependent variable matrix. Each column corresponds to one cell assembly
%       (or one unit/feature) whose activation will be modeled across observations.
%
%   Memage : vector, length nObservations
%       Memory age label per observation. Converted to categorical internally.
%
%   SeqEl : vector, length nObservations
%       Sequence element per observation (expected values 1:4). Converted to ordinal
%       internally using: ordinal(SeqEl, {'1','2','3','4'}).
%
%   TrialPeriod : vector, length nObservations
%       Trial period/epoch label per observation. Converted to categorical internally.
%
%   Outcome : vector, length nObservations
%       Outcome label per observation. Converted to categorical internally.
%
%   LinearVelocity : numeric vector, length nObservations
%       Continuous linear velocity covariate per observation.
%
%   AngularVelocity : numeric vector, length nObservations
%       Continuous angular velocity covariate per observation.
%
% NAME-VALUE OPTIONAL INPUTS
%   'Session' : any (default: [])
%       Session identifier stored in output table (CellAssemblyGLMModelTable.Session).
%
%   'AnimalID' : any (default: [])
%       Animal identifier stored in output table (CellAssemblyGLMModelTable.AnimalID).
%
%   'Sess' : numeric scalar (default: 1)
%       Index at which results are stored inside CellAssemblyGLMModelTable fields.
%       This preserves your original pattern: CellAssemblyGLMModelTable.<field>{Sess} = ...
%
%   'TaskUnstableFR' : [] or logical/numeric vector (default: [])
%       Stored in output as CellAssemblyGLMModelTable.TaskUnstableFR{Sess}.
%       (No transformation is applied.)
%
%   'modelFormula' : char/string (default matches your code)
%       Model formula passed to fitglm().
%
%   'Distribution' : char/string (default: 'normal')
%       GLM distribution argument passed to fitglm().
%
%   'Link' : char/string (default: 'identity')
%       GLM link function passed to fitglm().
%
% OUTPUTS
%   CellAssemblyGLMModelTable : struct
%       Fields (cell arrays indexed by Sess):
%         .Session{Sess}        = Session (if provided)
%         .AnimalID{Sess}       = AnimalID (if provided)
%         .Model{Sess}          = column cell array of fitted GLM model objects (nCells x 1)
%         .TaskUnstableFR{Sess} = TaskUnstableFR (if provided)
%
% -------------------------------------------------------------------------

% ----------------------------- Parse inputs ------------------------------
p = inputParser;
p.FunctionName = mfilename;

addRequired(p, 'Activation', @(z) isnumeric(z) && ismatrix(z) && ~isempty(z));
addRequired(p, 'Memage', @(z) ~isempty(z) && numel(z) == size(Activation,1));
addRequired(p, 'SeqEl', @(z) ~isempty(z) && numel(z) == size(Activation,1));
addRequired(p, 'TrialPeriod', @(z) ~isempty(z) && numel(z) == size(Activation,1));
addRequired(p, 'Outcome', @(z) ~isempty(z) && numel(z) == size(Activation,1));
addRequired(p, 'LinearVelocity', @(z) isnumeric(z) && isvector(z) && numel(z) == size(Activation,1));
addRequired(p, 'AngularVelocity', @(z) isnumeric(z) && isvector(z) && numel(z) == size(Activation,1));

addParameter(p, 'Session', [], @(z) true);
addParameter(p, 'AnimalID', [], @(z) true);
addParameter(p, 'Sess', 1, @(z) isnumeric(z) && isscalar(z) && z>=1 && mod(z,1)==0);
addParameter(p, 'TaskUnstableFR', [], @(z) isempty(z) || islogical(z) || (isnumeric(z) && isvector(z)));

defaultFormula = "SignificantActivation ~ Memage + SeqEl + TrialPeriod + Outcome + " + ...
                 "LinearVelocity + AngularVelocity + Memage:SeqEl + Memage:TrialPeriod";
addParameter(p, 'modelFormula', defaultFormula, @(z) ischar(z) || isstring(z));
addParameter(p, 'Distribution', 'normal', @(z) ischar(z) || isstring(z));
addParameter(p, 'Link', 'identity', @(z) ischar(z) || isstring(z));

parse(p, Activation, Memage, SeqEl, TrialPeriod, Outcome, LinearVelocity, AngularVelocity, varargin{:});

Sess          = p.Results.Sess;
Session       = p.Results.Session;
AnimalID      = p.Results.AnimalID;
TaskUnstableFR= p.Results.TaskUnstableFR;
modelFormula  = char(p.Results.modelFormula);
distName      = char(p.Results.Distribution);
linkName      = char(p.Results.Link);

% ----------------------------- Prepare table struct ----------------------
CellAssemblyGLMModelTable = struct();
CellAssemblyGLMModelTable.Session        = cell(Sess, 1);
CellAssemblyGLMModelTable.AnimalID       = cell(Sess, 1);
CellAssemblyGLMModelTable.Model          = cell(Sess, 1);
CellAssemblyGLMModelTable.TaskUnstableFR = cell(Sess, 1);

CellAssemblyGLMModelTable.Session{Sess}  = Session;
CellAssemblyGLMModelTable.AnimalID{Sess} = AnimalID;

% ----------------------------- Type conversions --------------------------
% Match your code behavior:
Memage      = categorical(Memage);
SeqEl       = ordinal(SeqEl, {'1','2','3','4'});   % expects levels 1..4
TrialPeriod = categorical(TrialPeriod);
Outcome     = categorical(Outcome);

% ----------------------------- Fit models --------------------------------
nCells = size(Activation, 2);
Model  = cell(nCells, 1);

for CellNum = 1:nCells
    SignificantActivation = Activation(:, CellNum);

    % Construct the data table for fitglm
    dataTable = table(Memage, SeqEl, TrialPeriod, LinearVelocity(:), AngularVelocity(:), Outcome, SignificantActivation);

    % Fit GLM
    glmModel = fitglm(dataTable, modelFormula, ...
        'Distribution', distName, ...
        'Link', linkName);

    % Store fitted model
    Model{CellNum} = glmModel;
end

CellAssemblyGLMModelTable.Model{Sess} = Model; % (nCells x 1) cell array
CellAssemblyGLMModelTable.TaskUnstableFR{Sess} = TaskUnstableFR;

end
