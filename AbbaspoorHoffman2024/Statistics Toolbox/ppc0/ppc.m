function ppc0 = ppc(phases, varargin)
%% PPC Compute pairwise phase consistency (PPC0 of Vinck et al, 2010)
% INPUT:
%  - phases: Vector of phases (radians, single-precision)
% OUTPUT:
%  - ppc0: pairwise phase consistency of the phases in `phases`.
% OPTIONS:
%  - UseParallel: When to use the parallel pool: 'always', 'never', or 'auto'
%      If auto, the parallel pool is used whenever there are at least
%      `ParallelThreshold` items in vector and a pool is available or
%      can be launched (see `LaunchPool`).
%  - ParallelThreshold: Minimum number of items to use parallel pool in (auto mode only)
%  - LaunchPool: If true, launch a parallel pool if necessary (auto mode only)
% NOTES:
% - Both versions are *much* faster than the naive implementation of the
%    PPC formula and a bit more numerically stable.
%
% Matt Krause <mrkrause@gmail.com>

ip = inputParser();
ip.addRequired('phases', @(x) isvector(x) && isnumeric(x));
ip.addParameter('UseParallel', 'always', @(x) ismember(x, {'never', 'always', 'auto'}));
ip.addParameter('ParallelThreshold', 8000, @(x) isnumeric(x) && isscalar(x));
ip.addParameter('LaunchPool', true, @(x) x || ~x);
ip.parse(phases, varargin{:});

auto_use_parallel = strcmpi(ip.Results.UseParallel, 'auto') && ...
    (length(phases) >= ip.Results.ParallelThreshold) && ...
    (~isempty(gcp('nocreate')) || ip.Results.LaunchPool);

N = length(phases);
scale_factor = (2/(N*(N-1)));

% [sphases, cphases] = sincos(phases);
cphases = cos(phases);
sphases = sin(phases);

sphases = double(sphases); %% Force to double or you will rapidly get overflow
cphases = double(cphases);

if strcmpi(ip.Results.UseParallel, 'always') || auto_use_parallel
    ppc0 = nan(N-1, 1);
    parfor ii=1:(N-1) % Break-even point for parfor is about 10k samples, turn it off if your data is usually smaller
        ppc0(ii) = (sum((cphases(ii) .* cphases((ii+1):N)) + (sphases(ii).* sphases((ii+1):N))));
    end
    ppc0 = sum(ppc0) * scale_factor;
else
    ppc0 = 0;
    for ii=1:(N-1) % Break-even point for parfor is about 10k samples, turn it on if your data is usually larger
        ppc0 = ppc0 + (scale_factor * (sum((cphases(ii) .* cphases((ii+1):N)) + (sphases(ii).*sphases((ii+1):N)))));
    end
end
