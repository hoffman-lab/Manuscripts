function [R, lags] = assembly_xcorr(ReactivationStrength, maxLag)
% ReactivationStrength: (nAsm x T)
% maxLag: scalar (e.g. 1000)
%
% R: (nAsm x nAsm x (2*maxLag+1))
% lags: 1 x (2*maxLag+1)

X = double(ReactivationStrength);
[nAsm, ~] = size(X);

% -------- preprocessing (recommended for assemblies) --------
X = X - mean(X, 2);   % remove baseline

nLag = 2*maxLag + 1;
R = zeros(nAsm, nAsm, nLag, 'single');

for i = 1:nAsm
    xi = X(i,:);

    % autocorrelation (optional, but cheap and useful)
    [R(i,i,:), lags] = xcorr(xi, xi, maxLag, 'coeff');

    for j = i+1:nAsm
        xj = X(j,:);

        % compute once
        r = xcorr(xi, xj, maxLag, 'coeff');

        % store
        R(i,j,:) = r;
        R(j,i,:) = flip(r);  % exploit symmetry
    end
end
end
