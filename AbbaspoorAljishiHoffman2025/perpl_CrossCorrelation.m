function [lags, correlations] = perpl_CrossCorrelation(x, y, maxLag)
    % Validate input lengths
    if length(x) ~= length(y)
        error('Input vectors x and y must be of the same length.');
    end

    % Initialize variables
    n = length(x);
    correlations = zeros(2 * maxLag + 1, 1);
    lags = -maxLag:maxLag;
    
    % Calculate correlation for each lag
    for i = -maxLag:maxLag
        if i < 0
            x_lagged = x(1:end+i);
            y_lagged = y(1-i:end);
        else
            x_lagged = x(i+1:end);
            y_lagged = y(1:end-i);
        end
        
        % Compute Pearson correlation coefficient
        correlations(i + maxLag + 1) = corr(x_lagged', y_lagged');
    end
end
