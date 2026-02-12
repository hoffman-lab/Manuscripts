function pvalue_matrix = perpl_fisherZtest(task_matrix, sleep_matrix, n1, n2)
    % Check if the matrices are of the same size
    if size(task_matrix) ~= size(sleep_matrix)
        error('Task and Sleep matrices must be of the same size');
    end
    
    % Number of elements in the matrices
    [num_rows, num_cols] = size(task_matrix);
    
    % Initialize p-value matrix
    pvalue_matrix = zeros(num_rows, num_cols);
    
    % Fisher's Z transformation
    z_task_matrix = atanh(task_matrix);
    z_sleep_matrix = atanh(sleep_matrix);
    
    % Calculate the standard error
    se = sqrt(1/(n1 - 3) + 1/(n2 - 3));
    
    % Perform the Z-Test for each pairwise correlation
    for i = 1:num_rows
        for j = 1:num_cols
            % Calculate the Z-score for the difference between correlations
            z_diff = abs(z_task_matrix(i,j) - z_sleep_matrix(i,j)) / se;
            
            % Calculate the p-value from the Z-score
            pvalue_matrix(i,j) = 2 * (1 - normcdf(z_diff, 0, 1)); % Two-tailed test
            
            if i == j
                pvalue_matrix(i,j) = 1;
            end
        end
    end
end
