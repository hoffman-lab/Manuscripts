function Cs = perpl_CosineSimilarity(x,y)
%
% call:
%
%      Cs = perpl_CosineSimilarity(x)
%      Cs = perpl_CosineSimilarity(x,y)
%
% - If called with one argument (matrix x): returns cosine similarity
%   between all column pairs of x  ->  [nCols(x) x nCols(x)].
%
% - If called with two vectors x,y of same length: returns scalar cosine similarity.
%
% - If called with two matrices x,y (same number of rows): returns
%   cosine similarity between columns of x and columns of y
%   ->  [nCols(x) x nCols(y)].
%
% The interpretation of cosine similarity is analogous to that of a
% Pearson correlation.
%
% R.G. Bettinardi
% -----------------------------------------------------------------

% ----------------- Input checks & mode selection -------------------
if nargin == 1
    % Single input must be a matrix
    if ~ismatrix(x)
        error('Single input must be a matrix.');
    end
    Method = 1;  % self-similarity of columns in x

elseif nargin == 2
    % Two inputs: either both vectors OR both matrices
    if isvector(x) && isvector(y)
        if numel(x) ~= numel(y)
            error('Vectors x and y must have the same length.');
        end
        Method = 2;  % vector–vector similarity

    elseif ismatrix(x) && ismatrix(y)
        if size(x,1) ~= size(y,1)
            error('Matrices x and y must have the same number of rows.');
        end
        Method = 3;  % matrix–matrix (column-wise) similarity

    else
        error('x and y have to be either both vectors or both matrices.');
    end

else
    error('Function accepts either 1 or 2 input arguments.');
end

% ---------------------- Cosine similarity --------------------------
if Method == 1
    % Self-similarity between columns of x
    norms = sqrt(sum(x.^2, 1));          % 1 x nCols(x)
    normalizedMatrix = x ./ norms;       % M x nCols(x)
    Cs = normalizedMatrix' * normalizedMatrix;  % nCols(x) x nCols(x)

elseif Method == 2
    % Cosine similarity between two vectors
    x = x(:); y = y(:);                  % ensure column vectors
    xy   = dot(x,y);
    nx   = norm(x);
    ny   = norm(y);
    Cs   = xy / (nx * ny);

elseif Method == 3
    % Column-wise cosine similarity between two matrices
    % x: M x Nx, y: M x Ny  ->  Cs: Nx x Ny
    nx = sqrt(sum(x.^2, 1));             % 1 x Nx
    ny = sqrt(sum(y.^2, 1));             % 1 x Ny

    Xn = x ./ nx;                        % M x Nx (implicit expansion)
    Yn = y ./ ny;                        % M x Ny

    Cs = Xn.' * Yn;                      % Nx x Ny
end

end


% function Cs = perpl_CosineSimilarity(x,y)
% %
% % call:
% %
% %      Cs = getCosineSimilarity(x,y)
% %
% % Compute Cosine Similarity between vectors x and y.
% % x and y have to be of same length. The interpretation of
% % cosine similarity is analogous to that of a Pearson Correlation
% %
% % R.G. Bettinardi
% % -----------------------------------------------------------------
% if nargin == 1 && ismatrix(x)
%     Method = 1;
%     
% elseif nargin == 2 && isvector(x) && isvector(y)
%     Method = 2;
%     
% elseif ~ismatrix(x) || [~isvector(x) && ~isvector(y)]
%     error('x and y have to be 2 vectors or a matrix!')
%     
% end
% 
% 
% if Method == 1
%     
%     % Normalize each column to have unit norm
%     norms = sqrt(sum(x.^2, 1));
%     normalizedMatrix = x ./ norms;
%     
%     % Calculate the cosine similarity matrix
%     Cs = normalizedMatrix' * normalizedMatrix;
%     
% elseif Method == 2
%     xy   = dot(x,y);
%     nx   = norm(x);
%     ny   = norm(y);
%     nxny = nx*ny;
%     Cs   = xy/nxny;
%     
% end
% 
% end
% 
% 
