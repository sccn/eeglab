function [mappedX, landmarks, costs] = fast_tsne(X, no_dims, initial_dims, landmarks, perplexity)
%FAST_TSNE Runs the fast Intel (IPP) implementation of t-SNE
%
%   [mappedX, landmarks, costs] = fast_tsne(X, no_dims, initial_dims, landmarks, perplexity)
%
% Runs the fast implementation Diffusion Stochastic Neighbor Embedding 
% algorithm. The high-dimensional datapoints are specified by X. First, the
% dimensionality of the datapoints is reduced to initial_dims dimensions
% using PCA (default = 30). Then, Diffusion SNE reduces the points to 
% no_dims dimensions (default = 2). The percentage of points to use as 
% landmarks may be specified through landmarks (0 <= landmarks <= 1), or 
% you may specify the indices of the landmark points in a vector. If not 
% specified, both the number and the amount of landmarks points is 
% determined automatically.  The used perplexity in the Gaussian kernel can 
% be set through the perplexity variable (default = 30). Note that the 
% perplexity is mainly of influence on small datasets where landmarks are 
% not employed.
% The function returns the low-dimensional datapoints in mappedX, the used 
% landmark points in landmarks, and the cost contribution per datapoint in
% costs.
%
%
% (C) Laurens van der Maaten
% Maastricht University, 2008

    if ~exist('no_dims', 'var') || isempty(no_dims)
        no_dims = 2;
    end
    if ~exist('initial_dims', 'var') || isempty(initial_dims)
        initial_dims = 30;
    end
    if ~exist('landmarks', 'var') || isempty(landmarks)
        if size(X, 1) < 6000
            landmarks = 1;
        else
            landmarks = 6000 / size(X, 1);
        end
    end
    if ~exist('perplexity', 'var')
        perplexity = 30;
    end
    
    % Perform the initial dimensionality reduction using PCA
    X = bsxfun(@minus, X, mean(X, 1));
    covX = X' * X;
    [M, lambda] = eig(covX);
    [lambda, ind] = sort(diag(lambda), 'descend');
    if initial_dims > size(M, 2)
        initial_dims = size(M, 2);
    end
    M = M(:,ind(1:initial_dims));
    X = X * M;
    clear covX M lambda
    
    % Run the fast diffusion SNE implementation
    write_data(X, no_dims, landmarks, perplexity);
    if ispc
        system('tSNE.exe');
    elseif ismac
        system('./tSNE_maci');
    else
      if strcmp(computer,'GLNX86')
        system('./tSNE_linux_32');
      else % assume 64 bit
        system('./tSNE_linux_64');
      end
    end
    [mappedX, landmarks, costs] = read_data;   
    landmarks = landmarks + 1;              % correct for Matlab indexing
end


% Writes the datafile for the fast diffusion SNE implementation
function write_data(X, no_dims, landmarks, perplexity)
    [n, d] = size(X);
    h = fopen('data.dat', 'wb');
	fwrite(h, n, 'integer*4');
	fwrite(h, d, 'integer*4');
    fwrite(h, no_dims, 'integer*4');
    fwrite(h, perplexity, 'double');
    if numel(landmarks) == 1
        fwrite(h, landmarks, 'double');
    else
        fwrite(h, numel(landmarks), 'double');
    end
	fwrite(h, X', 'double');
    if numel(landmarks) > 1
        fwrite(h, landmarks - 1, 'integer*4');
    end
	fclose(h);
end


% Reads the resultfile from the fast diffusion SNE implementation
function [X, landmarks, costs] = read_data
    h = fopen('result.dat', 'rb');
	n = fread(h, 1, 'integer*4');
	d = fread(h, 1, 'integer*4');
	X = fread(h, n * d, 'double');
    landmarks = fread(h, n, 'integer*4');
    costs = fread(h, n, 'double');
    X = reshape(X, [d n])';
	fclose(h);
end
