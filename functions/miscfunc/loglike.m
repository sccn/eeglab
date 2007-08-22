% loglike() - log likehood function to estimate dependence between components
%
% Usage: f = loglike(W, S);
% 
% Computation of the log-likelihood function under the model
% that the ICs are 1/cosh(s) distributed (according to the tanh
% nonlinearity in ICA). It does not exactly match for the logistic
% nonlinearity, but should be a decent approximation
%
% negative log likelihood function
% f = -( log(abs(det(W))) - sum(sum(log( cosh(S) )))/N - M*log(pi) );
%
% With these meanings:
% W: total unmixing matrix, ie, icaweights*icasphere
% S: 2-dim Matrix of source estimates
% N: number of time points
% M: number of components
%
% Author: Arnaud Delorme and Jorn Anemuller

function f=loglike(W, S);

    W = double(W);
    S = double(S);
    % normalize activities
    stds = std(S, [], 2);
    S = S./repmat(stds, [1 size(S,2)]);
    W = W./repmat(stds, [1 size(W,2)]);
    
    M = size(W,1);
    if ndims(S) == 3
        S = reshape(S, size(S,1), size(S,3)*size(S,2)); 
    end;
    N = size(S,2);
    
    % detect infinite and remove them
    tmpcoh = log( cosh(S) );
    tmpinf = find(isinf(tmpcoh));
    tmpcoh(tmpinf) = [];
    N = (N*M-length(tmpinf))/M;
    
    f=-( log(abs(det(W))) - sum(sum(tmpcoh))/N - M*log(pi) );

