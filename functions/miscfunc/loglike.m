% loglike() - log likehood function to estimate dependence between components
%
% Usage: f = loglike(W, S);
% 
% Computation of the log-likelihood function under the model
% that the ICs are 1/cosh(s) distributed (according to the tanh
% nonlinearity in ICA). It does not exactly match for the logistic
% nonlinearity, but should be a decent approximation -jorn email
%
% negative log likelihood function
% f = -( log(abs(det(W))) - sum(sum(log( cosh(S) )))/N - M*log(pi) );
%
% With these meanings:
% W: total unmixing matrix, ie, icaweights*icasphere
% S: 2-dim Matrix of source estimates
% N: number of time points
% M: number of components

function f=loglike(W, S);
    
    M = size(W,1);
    N = size(S,2);
    if ndims(S) == 3
        S = reshape(S, size(S,1), size(S,3)*size(S,2)); 
    end;
    f=-( log(abs(det(W))) - sum(sum(log( cosh(S) )))/N - M*log(pi) );

