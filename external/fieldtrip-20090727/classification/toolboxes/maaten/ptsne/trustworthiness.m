function T = trustworthiness(X, mappedX, k)
%TRUSTWORTHINESS Computes the trustworthiness of a low-D embedding
%
%   T1 = trustworthiness(X, mappedX, k)
%
% Computes the trustworthiness values T1 and T2 of the low-dimensional 
% embedding specified in mappedX. The original data should be specified in
% X.


    % Compute pairwise distance matrices
    sumX = sum(X .^ 2, 2);
    hD = bsxfun(@plus, sumX, bsxfun(@plus, sumX', -2 * X * X'));
    sumX = sum(mappedX .^ 2, 2);
    lD = bsxfun(@plus, sumX, bsxfun(@plus, sumX', -2 * mappedX * mappedX'));
    
    % Compute neighborhood indices
    [hD, ind1] = sort(hD, 2, 'ascend');
    clear hD
    [lD, ind2] = sort(lD, 2, 'ascend');
    clear lD

    % Compute thrustworthiness values
    n = size(X, 1);
    T = 0;
    ranks = zeros(k, 1);
    for i=1:n
        for j=1:k
            ranks(j) = find(ind1(i,:) == ind2(i, j + 1));
        end
        ranks = ranks - k - 1;
        T = T + sum(ranks(ranks > 0));
    end
    T = 1 - ((2 / (n * k * (2 * n - 3 * k - 1))) * T);
    