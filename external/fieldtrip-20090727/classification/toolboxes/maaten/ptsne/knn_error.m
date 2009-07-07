function err = knn_error(train_X, train_labels, test_X, test_labels, k)
%KNN_ERROR Compute k-nearest neighbor error
%
%   err = knn_error(train_X, train_labels, test_X, test_labels, k)
%
% Measures the generalization error of a k-nearest neighbor classifier that
% is trained train_X and train_labels, and tested on test_X and 
% test_labels. The number of nearest neighbors can be specified through k 
% (default = 1).
%
%
% (C) Laurens van der Maaten
% Maastricht University, 2008


    if ~exist('k', 'var') || isempty(k)
        k = 1;
    end

    % Compute pairwise distance matrix
    sum_train = sum(train_X .^ 2, 2);
    sum_test  = sum(test_X  .^ 2, 2);
    D = bsxfun(@plus, sum_train, bsxfun(@plus, sum_test', -2 * train_X * test_X'));

    % Perform labeling
    classification = zeros(size(test_X, 1), 1);
    for j=1:size(test_X, 1)
	[foo, ii] = sort(D(:,j), 1, 'ascend');
        tmp1 = train_labels(ii(1:k));
        tmp2 = unique(tmp1);
        counts = zeros(length(tmp2), 1);
        for h=1:length(tmp2)
            counts(h) = sum(ismember(tmp1, tmp2(h))) + sum(.01 ./ find(tmp1 == tmp2(h)));
        end
        [foo, index] = max(counts);
        classification(j) = tmp2(index);
    end

    % Evaluate error
    err = sum(test_labels ~= classification) ./ numel(test_labels);
    