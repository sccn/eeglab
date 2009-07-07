function [network, err] = tsne_v_backprop(network, train_X, train_labels, test_X, test_labels, max_iter, perplexity)
%TSNE_V_BACKPROP Perform finetuning of network using the t-SNE gradient
%
%   [network, err] = tsne_v_backprop(network, train_X, train_labels,
%   test_X, test_labels, max_iter, perplexity)
%
% Perform finetuning of the specified network using the t-SNE gradient. The
% finetuning is performed for max_iter iterations (default = 10). The
% Gaussians employed in the high-dimensional space have the specified
% perplexity (default = 30).
%
%
% (C) Laurens van der Maaten
% Maastricht University, 2008


    if ~exist('max_iter', 'var') || isempty(max_iter)
        max_iter = 30;        
    end
    if ~exist('perplexity', 'var') || isempty(perplexity)
        perplexity = 30;
    end
    
    % Initialize some variables
    n = size(train_X, 1);
    batch_size = min([5000 n]);
    ind = randperm(n);
    v = 1;
    
    % Precompute joint probabilities for all batches
    disp('Precomputing P-values...');
    curX = cell(floor(n ./ batch_size), 1);
    P = cell(floor(n ./ batch_size), 1);
    i = 1;
    for batch=1:batch_size:n          
        if batch + batch_size - 1 <= n
            curX{i} = double(train_X(ind(batch:min([batch + batch_size - 1 n])),:));    % select batch
            P{i} = x2p(curX{i}, perplexity, 1e-5);                                      % compute affinities using fixed perplexity
            P{i}(isnan(P{i})) = 0;                                                      % make sure we don't have NaN's
            P{i} = (P{i} + P{i}') / 2;                                                  % make symmetric
            P{i} = P{i} ./ sum(P{i}(:));                                                % obtain estimation of joint probabilities
            P{i} = max(P{i}, eps);
            i = i + 1;
        end
    end

    % Run the optimization
    err = zeros(max_iter, 1);
    for iter=1:max_iter
        
        % Run for all mini-batches
        disp(['Iteration ' num2str(iter) '...']);
        b = 1;
        for batch=1:batch_size:n
            if batch + batch_size - 1 <= n

                % Construct current solution
                x = v;
                for i=1:length(network)
                    x = [x; network{i}.W(:); network{i}.bias_upW(:)];
                end

                % Perform conjugate gradient using three linesearches
                x = minimize(x, 'tsne_v_grad', 3, curX{b}, P{b}, network);
                b = b + 1;                
                
                % Store new solution
                v = x(1);
                ii = 2;
                for i=1:length(network)
                    network{i}.W = reshape(x(ii:ii - 1 + numel(network{i}.W)), size(network{i}.W)); 
                    ii = ii + numel(network{i}.W);
                    network{i}.bias_upW = reshape(x(ii:ii - 1 + numel(network{i}.bias_upW)), size(network{i}.bias_upW));
                    ii = ii + numel(network{i}.bias_upW);
                end
            end
        end        
                    
        % Estimate the current error
        activations = run_data_through_network(network, curX{1});
        sum_act = sum(activations .^ 2, 2);
        Q = (1 + (bsxfun(@plus, sum_act, bsxfun(@plus, sum_act', -2 * activations * activations')) ./ v)) .^ -((v + 1) / 2);
        Q(1:n+1:end) = 0;
        Q = Q ./ sum(Q(:));
        Q = max(Q, eps);
        C = sum(sum(P{1} .* log((P{1} + eps) ./ (Q + eps))));
        disp(['Degrees of freedom: ' num2str(v)]);
        disp(['t-SNE error: ~' num2str(C)]);
        
        % Compute current 1-NN error        
        err(iter) = knn_error(run_data_through_network(network, train_X), train_labels, ...
                              run_data_through_network(network,  test_X), test_labels, 1);
        disp(['1-NN  error: ' num2str(err(iter))]);
    end
    