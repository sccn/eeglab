function ydata = tsne_d(D, labels, no_dims, perplexity)
%TSNE_D Performs symmetric t-SNE on the pairwise Euclidean distance matrix D
%
%   mappedX = tsne_d(D, labels, no_dims, perplexity)
%   mappedX = tsne_d(D, labels, initial_solution, perplexity)
%
% The function performs symmetric t-SNE on the NxN pairwise Euclidean 
% distance matrix D to construct an embedding with no_dims dimensions 
% (default = 2). An initial solution obtained from an other dimensionality 
% reduction technique may be specified in initial_solution. 
% The perplexity of the Gaussian kernel that is employed can be specified 
% through perplexity (default = 30). The labels of the data are not used 
% by t-SNE itself, however, they are used to color intermediate plots. 
% Please provide an empty labels matrix [] if you don't want to plot 
% results during the optimization.
% The data embedding is returned in mappedX.
%
%
% (C) Laurens van der Maaten, 2008
% Maastricht University


    if ~exist('labels', 'var')
        labels = [];
    end
    if ~exist('no_dims', 'var') || isempty(no_dims)
        no_dims = 2;
    end
    if ~exist('perplexity', 'var') || isempty(perplexity)
        perplexity = 30;
    end
    
    % First check whether we already have an initial solution
    if numel(no_dims) > 1
        initial_solution = true;
        ydata = no_dims;
        no_dims = size(ydata, 2);
        perplexity = no_dims;
    else
        initial_solution = false;
    end
    
    % Initialize some variables
    n = size(D, 1);                                     % number of instances
    momentum = 0.5;                                     % initial momentum
    final_momentum = 0.8;                               % value to which momentum is changed
    mom_switch_iter = 250;                              % iteration at which momentum is changed
    max_iter = 1000;                                    % maximum number of iterations
    epsilon = 500;                                      % initial learning rate
    min_gain = .01;
    
    % Compute joint probabilities
    P = d2p(D .^ 2, perplexity, 1e-5);                                          % compute affinities using fixed perplexity
    P = 0.5 * (P + P');                                                         % make symmetric
    P = P ./ sum(P(:));                                                         % obtain estimation of joint probabilities
    P = max(P, eps);
    P = P * 4;                                                                  % prevent local minima by lying about P-vals
    
    % Initialize the solution
    if ~initial_solution
        ydata = .0001 * randn(n, no_dims);
    end
    y_grads = zeros(size(ydata));
    y_incs  = zeros(size(ydata));
    gains = ones(size(ydata));
    
    % Run the iterations
    for iter=1:max_iter
        
        % Compute joint probability that point i and j are neighbors
        sum_ydata = sum(ydata .^ 2, 2);
        num = 1 ./ (1 + bsxfun(@plus, sum_ydata, bsxfun(@plus, sum_ydata', -2 * ydata * ydata')));            % Student-t distribution
        num(1:n+1:end) = 0;                                                     % set diagonal to zero
        Q = num ./ sum(num(:));                                                 % normalize to get probabilities
        Q = max(Q, eps); 
        
        % Compute the gradients
        stiffnesses = 4 * (P - Q) .* num;
        for i=1:n
            y_grads(i,:) = sum(bsxfun(@times, bsxfun(@minus, ydata(i,:), ydata), stiffnesses(:,i)), 1);
        end
            
        % Update the solution
        gains = (gains + .2) .* (sign(y_grads) ~= sign(y_incs)) ...            % note that the y_grads are actually -y_grads
              + (gains * .8) .* (sign(y_grads) == sign(y_incs));
        gains(gains < min_gain) = min_gain;
        y_incs = momentum * y_incs - epsilon * (gains .* y_grads);
        ydata = ydata + y_incs;
        ydata = bsxfun(@minus, ydata, mean(ydata, 1));
        
        % Update the momentum if necessary
        if iter == mom_switch_iter
            momentum = final_momentum;
        end
        if iter == 100
            P = P ./ 4;
        end
        
        % Print out progress
        if ~rem(iter, 10)
            cost = sum(sum(P .* log((P + eps) ./ (Q + eps)), 2));
            disp(['Iteration ' num2str(iter) ': error is ' num2str(cost)]);
        end
        
        % Display scatter plot (maximally first three dimensions)
        if ~rem(iter, 10) && ~isempty(labels)
            if no_dims == 1
                scatter(ydata, ydata, 9, labels, 'filled');
            elseif no_dims == 2
                scatter(ydata(:,1), ydata(:,2), 9, labels, 'filled');
            else
                scatter3(ydata(:,1), ydata(:,2), ydata(:,3), 40, labels, 'filled');
            end
            axis tight
            axis off
            drawnow
        end
    end
    