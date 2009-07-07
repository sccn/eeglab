function ydata = tsne_p(P, labels, no_dims)
%TSNE_P Performs symmetric t-SNE on affinity matrix P
%
%   mappedX = tsne_p(P, labels, no_dims)
%
% The function performs symmetric t-SNE on pairwise similarity matrix P 
% to create a low-dimensional map of no_dims dimensions (default = 2).
% The matrix P is assumed to be symmetric, sum up to 1, and have zeros
% on the diagonal.
% The labels of the data are not used by t-SNE itself, however, they 
% are used to color intermediate plots. Please provide an empty labels
% matrix [] if you don't want to plot results during the optimization.
% The low-dimensional data representation is returned in mappedX.
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
    
    % Initialize some variables
    n = size(P, 1);                                     % number of instances
    momentum = 0.5;                                     % initial momentum
    final_momentum = 0.8;                               % value to which momentum is changed
    mom_switch_iter = 250;                              % iteration at which momentum is changed
    max_iter = 1000;                                    % maximum number of iterations
    epsilon = 500;                                      % initial learning rate
    min_gain = .01;
    
    % Make sure P-vals are set properly
    P(1:n + 1:end) = 0;
    P = P ./ sum(P(:));
    P = max(P, eps);
    P = P * 4;		% lie about the P-vals
    
    % Initialize the solution
    ydata = .0001 * randn(n, no_dims);
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
    