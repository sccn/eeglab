function machine = train_lin_rbm(X, h, eta, max_iter, weight_cost)
%TRAIN_RBM Trains a Restricted Boltzmann Machine with Gaussian hiddens
%
%   machine = train_lin_rbm(X, h, eta, max_iter, weight_cost)
%
% Trains a first-order Restricted Boltzmann Machine on dataset X. The RBM
% has h hidden nodes (default = 20). The training is performed by means of
% the contrastive divergence algorithm. The activation functions that
% is applied in the visible and hidden layers are binary stochastic.
% In the training of the RBM, the learning rate is determined by eta 
% (default = 0.001). The maximum number of iterations can be specified 
% through max_iter (default = 30). The variable weight_cost sets the amount
% of weight decay that is employed (default = 0.0002).
% The trained RBM is returned in the machine struct.
%
%
% (C) Laurens van der Maaten
% Maastricht University, 2008


    % Process inputs
    if ischar(X)
        load(X);
    end
    if ~exist('h', 'var') || isempty(h)
        h = 20;
    end
    if ~exist('eta', 'var') || isempty(eta)
        eta = 0.001;       
    end
    if ~exist('max_iter', 'var') || isempty(max_iter)
        max_iter = 50;        
    end
    if ~exist('weight_cost', 'var') || isempty(weight_cost)
        weight_cost = 0.0002;
    end  
    
    % Other parameters
    initial_momentum = 0.5;                 % momentum for first five iterations
    final_momentum = 0.9;                   % momentum for remaining iterations
    
    % Initialize some variables
    [n, v] = size(X);
    batch_size = 100;
    W = randn(v, h) * 0.1;
    bias_upW = zeros(1, h);
    bias_downW = zeros(1, v);
    deltaW = zeros(v, h);
    deltaBias_upW = zeros(1, h);
    deltaBias_downW = zeros(1, v);   
    
    % Main loop
    for iter=1:max_iter
        
        % Set momentum
        err = 0;
        ind = randperm(n); 
        if iter <= 5
            momentum = initial_momentum;
        else
            momentum = final_momentum;
        end
        
        % Run for all mini-batches
        for batch=1:batch_size:n          
            if batch + batch_size <= n
            
                % Set values of visible nodes
                vis1 = double(X(ind(batch:min([batch + batch_size - 1 n])),:));

                % Compute probabilities for hidden nodes
                hid1 = bsxfun(@plus, vis1 * W, bias_upW);

                % Sample states for hidden nodes
                hid_states = hid1 + randn(size(hid1));
                    
                % Compute probabilities for visible nodes
                vis2 = 1 ./ (1 + exp(-(bsxfun(@plus, hid_states * W', bias_downW))));

                % Compute probabilities for hidden nodes
                hid2 = bsxfun(@plus, vis2 * W, bias_upW);

                % Now compute the weights update
                posprods = vis1' * hid1;
                negprods = vis2' * hid2;
                deltaW = momentum * deltaW + eta * (((posprods - negprods) / batch_size) - (weight_cost * W));
                deltaBias_upW   = momentum * deltaBias_upW   + (eta / batch_size) * (sum(hid1, 1) - sum(hid2, 1));
                deltaBias_downW = momentum * deltaBias_downW + (eta / batch_size) * (sum(vis1, 1) - sum(vis2, 1));
                              
                % Update the network weights
                W          = W          + deltaW;
                bias_upW   = bias_upW   + deltaBias_upW;
                bias_downW = bias_downW + deltaBias_downW;
                
                % Estimate the error
                err = err + sum(sum((vis1 - vis2) .^ 2));
            end
        end 
        
        % Estimation of reconstruction error
        disp(['Iteration ' num2str(iter) ': error ~' num2str(err ./ n)]);
    end
    
    % Return RBM
    machine.W = W;
    machine.bias_upW = bias_upW;
    machine.bias_downW = bias_downW;
    disp(' ');
