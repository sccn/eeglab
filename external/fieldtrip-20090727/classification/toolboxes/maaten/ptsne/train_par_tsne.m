function [network, err] = train_par_tsne(train_X, train_labels, test_X, test_labels, layers, training)
%TRAIN_PAR_TSNE Trains a parametric t-SNE embedding
%
%   [network, err] = train_par_tsne(train_X, train_labels, test_X,
%   test_labels, layers, training)
%
% Trains up a parametric t-SNE embedding with the structure that is
% specified in layers. The used training technique is specified in
% training. Possible values are 'CD1' and 'PCD' (default = 'CD1').
%
%
% (C) Laurens van der Maaten
% Maastricht University, 2008


    if ~exist('training', 'var') || isempty(training)
        training = 'CD1';
    end

    % Pretrain the network
    origX = train_X;
    no_layers = length(layers);
    network = cell(1, no_layers);
    for i=1:no_layers

        % Print progress
        disp(['Training layer ' num2str(i) ' (size ' num2str(size(train_X, 2)) ' -> ' num2str(layers(i)) ')...']);
        
        if i ~= no_layers
          
            % Train layer using binary units
            if strcmp(training, 'CD1')
                network{i} = train_rbm(train_X, layers(i));
            elseif strcmp(training, 'PCD')
                network{i} = train_rbm_pcd(train_X, layers(i));
            elseif strcmp(training, 'None')
                v = size(train_X, 2);
                network{i}.W = randn(v, layers(i)) * 0.1;
                network{i}.bias_upW = zeros(1, layers(i));
                network{i}.bias_downW = zeros(1, v);
            else
                error('Unknown training procedure.');
            end
                
            % Transform data using learned weights
            train_X = 1 ./ (1 + exp(-(bsxfun(@plus, train_X * network{i}.W, network{i}.bias_upW))));
        else
            
            % Train layer using Gaussian hidden units
            if ~strcmp(training, 'None')
                network{i} = train_lin_rbm(train_X, layers(i));
            else
                v = size(train_X, 2);
                network{i}.W = randn(v, layers(i)) * 0.1;
                network{i}.bias_upW = zeros(1, layers(i));
                network{i}.bias_downW = zeros(1, v);
            end
        end
    end
        
    % Perform backpropagation of the network using t-SNE gradient
    [network, err] = tsne_backprop(network(1:no_layers), origX, train_labels, test_X, test_labels, 30, 30, 1);
    