function mappedX = run_data_through_network(network, X)
%RUN_DATA_THROUGH_NETWORK Run data through the network
%
%   mappedX = run_data_through_network(network, X)
%
% Runs the dataset X through the parametric t-SNE embedding defined in
% network. The result is returned in mappedX.
%
%
% (C) Laurens van der Maaten
% Maastricht University, 2008

    if isstruct(network)
        network = {network};
    end

    % Run the data through the network
    n = size(X, 1);
    mappedX = [X ones(n, 1)];
    for i=1:length(network) - 1
        mappedX = [1 ./ (1 + exp(-(mappedX * [network{i}.W; network{i}.bias_upW]))) ones(n, 1)];
    end
    mappedX = mappedX * [network{end}.W; network{end}.bias_upW];
    