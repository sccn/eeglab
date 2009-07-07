classdef pnn < classifier
%PNN probabilistic neural network classifier
%
%   Options
%   'spread' : newpnn parameter that will be optimized if it is a vector
%
%   SEE ALSO:
%   newpnn.m
%
%   NOTE:
%   The FieldTrip replacement of dist.m does not work in this context!
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: not supported by cvs2svn $
%

    properties
        net; % the neural network
        spread = 0.1; % modifies nn behaviour; if a vector, will be optimized
    end

    methods
       function obj = pnn(varargin)
                  
           obj = obj@classifier(varargin{:});
           
            % check availability
           if ~license('test','neural_network_toolbox')
               error('requires Matlab neural network toolbox');
           end
           
       end
       function obj = train(obj,data,design)
                                  
           if iscell(data), error('PNN does not take multiple datasets as input'); end                       
                                
            if isnan(obj.nclasses), obj.nclasses = max(design(:,1)); end
            
            obj.net = newpnn(data',design(:,1)',obj.spread);
                     
       end
       function post = test(obj,data)       
           
           if iscell(data), error('NB does not take multiple datasets as input'); end

           post = sim(obj.net,data')';
       end

    end
end 
