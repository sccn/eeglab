classdef randomclassifier < classifier
%RANDOMCLASSIFIER returns a random classification
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: not supported by cvs2svn $
%

    methods
        function obj = randomclassifier(varargin)                  
           obj = obj@classifier(varargin{:});                      
        end
        function obj = train(obj,data,design)
            % does nothing
            
            if iscell(data), error('classifier does not take multiple datasets as input'); end

            if isnan(obj.nclasses), obj.nclasses = max(design(:,1)); end
            
       end
       function post = test(obj,data)       
           
           if iscell(data), error('classifier does not take multiple datasets as input'); end
                       
            % random classification
           post = rand(size(data,1),obj.nclasses);
           post = post ./ repmat(sum(post,2),[1 size(post,2)]);
       end

    end
end 
