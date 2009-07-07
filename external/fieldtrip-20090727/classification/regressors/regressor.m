classdef regressor < clfbase
%REGRESSOR abstract regression method class
%
% A regressor takes a variable number of arguments upon construction. 
% During operation, the regressor takes data and
% produces predicted responses as an N x 1 matrix for N examples.
% 
% Subclasses should implement the train and test functions.
%
% OPTIONS
%   'nfeatures' : number of features (determined from data)
%   'nexamples' : number of examples (determined from data)
%
% SEE ALSO
%   circreg.m
%   linreg.m
%   mulreg.m
%   blinreg.m
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: not supported by cvs2svn $
%    
    
    properties
             
        nfeatures = nan;
        nexamples = nan;
        
    end
    
    methods
        
        function obj = regressor(varargin)
               
          % parse options
          for i=1:2:length(varargin)
            obj.(varargin{i}) = varargin{i+1};
          end

        end
        
    end

    methods(Abstract)
        
        obj = train(obj,data,design);
        res = test(obj,data);    
    end
    
end 
