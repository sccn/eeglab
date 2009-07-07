classdef loocrossvalidator < crossvalidator
%LOOCROSSVALIDATOR leave-one-out crossvalidation class
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: not supported by cvs2svn $
%

   methods
       function obj = loocrossvalidator(varargin)
                                 
           obj = obj@crossvalidator(varargin{:});
           
           % leave one out crossvalidation
           obj.cvfolds = inf;           
       end
       
   end
end
