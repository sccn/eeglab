classdef classifier < predictor
%CLASSIFIER abstract classifier method class
%
% A classifier takes a variable number of arguments upon construction. 
% During operation, the classifier takes data and
% produces posterior probabilities of class labels as an N x C matrix for N
% examples and C classes.
% 
% Subclasses should implement the train and test functions and optionally
% the getmodel function which reshapes parameters to something
% interpretable
%
% PROPERTIES
%   'nclasses'  : number of classes  (determined from data)
%   'nfeatures' : number of features (determined from data)
%   'nexamples' : number of examples (determined from data)
%
% SEE ALSO:
%   combiner.m
%   da.m
%   dynamic_classifier.m
%   ensemble.m
%   gnb.m
%   gp.m
%   gslr.m
%   gslr_transfer.m
%   hgnb_transfer.m
%   hmm.m
%   kernelmethod.m
%   libsvm.m
%   lr.m
%   lvq.m
%   mixtureclassifier.m
%   nb.m
%   nearestneighbour.m
%   one_against_one.m
%   one_against_rest.m
%   pnn.m
%   randomclassifier.m
%   regdif.m
%   regressor.m
%   rfda.m
%   rlda.m
%   rnb.m
%   static_classifier.m
%   svmmethod.m
%   transfer_classifier.m
%   blogreg.m
%
% Copyright (c) 2008, Marcel van Gerven
%
% $Log: not supported by cvs2svn $
%    
    
    methods
        
        function obj = classifier(varargin)
     
          % parse options
          for i=1:2:length(varargin)
            obj.(varargin{i}) = varargin{i+1};
          end

        end        
        
        function clf = predict(obj,data)
           % convert posteriors into classifications
           
           [tmp,clf] = max(obj.test(data),[],2);
        end
        
    end

    methods(Access = protected)
      
      function obj = parse_input(obj,data,design)
        
        obj.nexamples = size(design,1);
        obj.nfeatures = size(data,2);
        
      end
      
    end
    
end 
