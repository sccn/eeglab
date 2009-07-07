classdef clfproc < clfbase
%CLFPROC classification procedure class
%   
%   A classification procedure is just a sequence of classification methods {method1 method2 method3 ...} that are
%   called in this order and where the output of the previous method acts
%   as input to the next method.
%
%   A classification method always has a train and test function for the
%   two different modes of operation. A classification procedure is
%   just a sequence of methods {method1 method2 method3 ...} that are
%   called in this order and where the output of the previous method acts
%   as input to the next method. Note that method derives from handle
%   so any changes are permanent. Methods:
%
%   PREPROCESSOR: produces transformed data
%   FEATURESELECTOR: selects and/or extracts features and produces data
%   CLASSIFIER: produces posterior probabilities
%   REGRESSOR: produces real-valued predictions
%   OPTIMIZER: wrapper class that optimizes parameters of another method
%   
%   the output of a classifier can be used as input to:
%
%   evaluate: takes posteriors and labels and produces an evaluation measure
%
%   parameters of a method should be set by means of a constructor call
%   method(obj,varargin)
%
%   repeated calls to train should allow for online learning as much as
%   possible
%   
%   data is given by an N x M matrix with N examples/times and 
%   M features
%   design is an N x 1 matrix with the N class labels 1,2,3,...
%
%   multiple datasets are entered using cell arrays of data/design
%
%   Options:
%   'verbose' : verbose output [false]
%   'randomize' : randomize trials before training [false]
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: not supported by cvs2svn $
%

    properties

        clfmethods; % the methods that specify the classification procedure
        
        randomize = false;
    end
    
    methods
       function obj = clfproc(clfmethods,varargin)
       % constructor expects classification methods

        if ~nargin
            % if unspecified use a random classifier
            
            if obj.verbose
              fprintf('classifier not specified; using random classifier');
            end
            
            clfmethods = { randomclassifier() };
        end
        
        % cast to cell array if only a classifier/regressor is specified
        if ~iscell(clfmethods)
          clfmethods = { clfmethods };
        end

        if ~isempty(clfmethods)

            % check if this is a valid classification procedure:
            if ~(all(cellfun(@(x)(isa(x,'clfbase')),clfmethods)))
                error('invalid classification procedure; check contents of clfproc');
            end

            % at least one classification or regression method at the end            
            if ~(isa(clfmethods{end},'classifier') || isa(clfmethods{end},'regressor') || isa(clfmethods{end},'optimizer'))
              error('classification procedure should end with a classifier or a regressor');
            end
            
        end

        obj.clfmethods = clfmethods;
       
       end
       
       function obj = train(obj,data,design)
            % train just calls the methods' train functions in order to produce a posterior

            % check for missing data
            if obj.verbose && obj.missing(data)
                warning('data contains missing values');
            end
            
            % check for missing labels
            if obj.verbose && obj.missing(design)
              warning('design contains missing values');
            end
            
            % make sure data is a (cell array of) matrix
            if (iscell(data) && ndims(data{1}) > 2) || ndims(data) > 2
                data = obj.collapse(data);
            end
            
            if obj.randomize
              [data,design] = obj.shuffle(data,design);
            end
                        
            for c=1:length(obj.clfmethods)      
                obj.clfmethods{c} = obj.clfmethods{c}.train(data,design);
                data = obj.clfmethods{c}.test(data);                
            end        
                   
       end
       
       function data = test(obj,data)
           % test just calls the methods' test functions in order to produce a posterior

           if isempty(data)
             data = [];
             return;
           end

           % check for missing data
            if obj.verbose && obj.missing(data)
                warning('data contains missing values');
            end
                       
           % make sure data is a (cell array of) matrix
           if (iscell(data) && ndims(data{1}) > 2) || ndims(data) > 2
               data = obj.collapse(data);
           end

           for c=1:length(obj.clfmethods)               
               data = obj.clfmethods{c}.test(data);
           end

       end       
       
       function clf = classify(obj,data)
           % calls test functions and converts posterior into
           % classifications
          
           [tmp,clf] = max(obj.test(data),[],2);
       end
       
       function s = name(obj)
         % returns classification procedure as a string
         s = '{ ';
         for j=1:length(obj.clfmethods)
           s = [s class(obj.clfmethods{j}) ' '];
         end
         s = strcat(s,' }');
         
       end
       
       function m = getmodel(obj,label,dims)
         % return the model implied by this classification procedure
    
         if nargin < 2, label = 1; end      
         if nargin < 3, dims = []; end
    
         % get the final model
         m = obj.clfmethods{end}.getmodel(label,dims);
         
         % if not yet handled by the method itself
         if ~isempty(dims) && ~all(size(m) == dims) 
           
           if iscell(m) % multiple models for multiple datasets

             for c=1:length(m)
               if prod(dims) == numel(m{c});
                 m{c} = reshape(m{c},dims);
               else
                 if obj.verbose
                   fprintf('original dimensions and new dimensions do not match\n');
                 end
               end
             end
             
           else
             
             if prod(dims) == numel(m)
               m = reshape(m,dims);
             else
               if obj.verbose
                 fprintf('original dimensions and new dimensions do not match\n');
               end
             end
             
           end
           
         end

       end
                    
    end
   
    methods(Static)
    
    function b = missing(data)
      % check if there are nans in the data
      
      if iscell(data)
        b = false;
        for c=1:length(data)
          if any(isnan(data{c}(:)))
            b = true;
            return;
          end
        end
      else
        b = any(isnan(data(:)));
      end
    end
    
  end
    
end 

