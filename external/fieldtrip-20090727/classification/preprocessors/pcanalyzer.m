classdef pcanalyzer < standardizer
%PCANALYZER performs a principal component analysis
%
%   Options:
%   'proportion' : proportion of pc's or number of pc's. If < 1 then
%                 interpreted as a proportion of accounted variance; 
%                 otherwise as an absolute number; if empty 
%                 then all components are used (default = 0.80);
%
%   NOTE:
%   pcanalyzer inherits from standardizer since we normalize the data first
%
%   SEE ALSO:
%   princomp.m
%
%   REQUIRES:
%   statistics toolbox
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: not supported by cvs2svn $
%

    properties
       
        proportion = 0.80; % proportion of variance accounted for
        
        accvar; % cumulative variance accounted for per component
        pc; % principal components as column vectors
        ev; % eigenvalues of principal components
    end

    methods
    
        function obj = pcanalyzer(varargin)
           
          % check availability
          if ~license('test','statistics_toolbox')
            error('requires Matlab statistics toolbox');
          end

            obj = obj@preprocessor(varargin{:});                        
        end
        
        function obj = train(obj,data,design)

          obj = train@standardizer(obj,data,design);
          
          if iscell(data)

            obj.pc = cell(1,length(data));
            obj.accvar = cell(1,length(data));

            for c=1:length(data)

              [obj.pc{c},score,obj.ev{c}] = princomp(data{c});

              % proportion of the variance that is accounted for
              obj.accvar{c} = cumsum(obj.ev{c}/sum(obj.ev{c}));

              % determine how many principal components to use
              if ~isempty(obj.proportion)
                if obj.proportion >= 1
                  obj.pc{c} = obj.pc{c}(:,1:obj.proportion);
                else
                  obj.pc{c} = obj.pc{c}(:,1:min(1,find(obj.accvar{c} <= obj.proportion,1,'last')));
                end
              end
            end
          else

            [obj.pc,score,obj.ev] = princomp(data);

            % proportion of the variance that is accounted for
            obj.accvar = cumsum(obj.ev/sum(obj.ev));

            % determine how many principal components to use
            if ~isempty(obj.proportion)
              if obj.proportion >= 1
                obj.pc = obj.pc(:,1:obj.proportion);
              else
                obj.pc = obj.pc(:,1:min(1,find(obj.accvar <= obj.proportion,1,'last')));
              end
            end

          end

        end

        function data = test(obj,data)

          data = test@standardizer(obj,data);

          if iscell(data)

            for c=1:length(data)
              data{c} = (obj.pc{c}' * data{c}')';
            end

          else
            data = (obj.pc' * data')';
          end
        end

    end
end
