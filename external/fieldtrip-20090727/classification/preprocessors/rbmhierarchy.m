classdef rbmhierarchy < preprocessor
%RBMHIERARCHY creates a hierarchy of restricted Boltzmann machines
%
% OPTIONS
% 'rbms'      : the RBM layers as a cell array
% 'nbatches'  : number of splits of the data when training the RBM 
%
% EXAMPLE
%  myproc = { ...    
%    standardizer() ...
%    rbmhierarchy('nbatches',10,'rbms',{rbmlayer('numhid',100,'maxepoch',3) rbmlayer('numhid',10,'maxepoch',3)}) ...
%    gslr('maxgroup',10) ...
%    };
%
% SEE ALSO
%     rbmlayer.m
%
%   Copyright (c) 2008, Marcel van Gerven
%
%   $Log: not supported by cvs2svn $
%

    properties
        rbms            % set of rbm layers

        nbatches = 10;  % number of batches        
    end

    methods
        function obj = rbmhierarchy(varargin)
            
            obj = obj@preprocessor(varargin{:});           
            if isempty(obj.rbms), error('rbm layers must be specified'); end   
            
        end
        
        function obj = train(obj,data,design)
            % create hierarchy of rbms

            if obj.nbatches > size(data,1)
                nbatches = size(data,1);
            else
                nbatches = obj.nbatches;
            end
            
            % convert to batch data
            bsize = floor(size(data,1)./nbatches);
            batchdata = zeros(bsize,size(data,2),nbatches);
            for j=1:nbatches
                batchdata(:,:,j) = data(((j-1)*bsize+1):(j*bsize),:);
            end

            % pretrain layers
            if obj.verbose, fprintf('pretraining\n'); end

            for c=1:length(obj.rbms)
                [obj.rbms{c},batchdata] = obj.rbms{c}.train(batchdata);
            end
        end
        
        function data = test(obj,data)            
            % propagate activations and save features as examples
          
            features = cell(1,length(obj.rbms)+1);
            sz = zeros(1,length(obj.rbms)+1);
            
            features{1} = data; sz(1) = size(data,2);         
            for c=1:length(obj.rbms)
                features{c+1} = obj.rbms{c}.test(features{c});
                sz(c+1) = size(features{c+1},2);
            end
            sz = cumsum(sz);
            
            newdata = zeros(size(data,1),sz(end));
            newdata(:,1:sz(1)) = data;
            for c=2:length(obj.rbms)
               newdata(:,(sz(c-1)+1):sz(c)) = features{c}; 
            end
        end
        
    end
end 
