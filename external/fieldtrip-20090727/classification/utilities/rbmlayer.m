classdef rbmlayer
%
% Code provided by Ruslan Salakhutdinov and Geoff Hinton
% Permission is granted for anyone to copy, use, modify, 
% or distribute this program and accompanying programs and 
% documents for any purpose, provided this copyright notice
% is retained and prominently displayed, along with a note 
% saying that the original programs are available from our web page. 
% The programs and documents are distributed without any warranty, 
% express or implied. As the programs were written for research 
% purposes only, they have not been tested to the degree that would 
% be advisable in any important application. All use of these 
% programs is entirely at the user's own risk.
%
% http://www.cs.toronto.edu/~hinton/MatlabForSciencePaper.html
%

    properties

        numhid = 10;            % number of hidden units
        
        maxepoch = 5;           % maximum number of epochs
        epsilonw      = 0.1;    % Learning rate for weights
        epsilonvb     = 0.1;    % Learning rate for biases of visible units
        epsilonhb     = 0.1;    % Learning rate for biases of hidden units
        weightcost  = 0.0002;
        initialmomentum  = 0.5;
        finalmomentum    = 0.9;

        visbiases               % biases of visible (layer 1) units
        hidbiases               % biases of hidden (layer 2) units

        vishid                  % symmetric weights
        
        verbose = false;
    end
    
    properties (Hidden = true)

        % momentum terms
        vishidinc;
        hidbiasinc;
        visbiasinc;        
    end
    
    methods
        
        function obj = rbmlayer(varargin)
           
            % parse options
            for i=1:2:length(varargin)
                obj.(varargin{i}) = varargin{i+1};
            end
            
        end
        
        function [obj,batchposhidprobs] = train(obj,batchdata)
            % train RBM on batchdata
            
            % This program trains Restricted Boltzmann Machine in which
            % visible, binary, stochastic pixels are connected to
            % hidden, binary, stochastic feature detectors using symmetrically
            % weighted connections. Learning is done with 1-step Contrastive Divergence.

            % The function returns:
            % visbiases -- the biases of the visible (lower layer) units
            % hidbiases -- the biases of the hidden (upper layer) units
            % vishid -- the visible-to-hidden weights
            % batchposhidprobs -- the hidden layer activations which can act as
            %   data on a subsequent RBM.
            
            [numcases numdims numbatches]=size(batchdata);
            batchposhidprobs=zeros(numcases,obj.numhid,numbatches);

            if obj.verbose, fprintf('pretraining layer with RBM: %d-%d\n',numdims,obj.numhid); end

            %%%% initialization 
            
            numhid = obj.numhid;

            obj.vishid     = 0.1*randn(numdims, numhid);
            obj.hidbiases  = zeros(1,numhid);
            obj.visbiases  = zeros(1,numdims);

            obj.vishidinc  = zeros(numdims,numhid);
            obj.hidbiasinc = zeros(1,numhid);
            obj.visbiasinc = zeros(1,numdims);

            poshidprobs = zeros(numcases,numhid);
            neghidprobs = zeros(numcases,numhid);
            posprods    = zeros(numdims,numhid);
            negprods    = zeros(numdims,numhid);

            for epoch = 1:obj.maxepoch

                if obj.verbose, fprintf(1,'epoch %d\r',epoch); end
                
                errsum=0;
                for batch = 1:numbatches,
                
                    if obj.verbose, fprintf(1,'epoch %d batch %d\r',epoch,batch); end

                    %%%%%%%%% START POSITIVE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
                    data = batchdata(:,:,batch);
                    poshidprobs = 1./(1 + exp(-data*obj.vishid - repmat(obj.hidbiases,numcases,1)));
                    batchposhidprobs(:,:,batch)=poshidprobs;
                    posprods    = data' * poshidprobs;
                    poshidact   = sum(poshidprobs);
                    posvisact   = sum(data);

                    %%%%%%%%% END OF POSITIVE PHASE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    poshidstates = poshidprobs > rand(numcases,numhid);

                    %%%%%%%%% START NEGATIVE PHASE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    negdata = 1./(1 + exp(-poshidstates*obj.vishid' - repmat(obj.visbiases,numcases,1)));
                    neghidprobs = 1./(1 + exp(-negdata*obj.vishid - repmat(obj.hidbiases,numcases,1)));
                    negprods  = negdata'*neghidprobs;
                    neghidact = sum(neghidprobs);
                    negvisact = sum(negdata);

                    %%%%%%%%% END OF NEGATIVE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    err= sum(sum( (data-negdata).^2 ));
                    errsum = err + errsum;

                    if epoch>5,
                        momentum=obj.finalmomentum;
                    else
                        momentum=obj.initialmomentum;
                    end;

                    %%%%%%%%% UPDATE WEIGHTS AND BIASES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    obj.vishidinc = momentum*obj.vishidinc + ...
                        obj.epsilonw*( (posprods-negprods)/numcases - obj.weightcost*obj.vishid);
                    obj.visbiasinc = momentum*obj.visbiasinc + (obj.epsilonvb/numcases)*(posvisact-negvisact);
                    obj.hidbiasinc = momentum*obj.hidbiasinc + (obj.epsilonhb/numcases)*(poshidact-neghidact);

                    obj.vishid = obj.vishid + obj.vishidinc;
                    obj.visbiases = obj.visbiases + obj.visbiasinc;
                    obj.hidbiases = obj.hidbiases + obj.hidbiasinc;

                    %%%%%%%%%%%%%%%% END OF UPDATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                end

                if 1||obj.verbose, fprintf(1, 'epoch %4i error %6.1f  \n', epoch, errsum); end
            
            end

        end
        
        function hidprobs = test(obj,data)
            % propagate activations to get hidden layer activations

            hidprobs = 1./(1 + exp(-data*obj.vishid - repmat(obj.hidbiases,size(data,1),1)));
        end

    end
end
