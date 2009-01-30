function [R]=test_sc(CC,D,mode,classlabel)
% TEST_SC: apply statistical and SVM classifier to test data 
%
%  R = test_sc(CC,D,TYPE [,target_Classlabel]) 
%       R.output     output distance for each class
%       R.classlabel class for output data
%  The target class is optional. If it is provided, the following values are returned. 
%       R.kappa Cohen's kappa coefficient
%       R.ACC   Classification accuracy 
%       R.H     Confusion matrix 
%
% The classifier CC is typically obtained by TRAIN_SC. If a statistical 
% classifier is used, TYPE can be used to modify the classifier. 
%    TYPE = 'MDA'    mahalanobis distance based classifier
%    TYPE = 'MD2'    mahalanobis distance based classifier
%    TYPE = 'MD3'    mahalanobis distance based classifier
%    TYPE = 'GRB'    Gaussian radial basis function 
%    TYPE = 'QDA'    quadratic discriminant analysis
%    TYPE = 'LD2'    linear discriminant analysis (see LDBC2)
%    TYPE = 'LD3'    linear discriminant analysis (see LDBC3)
%    TYPE = 'LD4'    linear discriminant analysis (see LDBC4)
%    TYPE = 'GDBC'   general distance based classifier
% 
% see also: TRAIN_SC
%
% References: 
% [1] R. Duda, P. Hart, and D. Stork, Pattern Classification, second ed. 
%       John Wiley & Sons, 2001. 

%	$Id: test_sc.m,v 1.1 2009-01-30 06:04:49 arno Exp $
%	Copyright (C) 2005,2006,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the  License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

if nargin<3, 
        mode = [];
end;
[t1,t] = strtok(CC.datatype,':');
[t2,t] = strtok(t,':');
[t3,t] = strtok(t,':');
if ~strcmp(t1,'classifier'), return; end; 

if isfield(CC,'prewhite')
        D = D*CC.prewhite(2:end,:) + CC.prewhite(ones(size(D,1),1),:);
        CC = rmfield(CC,'prewhite');
end;

POS1 = [strfind(CC.datatype,'/gsvd'),strfind(CC.datatype,'/sparse')];

if 0,


elseif strcmp(CC.datatype,'classifier:nbpw')
	error('NBPW not implemented yet')
	%%%% Naive Bayesian Parzen Window Classifier %%%%
        d = repmat(NaN,size(D,1),size(CC.MEAN,1));
	for k = 1:size(CC.MEAN,1)
		z = (D - CC.MEAN(repmat(k,size(D,1),1),:)).^2 ./ (CC.VAR(repmat(k,size(D,1),1),:));
		z = z + log(CC.VAR(repmat(k,size(D,1),1),:)); % + log(2*pi);
		d(:,k) = sum(-z/2, 2) + log(mean(CC.N(k,:)));
	end; 
	d = exp(d-log(mean(sum(CC.N,1)))-log(2*pi)/2);


elseif strcmp(CC.datatype,'classifier:nbc')
	%%%% Naive Bayesian Classifier %%%%
        d = repmat(NaN,size(D,1),size(CC.MEAN,1));
	for k = 1:size(CC.MEAN,1)
		z = (D - CC.MEAN(repmat(k,size(D,1),1),:)).^2 ./ (CC.VAR(repmat(k,size(D,1),1),:));
		z = z + log(CC.VAR(repmat(k,size(D,1),1),:)); % + log(2*pi);
		d(:,k) = sum(-z/2, 2) + log(mean(CC.N(k,:)));
	end; 
	d = exp(d-log(mean(sum(CC.N,1)))-log(2*pi)/2);


elseif strcmp(CC.datatype,'classifier:anbc')
	%%%% Augmented Naive Bayesian Classifier %%%%
        d = repmat(NaN,size(D,1),size(CC.MEAN,1));
	for k = 1:size(CC.MEAN,1)
		z = (D*CC.V - CC.MEAN(repmat(k,size(D,1),1),:)).^2 ./ (CC.VAR(repmat(k,size(D,1),1),:));
		z = z + log(CC.VAR(repmat(k,size(D,1),1),:)); % + log(2*pi);
		d(:,k) = sum(-z/2, 2) + log(mean(CC.N(k,:)));
	end; 
	d = exp(d-log(mean(sum(CC.N,1)))-log(2*pi)/2);


elseif strcmp(CC.datatype,'classifier:statistical:rda')
	% Friedman (1989) Regularized Discriminant analysis
	if isfield(CC,'hyperparameters') && isfield(CC.hyperparameters,'lambda')  && isfield(CC.hyperparameters,'gamma')
	        D = [ones(size(D,1),1),D];  % add 1-column
		lambda = CC.hyperparameters.lambda;
		gamma  = CC.hyperparameters.gamma;
	        d = repmat(NaN,size(D,1),size(CC.MD,1));
                ECM = CC.MD./CC.NN; 
                NC = size(ECM); 
                ECM0 = squeeze(sum(ECM,1));  %decompose ECM
                [M0,sd,COV0,xc,N,R2] = decovm(ECM0);
                for k = 1:NC(1);
                        [M,sd,s,xc,N,R2] = decovm(squeeze(ECM(k,:,:)));
                	s = ((1-lambda)*N*s+lambda*COV0)/((1-lambda)*N+lambda);
                	s = (1-gamma)*s+gamma*(trace(s))/(NC(2)-1)*eye(NC(2)-1);
                        ir  = [-M;eye(NC(2)-1)]*inv(s)*[-M',eye(NC(2)-1)];  % inverse correlation matrix extended by mean
                        d(:,k) = -sum((D*ir).*D,2); % calculate distance of each data point to each class
                end;
	else 
		error('QDA: hyperparamters lambda and/or gamma not defined')
	end; 	         
	 


elseif strcmp(CC.datatype,'classifier:csp')
	d = (D*CC.csp_w).^2;
	d = filtfilt(CC.FiltB,CC.FiltA,(D*CC.csp_w).^2);
	R = test_sc(CC.CSP,log(d));	% LDA classifier of 
	d = R.output; 
	cl= R.classlabel; 


elseif strcmp(CC.datatype,'classifier:svm:lib:1vs1') | strcmp(CC.datatype,'classifier:svm:lib:rbf');
        
        [cl, accuracy] = svmpredict(classlabel, D, CC.model);   %Use the classifier

        %Create a pseudo tsd matrix for bci4eval
        d = zeros(size(cl,1), CC.model.nr_class);
        for i = 1:size(cl,1)
                d(i,cl(i)) = 1;
        end
        
        
elseif isfield(CC,'weights'); %strcmpi(t2,'svm') | (strcmpi(t2,'statistical') & strncmpi(t3,'ld',2)) ;
        % linear classifiers like: LDA, SVM, LPM 
        %d = [ones(size(D,1),1), D] * CC.weights;
        d = repmat(NaN,size(D,1),size(CC.weights,2));
        for k = 1:size(CC.weights,2),
                d(:,k) = D * CC.weights(2:end,k) + CC.weights(1,k);
        end;        
        if size(CC.weights,2)==1,
                d = [d, -d];
        end;

        
elseif ~isempty(POS1)	% GSVD & sparse
        CC.datatype = CC.datatype(1:POS1(1)-1);
        r = test_sc(CC,D*CC.G);
        d = r.output; 


elseif strcmp(t2,'statistical');
        if isempty(mode)
                mode.TYPE = upper(t3); 
        end;
        D = [ones(size(D,1),1),D];  % add 1-column

        if 0, 
        elseif strcmpi(mode.TYPE,'LD2'),
                %d = ldbc2(CC,D);
                ECM = CC.MD./CC.NN; 
                NC = size(ECM); 
                ECM0 = squeeze(sum(ECM,1));  %decompose ECM
                [M0,sd,COV0,xc,N,R2] = decovm(ECM0);
                for k = 1:NC(1);
                        ecm = squeeze(ECM(k,:,:));
                        [M1,sd,COV1,xc,N,R2] = decovm(ECM0-ecm);
                        [M2,sd,COV2,xc,N,R2] = decovm(ecm);
                        w     = (COV1+COV2)\(M2'-M1')*2;
                        w0    = -M0*w;
                        W(:,k) = [w0; w];
                end;
                d = D*W;
        elseif strcmpi(mode.TYPE,'LD3');
                %d = ldbc3(CC,D);
                ECM = CC.MD./CC.NN; 
                NC = size(ECM); 
                ECM0 = squeeze(sum(ECM,1));  %decompose ECM
                [M0,sd,COV0,xc,N,R2] = decovm(ECM0);
                for k = 1:NC(1);
                        ecm = squeeze(ECM(k,:,:));
                        [M1,sd,COV1,xc,N,R2] = decovm(ECM0-ecm);
                        [M2,sd,COV2,xc,N,R2] = decovm(ecm);
                        w     = COV0\(M2'-M1')*2;
                        w0    = -M0*w;
                        W(:,k) = [w0; w];
                end;
                d = D*W;
        elseif strcmpi(mode.TYPE,'LD4');
                %d = ldbc4(CC,D);
                ECM = CC.MD./CC.NN; 
                NC = size(ECM); 
                ECM0 = squeeze(sum(ECM,1));  %decompose ECM
                [M0,sd,COV0,xc,N,R2] = decovm(ECM0);
                for k = 1:NC(1);
                        ecm = squeeze(ECM(k,:,:));
                        [M1,sd,COV1,xc,N1,R2] = decovm(ECM0-ecm);
                        [M2,sd,COV2,xc,N2,R2] = decovm(ecm);
                        w     = (COV1*N1+COV2*N2)\((M2'-M1')*(N1+N2));
                        w0    = -M0*w;
                        W(:,k) = [w0; w];
                end;
                d = D*W;
        elseif strcmpi(mode.TYPE,'MDA');
                for k = 1:length(CC.IR);
                        d(:,k) = -sum((D*CC.IR{k}).*D,2); % calculate distance of each data point to each class
                end;
        elseif strcmpi(mode.TYPE,'MD2');
                for k = 1:length(CC.IR);
                        d(:,k) = sum((D*CC.IR{k}).*D,2); % calculate distance of each data point to each class
                end;
                d = -sqrt(d);
        elseif strcmpi(mode.TYPE,'GDBC');
                for k = 1:length(CC.IR);
                        d(:,k) = sum((D*CC.IR{k}).*D,2) + CC.logSF7(k); % calculate distance of each data point to each class
                end;
                d = exp(-d/2);
        elseif strcmpi(mode.TYPE,'MD3');
                for k = 1:length(CC.IR);
                        d(:,k) = sum((D*CC.IR{k}).*D,2) + CC.logSF7(k); % calculate distance of each data point to each class
                end;
                d = exp(-d/2);
                d = d./repmat(sum(d,2),1,size(d,2));  % Zuordungswahrscheinlichkeit [1], p.601, equ (18.39)
        elseif strcmpi(mode.TYPE,'QDA');     
                for k = 1:length(CC.IR);
                        % [1] (18.33) QCF - quadratic classification function  
                        d(:,k) = -(sum((D*CC.IR{k}).*D,2) - CC.logSF5(k)); 
                end;
        elseif strcmpi(mode.TYPE,'GRB');     % Gaussian RBF
                for k = 1:length(CC.IR);
                        d(:,k) = sum((D*CC.IR{k}).*D,2); % calculate distance of each data point to each class
                end;
                d = exp(-sqrt(d)/2);
        elseif strcmpi(mode.TYPE,'GRB2');     % Gaussian RBF
                for k = 1:length(CC.IR);
                        d(:,k) = sum((D*CC.IR{k}).*D,2); % calculate distance of each data point to each class
                end;
                d = exp(-d);
        elseif strcmpi(mode.TYPE,'MQU');     % Multiquadratic 
                for k = 1:length(CC.IR);
                        d(:,k) = sum((D*CC.IR{k}).*D,2); % calculate distance of each data point to each class
                end;
                d = -sqrt(1+d);
        elseif strcmpi(mode.TYPE,'IMQ');     % Inverse Multiquadratic 
                for k = 1:length(CC.IR);
                        d(:,k) = sum((D*CC.IR{k}).*D,2); % calculate distance of each data point to each class
                end;
                d = (1+d).^(-1/2);
        elseif strcmpi(mode.TYPE,'Cauchy');     % Cauchy RBF
                for k = 1:length(CC.IR);
                        d(:,k) = sum((D*CC.IR{k}).*D,2); % calculate distance of each data point to each class
                end;
                d = 1./(1+d);
        end;
else
        fprintf(2,'Error TEST_SC: unknown classifier\n');
        return;
end;

[tmp,cl] = max(d,[],2);
cl = CC.Labels(cl); 
cl(isnan(tmp)) = NaN; 

R.output = d; 
R.classlabel = cl; 

if nargin>3,
        [R.kappa,R.sd,R.H,z,R.ACC] = kappa(classlabel(:),cl(:));
end;
