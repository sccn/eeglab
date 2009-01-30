function [CC]=train_sc(D,classlabel,MODE)
% Train a (statistical) classifier
% 
%  CC = train_sc(D,classlabel)
%  CC = train_sc(D,classlabel,MODE)
%
% CC contains the model parameters of a classifier which can be applied 
%   to test data using test_sc. 
%   R = test_sc(CC,D,...) 
%
%  The following classifier types are supported MODE.TYPE
%    'MDA'      mahalanobis distance based classifier [1]
%    'MD2'      mahalanobis distance based classifier [1]
%    'MD3'      mahalanobis distance based classifier [1]
%    'GRB'      Gaussian radial basis function     [1]
%    'QDA'      quadratic discriminant analysis    [1]
%    'LD2'      linear discriminant analysis (see LDBC2) [1]
%               MODE.hyperparameter.gamma: regularization parameter [default 0] 
%    'LD3'      linear discriminant analysis (see LDBC3) [1]
%               MODE.hyperparameter.gamma: regularization parameter [default 0] 
%    'LD4'      linear discriminant analysis (see LDBC4) [1]
%               MODE.hyperparameter.gamma: regularization parameter [default 0] 
%    'LD5'      another LDA (motivated by CSP)
%               MODE.hyperparameter.gamma: regularization parameter [default 0] 
%    'RDA'      regularized discriminant analysis [7]
%               MODE.hyperparameter.gamma: regularization parameter 
%               MODE.hyperparameter.lambda =  
%		gamma = 0, lambda = 0 : MDA 
%		gamma = 0, lambda = 1 : LDA 
% 		Hint: hyperparameters are used only in test_sc.m, testing different 
%		the hyperparameters do not need repetitive calls to train_sc, 
%		it is sufficient to modify CC.hyperparameters before calling test_sc. 	
%    'GDBC'     general distance based classifier  [1]
%    ''         statistical classifier, requires Mode argument in TEST_SC	
%    '###/GSVD'	GSVD and statistical classifier [2,3], 
%    '###/sparse'  sparse  [5] 
%               '###' must be 'LDA' or any other classifier 
%    'SVM','SVM1r'  support vector machines, one-vs-rest
%               MODE.hyperparameter.c_value = 
%    'PSVM'	Proximal SVM [8] 
%               MODE.hyperparameter.nu  (default: 1.0)
%    'PLS'	(linear) partial least squares regression 
%    'REG'      regression analysis;
%    'WienerHopf'	Wiener-Hopf equation  
%    'NBC'	Naive Bayesian Classifier [6]     
%    'aNBC'	Augmented Naive Bayesian Classifier [6]
%    'NBPW'	Naive Bayesian Parzen Window [9]     
%    'SVM11'    support vector machines, one-vs-one + voting
%               MODE.hyperparameter.c_value = 
%    'RBF'      Support Vector Machines with RBF Kernel
%               MODE.hyperparameter.c_value = 
%               MODE.hyperparameter.gamma = 
%    'LPM'      Linear Programming Machine
%               MODE.hyperparameter.c_value = 
%    'CSP'	CommonSpatialPattern is very experimental and just a hack
%		uses a smoothing window of 50 samples.
%
% 
% CC contains the model parameters of a classifier. Some time ago,     
% CC was a statistical classifier containing the mean 
% and the covariance of the data of each class (encoded in the 
%  so-called "extended covariance matrices". Nowadays, also other 
% classifiers are supported. 
%
% see also: TEST_SC, COVM
%
% References: 
% [1] R. Duda, P. Hart, and D. Stork, Pattern Classification, second ed. 
%       John Wiley & Sons, 2001. 
% [2] Peg Howland and Haesun Park,
%       Generalizing Discriminant Analysis Using the Generalized Singular Value Decomposition
%       IEEE Transactions on Pattern Analysis and Machine Intelligence, 26(8), 2004.
%       dx.doi.org/10.1109/TPAMI.2004.46
% [3] http://www-static.cc.gatech.edu/~kihwan23/face_recog_gsvd.htm
% [4] Jieping Ye, Ravi Janardan, Cheong Hee Park, Haesun Park
%       A new optimization criterion for generalized discriminant analysis on undersampled problems.
%       The Third IEEE International Conference on Data Mining, Melbourne, Florida, USA
%       November 19 - 22, 2003
% [5] J.D. Tebbens and P. Schlesinger (2006), 
%       Improving Implementation of Linear Discriminant Analysis for the Small Sample Size Problem
%	Computational Statistics & Data Analysis, vol 52(1): 423-437, 2007
%       http://www.cs.cas.cz/mweb/download/publi/JdtSchl2006.pdf
% [6] H. Zhang, The optimality of Naive Bayes, 
%	 http://www.cs.unb.ca/profs/hzhang/publications/FLAIRS04ZhangH.pdf
% [7] J.H. Friedman. Regularized discriminant analysis. 
%	Journal of the American Statistical Association, 84:165–175, 1989.
% [8] G. Fung and O.L. Mangasarian, Proximal Support Vector Machine Classifiers, KDD 2001.
%        Eds. F. Provost and R. Srikant, Proc. KDD-2001: Knowledge Discovery and Data Mining, August 26-29, 2001, San Francisco, CA.
% 	p. 77-86.
% [9] Kai Keng Ang, Zhang Yang Chin, Haihong Zhang, Cuntai Guan.
%	Filter Bank Common Spatial Pattern (FBCSP) in Brain-Computer Interface.
%	IEEE International Joint Conference on Neural Networks, 2008. IJCNN 2008. (IEEE World Congress on Computational Intelligence). 
%	1-8 June 2008 Page(s):2390 - 2397


%	$Id: train_sc.m,v 1.1 2009-01-30 06:04:49 arno Exp $
%	Copyright (C) 2005,2006,2007,2008 by Alois Schloegl <a.schloegl@ieee.org>	
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

if nargin<3, MODE = 'LDA'; end;
if ischar(MODE) 
        tmp = MODE; 
        clear MODE; 
        MODE.TYPE = tmp;
elseif ~isfield(MODE,'TYPE')
        MODE.TYPE=''; 
end;        

sz = size(D);
if sz(1)~=length(classlabel),
        error('length of data and classlabel does not fit');
end;

%CC.Labels = unique(classlabel);
CC.Labels = 1:max(classlabel);

% remove all NaN's
ix = any(isnan([D,classlabel]),2);
D(ix,:)=[];
classlabel(ix,:)=[];

sz = size(D);
if sz(1)~=length(classlabel),
        error('length of data and classlabel does not fit');
end;
if ~isfield(MODE,'hyperparameter')
        MODE.hyperparameter = [];
end


if 0, 

elseif ~isempty(strfind(lower(MODE.TYPE),'nbpw'))	
	error('NBPW not implemented yet')
	%%%% Naive Bayesian Parzen Window Classifier. 
        for k = 1:length(CC.Labels),
                [d,CC.MEAN(k,:)] = center(D(classlabel==CC.Labels(k),:),1);
                [CC.VAR(k,:),CC.N(k,:)] = sumskipnan(d.^2,1);  
                h2_opt = (4./(3*CC.N(k,:))).^(2/5).*CC.VAR(k,:);
                %%% TODO 
        end;
	
	
elseif ~isempty(strfind(lower(MODE.TYPE),'nbc'))	
	%%%% Naive Bayesian Classifier. 
	if ~isempty(strfind(lower(MODE.TYPE),'anbc'))
		%%%% Augmented Naive Bayesian classifier. 
		[CC.V,L] = eig(covm(D,'M')); 
		D = D*CC.V;
	else 
		CC.V = eye(size(D,2)); 		
	end; 
        for k = 1:length(CC.Labels),
                [d,CC.MEAN(k,:)] = center(D(classlabel==CC.Labels(k),:),1);
                [CC.VAR(k,:),CC.N(k,:)] = sumskipnan(d.^2,1);  
        end;
        CC.VAR = CC.VAR./max(CC.N-1,0); 
        CC.datatype = ['classifier:',lower(MODE.TYPE)];


elseif ~isempty(strfind(lower(MODE.TYPE),'lpm'))
        % linear programming machine 
        % CPLEX optimizer: ILOG solver, ilog cplex 6.5 reference manual http://www.ilog.com
        MODE.TYPE = 'LPM';
        if ~isfield(MODE.hyperparameter,'c_value')
                MODE.hyperparameter.c_value = 1; 
        end

        M = length(CC.Labels);
        if M==2, M=1; end;   % For a 2-class problem, only 1 Discriminant is needed 
        for k = 1:M,
                %LPM = train_LPM(D,(classlabel==CC.Labels(k)),'C',MODE.hyperparameter.c_value);
                LPM = train_LPM(D',(classlabel'==CC.Labels(k)));
                CC.weights(:,k) = [-LPM.b; LPM.w(:)];
        end;
        CC.hyperparameter.c_value = MODE.hyperparameter.c_value; 
        CC.datatype = ['classifier:',lower(MODE.TYPE)];

        
elseif ~isempty(strfind(lower(MODE.TYPE),'pls')) || ~isempty(strfind(lower(MODE.TYPE),'reg'))
        % regression analysis, kann handle sparse data, too. 
        % Q: equivalent to LDA? 
        M = length(CC.Labels); 
	%X = sparse(1:length(classlabel),classlabel,1,length(classlabel),M);
	X = sparse(length(classlabel),M);
	for k = 1:M,
		X(find(classlabel==CC.Labels(k)),k) = 2;
	end;
	CC.weights = [ones(size(D,1),1),D]\X;
	CC.weights(1,:) = CC.weights(1,:)-1;
        CC.datatype = ['classifier:statistical:',lower(MODE.TYPE)];


elseif ~isempty(strfind(MODE.TYPE,'WienerHopf'))
        % Q: equivalent to LDA, Regression? 
        M = length(CC.Labels);
        %if M==2, M==1; end;
        CC.weights = repmat(NaN,size(D,2)+1,M);
        for k = 1:M,
		ix = ~any(isnan([classlabel,D]),2);
		w  = covm(D(ix,:),'E')\covm([ones(sum(ix),1),D(ix,:)],(classlabel(ix,:)==CC.Labels(k)),'M');
                CC.weights(:,k) = w;
	end;
        CC.datatype = ['classifier:statistical:',lower(MODE.TYPE)];


elseif ~isempty(strfind(lower(MODE.TYPE),'/gsvd'))
	% [2] Peg Howland and Haesun Park, 2004. 
        %       Generalizing Discriminant Analysis Using the Generalized Singular Value Decomposition
        %       IEEE Transactions on Pattern Analysis and Machine Intelligence, 26(8), 2004.
        %       dx.doi.org/10.1109/TPAMI.2004.46
        % [3] http://www-static.cc.gatech.edu/~kihwan23/face_recog_gsvd.htm

        Hw = zeros(size(D)+[length(CC.Labels),0]); 
        Hb = [];
	m0 = mean(D); 
        K = length(CC.Labels); 
	for k = 1:K,
		ix = find(classlabel==CC.Labels(k));
		N(k) = length(ix); 
		[Hw(ix,:), mu] = center(D(ix,:));
		%Hb(k,:) = sqrt(N(k))*(mu(k,:)-m0);
		Hw(size(D,1)+k,:) = sqrt(N(k))*(mu-m0);  % Hb(k,:)
	end;
        try
                [P,R,Q] = svd(Hw,'econ');
        catch   % needed because SVD(..,'econ') not supported in Matlab 6.x
                [P,R,Q] = svd(Hw,0);
        end;
        t = rank(R);

        clear Hw Hb mu; 
        %[size(D);size(P);size(Q);size(R)]
        R = R(1:t,1:t);
        %P = P(1:size(D,1),1:t); 
        %Q = Q(1:t,:);
        [U,E,W] = svd(P(1:size(D,1),1:t),0);
        %[size(U);size(E);size(W)]
        clear U E P;  
        %[size(Q);size(R);size(W)]
        
        %G = Q(1:t,:)'*[R\W'];
        G = Q(:,1:t)*[R\W'];   % this works as well and needs only 'econ'-SVD
        %G = G(:,1:t);  % not needed 
        
        % do not use this, gives very bad results for Medline database
        %G = G(:,1:K); this seems to be a typo in [2] and [3].

        CC = train_sc(D*G,classlabel,MODE.TYPE(1:find(MODE.TYPE=='/')-1));
        CC.G = G; 
        if isfield(CC,'weights')
                CC.weights = [CC.weights(1,:); G*CC.weights(2:end,:)];
                CC.datatype = ['classifier:statistical:',lower(MODE.TYPE)];
        else
                CC.datatype = [CC.datatype,'/gsvd'];
        end;


elseif ~isempty(strfind(lower(MODE.TYPE),'sparse'))
        % [5] J.D. Tebbens and P.Schlesinger (2006), 
        %       Improving Implementation of Linear Discriminant Analysis for the Small Sample Size Problem
        %       http://www.cs.cas.cz/mweb/download/publi/JdtSchl2006.pdf

        warning('sparse LDA is sensitive to linear transformations')
        M = length(CC.Labels); 
        G  = sparse([],[],[],size(D,1),M,size(D,1));
        for k = 1:M,
                G(classlabel==CC.Labels(k),k) = 1; 
        end;
        tol  = 1e-10;

        G    = train_lda_sparse(D,G,1,tol);
        CC.datatype = 'classifier:slda';
        POS1 = find(MODE.TYPE=='/'); 
        %G = v(:,1:size(G.trafo,2)).*G.trafo; 
        %CC.weights = s * CC.weights(2:end,:) + sparse(1,1:M,CC.weights(1,:),sz(2)+1,M); 
        G  = G.trafo; 
        CC = train_sc(D*G,classlabel,MODE.TYPE(1:POS1(1)-1));
        CC.G = G; 
        if isfield(CC,'weights')
                CC.weights = [CC.weights(1,:); G*CC.weights(2:end,:)];
                CC.datatype = ['classifier:statistical:',lower(MODE.TYPE)];
        else
                CC.datatype = [CC.datatype,'/sparse'];
        end;

        
elseif ~isempty(strfind(lower(MODE.TYPE),'rbf'))
        % Martin Hieden's RBF-SVM        
        if exist('svmpredict','file')==3,
                MODE.TYPE = 'SVM:LIB:RBF';
        else
                error('No SVM training algorithm available. Install LibSVM for Matlab.\n');
        end;
        if ~isfield(MODE.hyperparameter,'gamma')
                MODE.hyperparameter.gamma = 1; 
        end
        if ~isfield(MODE.hyperparameter,'c_value')
                MODE.hyperparameter.c_value = 1; 
        end
        CC.options = sprintf('-c %g -t 2 -g %g', MODE.hyperparameter.c_value, MODE.hyperparameter.gamma);  %use RBF kernel, set C, set gamma
        CC.hyperparameter.c_value = MODE.hyperparameter.c_value; 
        CC.hyperparameter.gamma = MODE.hyperparameter.gamma; 

        % pre-whitening
        [D,r,m]=zscore(D,1); 
        CC.prewhite = sparse(2:sz(2)+1,1:sz(2),r,sz(2)+1,sz(2),2*sz(2)); 
        CC.prewhite(1,:) = -m.*r; 

        CC.model = svmtrain(classlabel, D, CC.options);    % Call the training mex File     
        CC.datatype = ['classifier:',lower(MODE.TYPE)];


elseif ~isempty(strfind(lower(MODE.TYPE),'svm11'))
        % 1-versus-1 scheme 
        if ~isfield(MODE.hyperparameter,'c_value')
                MODE.hyperparameter.c_value = 1; 
        end
        %CC = train_svm11(D,classlabel,MODE.hyperparameter.c_value);

        CC.options=sprintf('-c %g -t 0',MODE.hyperparameter.c_value);  %use linear kernel, set C
        CC.hyperparameter.c_value = MODE.hyperparameter.c_value; 

        % pre-whitening
        [D,r,m]=zscore(D,1); 
        CC.prewhite = sparse(2:sz(2)+1,1:sz(2),r,sz(2)+1,sz(2),2*sz(2)); 
        CC.prewhite(1,:) = -m.*r; 

        CC.model = svmtrain(classlabel, D, CC.options);    % Call the training mex File
        
        FUN = 'SVM:LIB:1vs1';
        CC.datatype = ['classifier:',lower(FUN)];


elseif ~isempty(strfind(lower(MODE.TYPE),'psvm'))
        if isfield(MODE.hyperparameters,'nu')
	        nu = MODE.hyperparameter.nu;
	else 
		nu = 1;          
        end;
        [m,n] = size(D); 
        CC.weights = repmat(NaN,n+1,length(CC.Labels));
        for k = 1:length(CC.Labels),
		d = sparse(1:m,1:m,(classlabel==CC.Labels(k))*2-1);
		H = d * [-ones(m,1),D];
		r = sum(H,1)';
		r = (speye(n+1)/nu + H' * H)\r; %solve (I/nu+H’*H)r=H’*e
		u = nu*(1-(H*r)); 
		CC.weights(:,k) = u'*H;
        end;
        CC.hyperparameter.nu = nu; 
        CC.datatype = ['classifier:',lower(MODE.TYPE)];
        

elseif ~isempty(strfind(lower(MODE.TYPE),'svm'))
        if ~isfield(MODE.hyperparameter,'c_value')
                MODE.hyperparameter.c_value = 1; 
        end
        if any(MODE.TYPE==':'),
                % nothing to be done
        elseif exist('svmtrain','file')==3,
                MODE.TYPE = 'SVM:LIB';
        elseif exist('svmtrain','file')==2,
                MODE.TYPE = 'SVM:bioinfo';
        elseif exist('mexSVMTrain','file')==3,
                MODE.TYPE = 'SVM:OSU';
        elseif exist('svcm_train','file')==2,
                MODE.TYPE = 'SVM:LOO';
        elseif exist('svmclass','file')==2,
                MODE.TYPE = 'SVM:KM';
        elseif exist('svc','file')==2,
                MODE.TYPE = 'SVM:Gunn';
        else
                error('No SVM training algorithm available. Install OSV-SVM, or LOO-SVM, or libSVM for Matlab.\n');
        end;

        %%CC = train_svm(D,classlabel,MODE);
        M = length(CC.Labels);
        if M==2, M=1; end;
        CC.weights = repmat(NaN, sz(2)+1, M);

        % pre-whitening
        [D,r,m]=zscore(D,1); 
        s = sparse(2:sz(2)+1,1:sz(2),r,sz(2)+1,sz(2),2*sz(2)); 
        s(1,:) = -m.*r; 
        
        for k = 1:M,
                cl = sign((classlabel~=CC.Labels(k))-.5);
                if strcmp(MODE.TYPE, 'SVM:LIB');
                        if isfield(MODE,'options')
                                CC.options = MODE.options;
                        else
                                CC.options = sprintf('-s 0 -c %f -t 0 -d 1', MODE.hyperparameter.c_value);      % C-SVC, C=1, linear kernel, degree = 1,
                        end;
                        model = svmtrain(cl, D, CC.options);    % C-SVC, C=1, linear kernel, degree = 1,
                        w = -cl(1) * model.SVs' * model.sv_coef;  %Calculate decision hyperplane weight vector
                        % ensure correct sign of weight vector and Bias according to class label
                        Bias  = -model.rho * cl(1);

                elseif strcmp(MODE.TYPE, 'SVM:bioinfo');
                        CC.SVMstruct = svmtrain(D, cl,'AUTOSCALE', 0);    % 
                        Bias = CC.SVMstruct.Bias;
                        w = CC.SVMstruct.Alpha'*CC.SVMstruct.SupportVectors;

                elseif strcmp(MODE.TYPE, 'SVM:OSU');
                        [AlphaY, SVs, Bias, Parameters, nSV, nLabel] = mexSVMTrain(D', cl', [0 1 1 1 MODE.hyperparameter.c_value]);    % Linear Kernel, C=1; degree=1, c-SVM
                        w = -SVs * AlphaY'*cl(1);  %Calculate decision hyperplane weight vector
                        % ensure correct sign of weight vector and Bias according to class label
                        Bias = -Bias * cl(1);

                elseif strcmp(MODE.TYPE, 'SVM:LOO');
                        [a, Bias, g, inds, inde, indw]  = svcm_train(D, cl, MODE.hyperparameter.c_value); % C = 1;
                        w = D(inds,:)' * (a(inds).*cl(inds)) ;

                elseif strcmp(MODE.TYPE, 'SVM:Gunn');
                        [nsv, alpha, Bias,svi]  = svc(D, cl, 1, MODE.hyperparameter.c_value); % linear kernel, C = 1;
                        w = D(svi,:)' * alpha(svi) * cl(1);
                        Bias = mean(D*w);

                elseif strcmp(MODE.TYPE, 'SVM:KM');
                        [xsup,w1,Bias,inds,timeps,alpha] = svmclass(D, cl, MODE.hyperparameter.c_value, 1, 'poly', 1); % C = 1;
                        w = -D(inds,:)' * w1;

                else
                        fprintf(2,'Error TRAIN_SVM: no SVM training algorithm available\n');
                        return;
                end

                CC.weights(1,k) = -Bias;
                CC.weights(2:end,k) = w;
        end;
        CC.weights = s * CC.weights(2:end,:) + sparse(1,1:M,CC.weights(1,:),sz(2)+1,M); % include pre-whitening transformation
        CC.hyperparameter.c_value = MODE.hyperparameter.c_value; 
        CC.datatype = ['classifier:',lower(MODE.TYPE)];


elseif ~isempty(strfind(lower(MODE.TYPE),'csp'))
        CC.datatype = ['classifier:',lower(MODE.TYPE)];
        CC.MD = repmat(NaN,[length(CC.Labels),sz(2)+[1,1]]);
        CC.NN = CC.MD;
        for k = 1:length(CC.Labels),
                [CC.MD(k,:,:),CC.NN(k,:,:)] = covm(D(classlabel==CC.Labels(k),:),'E');
        end;
        ECM = CC.MD./CC.NN;
        NC  = size(ECM);
	W   = csp(ECM,'CSP3');
	%%% ### This is a hack ###
	CC.FiltA = 50; 
	CC.FiltB = ones(CC.FiltA,1); 
	d   = filtfilt(CC.FiltB,CC.FiltA,(D*W).^2);
	CC.csp_w = W; 
	CC.CSP = train_sc(log(d),classlabel);	


else          % Linear and Quadratic statistical classifiers 
        CC.datatype = ['classifier:statistical:',lower(MODE.TYPE)];
        CC.MD = repmat(NaN,[length(CC.Labels),sz(2)+[1,1]]);
        CC.NN = CC.MD;
        for k = 1:length(CC.Labels),
                [CC.MD(k,:,:),CC.NN(k,:,:)] = covm(D(classlabel==CC.Labels(k),:),'E');
        end;

        ECM = CC.MD./CC.NN;
        NC  = size(ECM);
        if strncmpi(MODE.TYPE,'LD',2);

                %if NC(1)==2, NC(1)=1; end;                % linear two class problem needs only one discriminant
                CC.weights = repmat(NaN,NC(2),NC(1));     % memory allocation
                type = MODE.TYPE(3)-'0';

                ECM0 = squeeze(sum(ECM,1));  %decompose ECM
                [M0,sd,COV0,xc,N,R2] = decovm(ECM0);
                for k = 1:NC(1);
                        ecm = squeeze(ECM(k,:,:));
                        [M1,sd,COV1,xc,N1,R2] = decovm(ECM0-ecm);
                        [M2,sd,COV2,xc,N2,R2] = decovm(ecm);
                        switch (type)
                                case 2          % LD2
                                        cov = (COV1+COV2)/2;
                                case 4          % LD4
                                        cov = (COV1*N1+COV2*N2)/(N1+N2);
                                case 5          % LD5
                                        cov = COV2;
                                otherwise       % LD3, LDA
                                        cov = COV0/2; 
                        end
	        	if isfield(MODE.hyperparameter,'gamma')
	        		cov = cov + mean(diag(cov))*eye(size(cov))*MODE.hyperparameter.gamma;
        		end	
                        w = cov\(M2-M1)';
                        w0    = -M0*w;
                        CC.weights(:,k) = [w0; w];
                end;
                
        elseif strcmpi(MODE.TYPE,'RDA');
		if isfield(MODE,'hyperparameters') && isfield(MODE.hyperparameters,'lambda')  && isfield(MODE.hyperparameters,'gamma')
		        CC.hyperparameters = MODE.hyperparameters;
		else 
			error('QDA: hyperparamters lambda and/or gamma not defined')
		end; 	         
        else
                c  = size(ECM,2);
                ECM0 = sum(ECM,1);
                nn = ECM0(1,1,1);	% number of samples in training set for class k
                XC = squeeze(ECM0(1,:,:))/nn;		% normalize correlation matrix
                M  = XC(1,2:NC(2));		% mean
                S  = XC(2:NC(2),2:NC(2)) - M'*M;% covariance matrix
                
                [v,d]=eig(S);               
                U0 = v(diag(d)==0,:);
                CC.iS2 = U0*U0';  
                
                %M  = M/nn; S=S/(nn-1);
                ICOV0 = inv(S);
                CC.iS0 = ICOV0; 
                ICOV1 = zeros(size(S));
                for k = 1:NC(1);
                        %[M,sd,S,xc,N] = decovm(ECM{k});  %decompose ECM
                        %c  = size(ECM,2);
                        nn = ECM(k,1,1);	% number of samples in training set for class k
                        XC = squeeze(ECM(k,:,:))/nn;		% normalize correlation matrix
                        M  = XC(1,2:NC(2));		% mean
                        S  = XC(2:NC(2),2:NC(2)) - M'*M;% covariance matrix
                        %M  = M/nn; S=S/(nn-1);

                        %ICOV(1) = ICOV(1) + (XC(2:NC(2),2:NC(2)) - )/nn

                        CC.M{k}   = M;
                        CC.IR{k}  = [-M;eye(NC(2)-1)]*inv(S)*[-M',eye(NC(2)-1)];  % inverse correlation matrix extended by mean
                        CC.IR0{k} = [-M;eye(NC(2)-1)]*ICOV0*[-M',eye(NC(2)-1)];  % inverse correlation matrix extended by mean
                        d = NC(2)-1;
                        CC.logSF(k)  = log(nn) - d/2*log(2*pi) - det(S)/2;
                        CC.logSF2(k) = -2*log(nn/sum(ECM(:,1,1)));
                        CC.logSF3(k) = d*log(2*pi) + log(det(S));
                        CC.logSF4(k) = log(det(S)) + 2*log(nn);
                        CC.logSF5(k) = log(det(S));
                        CC.logSF6(k) = log(det(S)) - 2*log(nn/sum(ECM(:,1,1)));
                        CC.logSF7(k) = log(det(S)) + d*log(2*pi) - 2*log(nn/sum(ECM(:,1,1)));
                        CC.SF(k) = nn/sqrt((2*pi)^d * det(S));
                        %CC.datatype='LLBC';
                end;
        end;
end;

