function [CC,KAPPA,tsd]=findclassifier(D,TRIG,cl,MODE1,t0,MODE)
% FINDCLASSIFIER
%   identifies and validates a classifier of a BCI systems [1-3]. 
%   Several evaluation criteria are obtained [4]. Several Cross-validation 
%   procedures are supported. 
%
% By default, a Trial-based Leave-One-Out-Method is used for Crossvalidation    
%	MODE.Segments = class_times; 
%	MODE.WIN   = t_ref; 
%       [CC,Q,TSD] = findclassifier(D,TRIG,Class,MODE,t_ref,TYPE);
% Also this will work but its use is discouraged (it might become obsolete). 
%       [CC,Q,TSD] = findclassifier(D,TRIG,Class,class_times,t_ref,TYPE);
% An K-fold cross-validation can be applied in this way: 
%       ng = floor([0:length(Class)-1]'/length(Class)*K);
%       [CC,Q,TSD] = findclassifier(D,TRIG,[Class,ng],...);
%
% D 	data, each row is one time point
% TRIG	trigger time points
% Class class information
% class_times	classification times, combinations of times must be in one row 
% t_ref	reference time for Class 0 (optional)
% TYPE  determines the type of classifier (see HELP TEST_SC for complete list)
%       bthe default method is 'LD3'
%
% CC 	contains the classifier and the validation results
% Q  	is a list of classification quality for each time segment (as defined by 'class_times')
% TSD 	returns the discrimination 
%
% Example: 
%  [CC,Q,TSD]=findclassifier(d,find(trig>0.5)-257,~mod(1:80,2),reshape(1:14*128,16,14*8)');
%
% see also: TRAIN_SC, TEST_SC, BCI4EVAL
%
% Reference(s): 
% [1] Schlögl A, Neuper C, Pfurtscheller G
% 	Estimating the mutual information of an EEG-based Brain-Computer-Interface
%  	Biomedizinische Technik 47(1-2): 3-8, 2002.
% [2] Schlögl A, Keinrath C, Scherer R, Pfurtscheller G,
%	Information transfer of an EEG-based Bran-computer interface.
%	Proceedings of the 1st International IEEE EMBS Conference on Neural Engineering, Capri, Italy, Mar 20-22, 2003 
% [3] Schlögl A, Lee FY, Bischof H, Pfurtscheller G
%	Characterization of Four-Class Motor Imagery EEG Data for the BCI-Competition 2005.
%	Journal of neural engineering 2 (2005) 4, S. L14-L22
% [4] Schlögl A, Kronegg J, Huggins JE, Mason SG;
%	Evaluation criteria in BCI research.
%	(Eds.) G. Dornhege, J.R. Millan, T. Hinterberger, D.J. McFarland, K.-R.Müller;
%	Towards Brain-Computer Interfacing, MIT Press, p327-342, 2007

%   $Id: findclassifier.m,v 1.1 2009-01-30 06:04:48 arno Exp $
%   Copyright (C) 1999-2006 by Alois Schloegl <a.schloegl@ieee.org>	
%   This is part of the BIOSIG-toolbox http://biosig.sf.net/


% BioSig is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
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

if ~isempty(t0)
	warning('arg5 (t_ref) should be empty. Use MODE.WIN=t_ref instead.') 
end; 

CC = []; Q = [];tsd=[];md=[];
if isstruct(MODE1)
	CC.T = MODE1; 
	T  = MODE1.Segments;
	t0 = MODE1.WIN; 
else 
	T = MODE1; 	
end; 	 

if nargin<6,
        MODE.TYPE='LD3';
elseif ischar(MODE);
        tmp = MODE; 
        clear MODE
        MODE.TYPE = tmp; 
end;
if nargin>4,
        if isempty(t0),
                t0=logical(ones(size(T,1),1));
        end;
end;
tmp=cl;tmp(isnan(tmp))=0;
if any(rem(tmp,1) & ~isnan(cl)),
        fprintf(2,'Error %s: class information is not integer\n',mfilename);
        return;
end;
if length(TRIG)~=length(cl);
        fprintf(2,'Warning FINDCLASSIFIER: number of Triggers do not match class information');
end;

if size(cl,1)~=length(TRIG);    
        fprintf(2,'Warning FINDCLASSIFIER: Classlabels must be a column vector');
        if length(TRIG)==size(cl,2),
                cl = cl';
        end;
end;

TRIG = TRIG(:);
TRIG(any(isnan(cl),2))=[];
cl(any(isnan(cl),2))=[];
if size(cl,2)>1,
        cl2 = cl(:,2);          % 2nd column contains the group definition, ( Leave-One (Group) - Out ) 
        cl  = cl(:,1); 
else
        cl2 = [1:length(cl)]';  % each trial is a group (important for cross-validation); Trial-based LOOM  
end;
[CL,i,cl] = unique(cl);
CL2 = unique(cl2);
%[CL,iCL] = sort(CL);

% add sufficient NaNs at the beginning and the end 
tmp = min(TRIG)+min(min(T))-1;
if tmp<0,
        TRIG = TRIG - tmp;
        D = [repmat(nan,[-tmp,size(D,2)]);D];
end;        
tmp = max(TRIG)+max(max(T))-size(D,1);
if tmp>0,
        D = [D;repmat(nan,[tmp,size(D,2)])];
end;        

% estimate classification result for all time segments in T - without crossvalidation 
CMX = zeros([size(T,1),length(CL)*[1,1]]);
KAPPA = repmat(NaN,size(t0));
for k = 1:size(T,1),
        if t0(k),
                c = []; d = [];
                for k1 = 1:length(CL), 
                        %t = perm(TRIG(cl==CL(k1)),T(k,:));
                        t = perm(TRIG(cl==k1),T(k,:));
                        d = [d; D(t(:),:)];
                        %c = [c; repmat(CL(k1),prod(size(t)),1)];
                        c = [c; repmat(k1,prod(size(t)),1)];
                end;
                cc{k} = train_sc(d,c,MODE);
                r     = test_sc(cc{k},d,MODE,c);
                KAPPA(k)  = r.kappa;
        end;
%        fprintf(1,'search for segment: %i-%i kappa=%4.2f\n',T(k,1),T(k,end),KAPPA(k));
end;	
[maxQ,TI] = max(KAPPA.*(t0~=0)); %d{K},
%[maxQ,TI] = max(q.*(t0~=0)); %d{K},
CC = cc{TI};
CC.KAPPA = KAPPA;
CC.TI = TI;
CC.TC = T(TI,:);
if isstruct(MODE1)
	CC.T = MODE1; 
end; 

if isnan(maxQ)
	fprintf(2,'ERROR FINDCLASSIFIER: no valid classifier available.\n'); 
	return; 
end;	

%% cross-validation with jackknife (group-based leave-one-out-method)
nc  = max(max(T))-min(min(T))+1;

M = length(CL);
tsd = repmat(nan,[nc*length(TRIG),M]);
tt  = tsd(:,1);
IX  = find(~isnan(cl(:)'));

for l = 1:length(CL2);        % XV based on "Leave-One(group)-Out-Method" 
        ix = find(cl2==CL2(l));         % identify members of l-th group 
        t  = perm(TRIG(ix), T(CC.TI,:));        % get samples of test set

%        fprintf(1,'\nX-V (%i/%i):',l, length(CL2));
        % decremental learning
        if 0, ~isempty(strfind(CC.datatype,'statistical')), 
                c  = repmat(cl(cl2==CL2(l))', size(T,2),1);     % classlabels of test set
                cc = untrain_sc(CC,c(:),D(t(:),:));             % untraining test set
                
        elseif 1, %~isempty(strfind(CC.datatype,'statistical')),
                t  = perm(TRIG(cl2~=CL2(l)), T(CC.TI,:));       % samples of training set
                c  = repmat(cl(cl2~=CL2(l))', size(T,2),1);     % classlabels of training set
                cc = train_sc(D(t(:),:),c(:),MODE);             % train classifier 

        end;
        
        t  = perm(TRIG(cl2==CL2(l)), min(min(T)):max(max(T)));  % samples of evaluation set (l-th group) 
        ix = perm((find(cl2==CL2(l))-1)*nc, 1:nc);              % save to ...     
        if any(~isnan(tt(ix(:))))
                fprintf(2,'WARNING FINDCLASSIFIER#: overlapping segments %i\n',sum(~isnan(tt(ix(:)))));
        end;
        tt(ix(:)) = t; 
        
        if ~isempty(strfind(CC.datatype,'svm:lib:1vs1')) | ~isempty(strfind(CC.datatype,'svm:lib:rbf')),
                c = repmat(cl(cl2==CL2(l))', size(t,1), 1);     % classlabels of test set
                r = test_sc(cc,D(t(:),:),MODE,c(:));            % evaluation of l-th group
        else
                r = test_sc(cc,D(t(:),:),MODE);                 % evaluation of l-th group
        end;
        tsd(ix(:),1:M) = r.output(:,1:M);                                % save results of l-th group
end; 

%CC.TSD  = bci4eval(tsd, (0:length(cl)-1)'*nc, cl, 1, nc);
CC.TSD  = bci4eval(tsd, (0:length(cl)-1)'*nc, CL(cl), 1, nc);
CC.TSD.T = CC.TSD.T - 1 + min(T(:));
