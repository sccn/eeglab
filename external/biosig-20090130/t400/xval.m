function [R,CC]=xval(D,classlabel,MODE,arg4)
% XVAL is used for crossvalidation 
%
%  [R,CC] = xval(D,classlabel)
%  .. = xval(D,classlabel,CLASSIFIER)
%  .. = xval(D,classlabel,CLASSIFIER,type)
%  .. = xval(D,[classlabel,NG],CLASSIFIER)
%
% Input:
%    D 	data features (one feature per column, one sample per row)
%    classlabel  classlabel (same number of rows)
%    CLASSIFIER can be any classifier supported by train_sc (default='LDA')
%    NG: used to define the type of cross-valdiation
% 	Leave-One-Out-Method (LOOM): NG = [1:length(classlabel)]' (default)
% 	Leave-K-Out-Method: NG = ceil([1:length(classlabel)]'/K)
%	K-fold XV:  NG = ceil([1:length(classlabel)]'*K/length(classlabel))
%	group-wise XV (if samples are not indepentent) can be also defined here
%	samples from the same group (dependent samples) get the same identifier
%	samples from different groups get different classifiers
%    TYPE  defines the type of cross-validation procedure if NG is not specified 
%	'LOOM'  leave-one-out-method
%       k	k-fold crossvalidation
%
% OUTPUT: 
%    R contains the resulting performance metric
%    CC contains the parameters of the classifier  
%
% see also: TRAIN_SC, TEST_SC, CLASSIFY
%
% References: 
% [1] R. Duda, P. Hart, and D. Stork, Pattern Classification, second ed. 
%       John Wiley & Sons, 2001. 

%	$Id: xval.m,v 1.1 2009-01-30 06:04:49 arno Exp $
%	Copyright (C) 2008 by Alois Schloegl <a.schloegl@ieee.org>	
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

if (nargin<3) || isempty(MODE), 
	MODE = 'LDA'; 
end;
if ischar(MODE) 
        tmp = MODE; 
        clear MODE; 
        MODE.TYPE = tmp;
elseif ~isfield(MODE,'TYPE')
        MODE.TYPE=''; 
end;        

sz = size(D);
if sz(1)~=size(classlabel,1),
        error('length of data and classlabel does not fit');
end;

if size(classlabel,2)>1,
	%% group-wise classvalidation
	[tmp,tmp1,NG] = unique(classlabel(:,2));
elseif nargin<4
	%% LOOM 
	NG = [1:sz(1)]';	
elseif isnumeric(arg4)
	if isscalar(arg4)  
	% K-fold XV
		NG = ceil([1:length(classlabel)]'*arg4/length(classlabel));
	elseif length(arg4)==2,
		NG = ceil([1:length(classlabel)]'*arg4(1)/length(classlabel));
	end; 	
	
elseif strcmpi(arg4,'LOOM')
	NG = [1:sz(1)]';	
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


for k = 1:sz(1),
 	ix = find(NG~=k);
	CC = train_sc(D(ix,:),classlabel(ix,1),MODE);
 	ix = find(NG==k);
	r  = test_sc(CC,D(ix,:));
	cl(ix,1) = r.classlabel;
end; 
R = kappa(cl,classlabel(:,1));
R.ERR = 1-R.ACC; 

if nargout>1,
	% final classifier 
	CC = train_sc(D,classlabel(ix),MODE);
end; 
