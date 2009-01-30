function [R]=hist2res(H,fun)
% Evaluates Histogram data
% [R]=hist2res(H)
%
% [y]=hist2res(H,fun)
%	estimates fun-statistic
%
% fun	'mean'	mean
%	'std'	standard deviation
%	'var'	variance
%	'sem'	standard error of the mean
%	'rms'	root mean square
%	'meansq' mean of squares
%	'sum'	sum
%	'sumsq'	sum of squares
%	'CM#'	central moment of order #
%	'skewness' skewness 
%	'kurtosis' excess coefficient (Fisher kurtosis)
%
% see also: NaN/statistic
%
% REFERENCES:
% [1] C.L. Nikias and A.P. Petropulu "Higher-Order Spectra Analysis" Prentice Hall, 1993.
% [2] C.E. Shannon and W. Weaver "The mathematical theory of communication" University of Illinois Press, Urbana 1949 (reprint 1963).
% [3] http://www.itl.nist.gov/
% [4] http://mathworld.wolfram.com/

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

%	$Id: hist2res.m,v 1.1 2009-01-30 06:04:43 arno Exp $
%	Copyright (c) 1996-2002,2006 by Alois Schloegl <a.schloegl@ieee.org>
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/


if strcmp(H.datatype,'HISTOGRAM')

elseif strcmp(H.datatype,'qc:histo')
	HDR = H; 
	if isfield(H,'THRESHOLD'),
		TH  = H.THRESHOLD;
	else
		TH = repmat([-inf,inf],HDR.NS,1); 
	end;
	HIS = H.HIS; 

	% remove overflowing samples
	HIS.N = sumskipnan(HIS.H); 
	for k = 1:size(HIS.H,2);
		t = HIS.X(:,min(k,size(HIS.X,2))); 
		HIS.H(xor(t<=min(TH(k,:)), t>=max(TH(k,:))),k) = 0; 
	end; 
	Nnew = sumskipnan(HIS.H); 
	R.ratio_lost = 1-Nnew./HIS.N;
	HIS.N = Nnew; 
	  
	% scale into physical values
	if H.FLAG.UCAL,
		%t = HIS.X;
		%for k=1:length(HDR.InChanSelect),
		%	HIS.X(:,k) = t(:,min(size(t,2),k))*HDR.Calib(k+1,k)+HDR.Calib(1,k);
		%end;
		HIS.X = [ones(size(HIS.X,1),1),repmat(HIS.X,1,size(HIS.H,2)./size(HIS.X,2))]*H.Calib;
	end; 	
	H = HIS; 
else
        fprintf(2,'ERROR: arg1 is not a histogram\n');
        return;
end;
if nargin<2, fun=[]; end;

global FLAG_implicit_unbiased_estimation; 
%%% check whether FLAG was already defined 
if exist('FLAG_implicit_unbiased_estimation')~=1,
	FLAG_implicit_unbiased_estimation=[];
end;
%%% set DEFAULT value of FLAG
if isempty(FLAG_implicit_unbiased_estimation),
	FLAG_implicit_unbiased_estimation=logical(1);
end;

sz 	= size(H.H)./size(H.X);
R.N 	= sumskipnan(H.H,1);
R.SUM 	= sumskipnan(H.H.*repmat(H.X,sz),1);
R.SSQ 	= sumskipnan(H.H.*repmat(H.X.*H.X,sz),1);
%R.S3P 	= sumskipnan(H.H.*repmat(H.X.^3,sz),1);	% sum of 3rd power
R.S4P 	= sumskipnan(H.H.*repmat(H.X.^4,sz),1);	% sum of 4th power
%R.S5P 	= sumskipnan(H.H.*repmat(H.X.^5,sz),1);	% sum of 5th power

R.MEAN	= R.SUM./R.N;
R.MSQ   = R.SSQ./R.N;
R.RMS	= sqrt(R.MSQ);
R.SSQ0  = R.SSQ-R.SUM.*R.MEAN;		% sum square of mean removed

if FLAG_implicit_unbiased_estimation,
    n1 	= max(R.N-1,0);			% in case of n=0 and n=1, the (biased) variance, STD and STE are INF
else
    n1	= R.N;
end;

R.VAR  	= R.SSQ0./n1;	     		% variance (unbiased) 
R.STD  	= sqrt(R.VAR);		     	% standard deviation
R.SEM  	= sqrt(R.SSQ0./(R.N.*n1)); 	% standard error of the mean
R.SEV	= sqrt(n1.*(n1.*R.S4P./R.N+(R.N.^2-2*R.N+3).*(R.SSQ./R.N).^2)./(R.N.^3)); % standard error of the variance
R.Coefficient_of_variation = R.STD./R.MEAN;

R.CM2	= R.SSQ0./n1;
x       = repmat(H.X,sz) - repmat(R.MEAN,size(H.X,1),1);
R.CM3 	= sumskipnan(H.H.*(x.^3),1)./n1;
R.CM4 	= sumskipnan(H.H.*(x.^4),1)./n1;
%R.CM5 	= sumskipnan(H.H.*(x.^5),1)./n1;

R.SKEWNESS = R.CM3./(R.STD.^3);
R.KURTOSIS = R.CM4./(R.VAR.^2)-3;
R.MAD = sumskipnan(H.H.*abs(x),1)./R.N; % mean absolute deviation

H.PDF = H.H./H.N(ones(size(H.H,1),1),:);
status=warning('off'); 
R.ENTROPY = -sumskipnan(H.PDF.*log2(H.PDF),1);
warning(status); 
R.QUANT = repmat(min(diff(H.X,[],1)),1,size(H.H,2)/size(H.X,2));
R.MAX = max(H.X); 
R.MIN = min(H.X); 
R.RANGE = R.MAX-R.MIN;

if ~isempty(fun),
        fun=upper(fun);
        if strncmp(fun,'CM',2) 
                oo = str2double(fun(3:length(fun)));
                R = sumskipnan(H.PDF.*(x.^oo),1);
    	else	            
		R = getfield(R,fun);
	end;
end;

