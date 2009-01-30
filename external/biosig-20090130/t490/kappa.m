function [kap,se,H,z,p0,SA,R]=kappa(d,c,arg3);
% KAPPA estimates Cohen's kappa coefficient
%   and related statistics 
%
% [...] = kappa(d1,d2);
%	NaN's are handled as missing values and are ignored
% [...] = kappa(d1,d2,'notIgnoreNAN');
%	NaN's are handled as just another Label.
% [kap,sd,H,z,ACC,sACC,MI] = kappa(...);
% X = kappa(...);
%
% d1    data of scorer 1 
% d2    data of scorer 2 
%
% kap	Cohen's kappa coefficient point
% se	standard error of the kappa estimate
% H	Concordance matrix, i.e. confusion matrix
% z	z-score
% ACC	overall agreement (accuracy) 
% sACC	specific accuracy 
% MI 	Mutual information or transfer information (in [bits])
% X 	is a struct containing all the fields above
%
% Reference(s):
% [1] Cohen, J. (1960). A coefficient of agreement for nominal scales. Educational and Psychological Measurement, 20, 37-46.
% [2] J Bortz, GA Lienert (1998) Kurzgefasste Statistik f|r die klassische Forschung, Springer Berlin - Heidelberg. 
%        Kapitel 6: Uebereinstimmungsmasze fuer subjektive Merkmalsurteile. p. 265-270.
% [3] http://www.cmis.csiro.au/Fiona.Evans/personal/msc/html/chapter3.html
% [4] Kraemer, H. C. (1982). Kappa coefficient. In S. Kotz and N. L. Johnson (Eds.), 
%        Encyclopedia of Statistical Sciences. New York: John Wiley & Sons.
% [5] http://ourworld.compuserve.com/homepages/jsuebersax/kappa.htm
%
%  

%	$Id: kappa.m,v 1.1 2009-01-30 06:04:50 arno Exp $
%	Copyright (c) 1997-2006,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/
%
%    BioSig is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BioSig is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with BioSig.  If not, see <http://www.gnu.org/licenses/>.


mode.ignoreNAN = 1; 
kk = [];
if nargin>2
	if ischar(arg3)
		if strcmpi(arg3,'notIgnoreNAN')
			mode.ignoreNAN = 0; 
		 end
	else 
		kk = arg3; 
	end
end; 		 


if nargin>1,
	d = d(:);
	c = c(:);
	
	if mode.ignoreNAN,
		if any(isnan([d;c]))
               		fprintf(2,'Warning KAPPA: some elements are NaN. These are handled as missing values and are ignored.\n');
               		fprintf(2,'If NaN should be handled as just another label, use kappa(..,''notIgnoreNaN'').\n');
			ix = find(~isnan(c) & ~isnan(d));
			d = d(ix); c=c(ix);
		end;
		[X.Label,i,j]   = unique([d;c]);
	else 
		tmp = [d;c]; 		 
		tmp(isnan(tmp))=max(tmp)+1;
		[X.Label,i,j]   = unique(tmp);
		X.Label(end)=NaN; 
	end;
	[X.Label,i,j]   = unique([d;c]);
	c = j(1+numel(d):end); 
	d = j(1:numel(d)); 
	
    	N  = length(d);
    	ku = max([d;c]); % upper range
    	kl = min([d;c]); % lower range
	
    	if isempty(kk),
            	kk = length(X.Label);  	% maximum element
    	else
            	if kk<ku;  	% maximum element
                    	fprintf(2,'Error KAPPA: some element is larger than arg3(%i)\n',kk);
            	end;
    	end;
    
    	if 0,
        	h = histo([d+c*kk; kk*kk+1; 1]); 
        	H = reshape(h(1:length(h)-1));
        	H(1,1) = H(1,1)-1;
    	else
		if 1;   % exist('OCTAVE_VERSION')>=5;
	        	H = zeros(kk);
    			for k = 1:N, 
    				if ~isnan(d(k)) & ~isnan(c(k)),
		    			H(d(k),c(k)) = H(d(k),c(k)) + 1;
		    		end;	
        		end;
		else
			H = full(sparse(d(1:N),c(1:N),1,kk,kk));
    		end;
	end;
else
	X.Label = 1:min(size(d));
    	H = d(X.Label,X.Label);
    	% if size(H,1)==size(H,2);	
	N = sum(sum(H));
    	% end;
end;
warning('off');
p0  = sum(diag(H))/N;  %accuracy of observed agreement, overall agreement 
%OA = sum(diag(H))/N);

p_i = sum(H); %sum(H,1);
pi_ = sum(H'); %sum(H,2)';

SA  = 2*diag(H)'./(p_i+pi_); % specific agreement 

pe  = (p_i*pi_')/(N*N);  % estimate of change agreement

px  = sum(p_i.*pi_.*(p_i+pi_))/(N*N*N);

%standard error 
kap = (p0-pe)/(1-pe);
sd  = sqrt((pe+pe*pe-px)/(N*(1-pe*pe)));

%standard error 
se  = sqrt((p0+pe*pe-px)/N)/(1-pe);
z = kap/se;

if ((1 < nargout) & (nargout<7)) return; end; 

% Nykopp's entropy
pwi = sum(H,2)/N;                       % p(x_i)
pwj = sum(H,1)/N;                       % p(y_j)
pji = H./repmat(sum(H,2),1,size(H,2));  % p(y_j | x_i) 
R   = - sumskipnan(pwj.*log2(pwj)) + sumskipnan(pwi'*(pji.*log2(pji)));

if (nargout>1) return; end; 

X.kappa = kap; 
X.kappa_se = se; 
X.H = H;
X.z = z; 
X.ACC = p0; 
X.sACC = SA;
X.MI = R;

if length(H)==2,
	X.FNR = H(2,1)/sum(H(2,:));
	X.FPR = H(1,2)/sum(H(1,:));
	X.TPR = H(2,2)/sum(H(2,:));
end;

kap = X;  

