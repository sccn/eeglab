function [LIM1,LIM2,LIM3,h0] = hist2limits(H,TH)
% HIST2LIMITS returns the threshold for detecting saturation artifacts. 
%
% Saturation thresholds can be obtained from the histogram [1]. This 
% routine tries to obtain the saturation threshold in an automated way. 
% 
% The routine was tested with the histograms of 528 recordings with 
% three respiratory channels each. 
%
%
% Reference(s): 
% [1] A. Schlögl, B. Kemp, T. Penzel, D. Kunz, S.-L. Himanen,A. Värri, G. Dorffner, G. Pfurtscheller.
%     Quality Control of polysomnographic Sleep Data by Histogram and Entropy Analysis.
%     Clin. Neurophysiol. 1999, Dec; 110(12): 2165 - 2170.


%       $Revision: 1.1 $
%	$Id: hist2limits.m,v 1.1 2009-01-30 06:04:43 arno Exp $
%	Copyright (C) 1999-2003 by Alois Schloegl <a.schloegl@ieee.org>

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




if isnumeric(H)
        H = histo2(H);
end;

for k = 1:size(H.H,2),
        h = H.H(:,k);
	x = H.X(:,min(k,size(H.X,2)));

        tmp = find(h>0);
        h(tmp([1,length(tmp)])) = 0;	% remove max and min values; (necessary for H recordings and some B) 
        N = sum(h);
	
%        Lim = H.X(find(h>0),1);	% calculate limit values of remaining Histogram
        Lim = x(find(h>0));		% calculate limit values of remaining Histogram
        Lim = [max(Lim),min(Lim)];
        LIM1(:,k) = sort(mean(Lim)+([1;-1]*TH/2)*abs(diff(Lim))); 	% take range between 10% and 90% of total range. 

	mu = x'*h/N;
	sd = sqrt(((x-mu)'.^2*h)/N);

	LIM2(:,k) = sort(sd*[1;-1]*norminv(2/N)+mu);
	
	h0 = sum(h)*normpdf(x,mu,sd);

	r0 = (h./max(h0,1e-2));
	[s,ix]  = sort(-r0);
	tmp = ix(x(ix)>mu);
	LIM3(:,k) = [NaN;NaN];
	if 1;r0(tmp(1))>1e2;
		LIM3(1,k) = x(tmp(1));
	end;
	tmp = ix(x(ix)<mu);
	if 1;r0(tmp(1))>1e3;
		LIM3(2,k) = x(tmp(1));
	end;
	
        LIM3(:,k) = sort(mean(LIM3(:,k))+([1;-1]*TH/2)*abs(diff(LIM3(:,k)))); 	% take range between 10% and 90% of total range. 


	LIM0(:,k) = [max(LIM1(1,k),LIM2(1,k));min(LIM1(2,k),LIM2(2,k))];
	LIM = LIM1;
end;


