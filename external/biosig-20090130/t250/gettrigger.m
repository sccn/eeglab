function TRIG = gettrigger(s,TH,rfp);
% GETTRIGGER identifies trigger points 
%
% TRIG = gettrigger( s [,TH [,rfp]]); % trigger at ascending edge
% TRIG = gettrigger(-s [,TH [,rfp]]); % trigger at decending edge
%
% input : s   	signal
%         TH	Threshold; default: (max+min)/2
%	  rfp   refractory period (default=0)	
% output: TRIG	TRIGGER time points
%
% see also: TRIGG

% 	$Id: gettrigger.m,v 1.1 2009-01-30 06:04:43 arno Exp $
%	Copyright (C) 2002-2003,2008 by Alois Schloegl <a.schloegl@ieee.org>		

% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 2 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.


if nargin<2,
	TH = (max(s)+min(s))/2;
end;
if nargin<3,
	rfp = 0; 
end; 	

TRIG = find(diff(sign(s-TH))>0)+1;
% perform check of trigger points

if (rfp<=0), return; end; 

% refractory period 
k0=1;
k1=2;
while k1<length(TRIG);
	T0 = TRIG(k0);
	T1 = TRIG(k1);
	if (T1-T0)<rfp,
		TRIG(k1)=NaN;
	else
		k0 = k1; 
	end
	k1 = k1+1;
end;	
TRIG=TRIG(~isnan(TRIG));


