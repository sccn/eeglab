function chEOG=identify_eog_channels(fn,x); 
% IDENTIFY_EOG_CHANNELS returns bipolar EOG channels for 
%  correcting of EOG artifacts using regression analysis
% 
%  EOGchan = IDENTIFY_EOG_CHANNELS(...) 
%
% EOGchan is a sparse matrix of size number_of_channels x 2. 
% The sparsity ensures that missing samples of unrelated channels 
% do not affect the data.  
%
%  [...] = IDENTIFY_EOG_CHANNELS(filename) 
%  [...] = IDENTIFY_EOG_CHANNELS(HDR) 
%	filename or HDR struct can be used
%  [...] = IDENTIFY_EOG_CHANNELS(...,'x') 
%     looks for EOG channels whos Label start with x
%
% see also: GET_REGRESS_EOG, SLOAD
%
% Reference(s):
% [1] Schlogl A, Keinrath C, Zimmermann D, Scherer R, Leeb R, Pfurtscheller G. 
%	A fully automated correction method of EOG artifacts in EEG recordings.
%	Clin Neurophysiol. 2007 Jan;118(1):98-104. Epub 2006 Nov 7.
% 	http://dx.doi.org/10.1016/j.clinph.2006.09.003
%       http://www.dpmi.tugraz.at/~schloegl/publications/schloegl2007eog.pdf

%	$Id: identify_eog_channels.m,v 1.1 2009-01-30 06:04:43 arno Exp $
%	Copyright (C) 2006,2007 by Alois Schloegl 
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

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


if ischar(fn), 
	HDR=sopen(fn); 
	HDR=sclose(HDR); 
elseif isstruct(fn)
	HDR=fn; 
end; 	

if any(any(isnan(HDR.THRESHOLD(:,1:2))))
	warning('missing threshold(s)'); 
end; 	

% graz 
g1 = strmatch('EOG-left',HDR.Label);
g2 = strmatch('EOG-central',HDR.Label);
g3 = strmatch('EOG-right',HDR.Label);
if isempty([g1,g2,g3])
	g1 = strmatch('EOG:ch01',HDR.Label);
	g2 = strmatch('EOG:ch02',HDR.Label);
	g3 = strmatch('EOG:ch03',HDR.Label);
end; 

% berlin
if nargin<2,
	v1 = strmatch('eogv1',lower(HDR.Label));
	v2 = strmatch('eogv2',lower(HDR.Label));
	v0 = strmatch('eogv',lower(HDR.Label));
	v3 = strmatch('eogvp',lower(HDR.Label));
	v4 = strmatch('eogvn',lower(HDR.Label));
	v5 = strmatch('veog',HDR.Label);

	h1 = strmatch('eogh1',lower(HDR.Label));
	h2 = strmatch('eogh2',lower(HDR.Label));
	h0 = strmatch('eogh' ,lower(HDR.Label));
	h3 = strmatch('eoghp',lower(HDR.Label));
	h4 = strmatch('eoghn',lower(HDR.Label));
else
	v1 = [];
	v2 = [];
	v0 = []; 
	v3 = strmatch('xeogvp',lower(HDR.Label));
	v4 = strmatch('xeogvn',lower(HDR.Label));

	h1 = [];
	h2 = [];
	h0 = [];
	h3 = strmatch('xeoghp',lower(HDR.Label));
	h4 = strmatch('xeoghn',lower(HDR.Label));
end;

g = [g1;g2;g3];
v = [v1,v2,v3,v4];
if isempty(v), v=v0; end; 
if isempty(v), v=v5; end; 
h = [h1,h2,h3,h4];
if isempty(h), h=h0; end; 

if length(g)==3,
	chEOG = sparse([g1,g2,g2,g3],[1,1,2,2],[1,-1,1,-1],HDR.NS,2);
elseif length(g)==2,
	chEOG = sparse([g1,g2,g3],[1,1],[1,-1],HDR.NS,1);
else 
	c = (length(v)>0);  
	sz2 = (length(v)>0) + (length(h)>0);  
	if length(v)==1, 
		chEOG = sparse(v,c,1,HDR.NS,sz2); 
	elseif length(v)==2, 
		chEOG = sparse(v,[c,c],[1,-1],HDR.NS,sz2); 
	else 
		chEOG = 0; 
	end;
	if length(h)==1, 
		chEOG = chEOG+sparse(h,1+c,1,HDR.NS,1+c); 
	elseif length(h)==2, 
		chEOG = chEOG+sparse(h,[1,1]+c,[1,-1],HDR.NS,1+c); 
	end;
end; 

if size(chEOG,2)<2, 
	warning('EOG channels are missing, or were not recognized'); 
end; 
