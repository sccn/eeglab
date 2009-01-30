function [x,sz] = trigg(s,TRIG,pre,post,gap)
% TRIGG cuts continous sequence into segments.
% Missing values (in case s is to short) are substituted by NaN's.  
%
% [X,sz] = trigg(s, TRIG, PRE, PST [, GAP])
%
% INPUT:
%  S	is the continous data sequence (1 channel per column)
%  TRIG	defines the trigger points
%  PRE 	offset of the start of each segment (relative to trigger) 
%  PST 	offset of the end of each segment (relative to trigger) 
%  GAP	number of NaN's to separate trials (default=0)
%  	TRIG, pre, post and gap are counted in samples
%
% OUTPUT:
%  X	is a matrix of size [sz(1), sz(2)*sz(3)]
%  sz	[size(s,2), post-pre+1+gap, length(TRIG)]   
%	sz(1) is the number of channels NS 
%	sz(2) is the number of samples per trial 
%	sz(3) is the number of trials i.e. length(TRIG)
%
% X3D = reshape(X,sz) returns a 3-dimensional matrix 
% X2D = reshape(X(K,:),sz(2:3)) returns channel K in a 2-dimensional matrix 
%
% see also: GETTRIGGER

% 	$Id: trigg.m,v 1.1 2009-01-30 06:04:43 arno Exp $
%	Copyright (c) 1999-2005 by Alois Schloegl <a.schloegl@ieee.org>
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

if nargin<5,
	gap = 0;
end;

post=round(post);

[nr,nc] = size(s);

% include leading nan's
off  = min(min([TRIG(:);+Inf])+pre-1,0);
% include following nan's
off2 = max(max([TRIG(:);-Inf])+post-length(s),0);
if ((off~=0) || (off2~=0))
	s    = [repmat(nan,-off,nc);s;repmat(nan,off2,nc)];
	TRIG = TRIG-off;
end; 

% devide into segments
N   = post-pre+1+gap;
sz  = [nc, post-pre+1+gap, length(TRIG)];
x   = repmat(NaN, [sz(1), sz(2)*sz(3)]);   
for m = 1:length(TRIG),
	x(:,m*N + (1-N:-gap)) = s(TRIG(m)+(pre:post)',:).';
end;

