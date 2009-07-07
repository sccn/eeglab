function [y]=logit(x)
% LOGIT
%
%    y = logit (x)
%       y = log(x/(1-x)) = log(x) - log(1-x); 
%
% features:
% - can deal with NaN's (missing values)
% - accepts dimension argument like in Matlab in Octave, too. 
% - compatible to Matlab and Octave 
%
% see also: PROBIT 
% 
% Reference: 
%    http://en.wikipedia.org/wiki/Logit

%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

%	$Revision: 1.1 $
%	$Id: logit.m,v 1.1 2009-07-07 02:23:49 arno Exp $
%	Copyright (C) 2000-2003,2006 by Alois Schloegl <a.schloegl@ieee.org>	



y = log(x) - log(1-x); 

