function [crit,t] = criteria2005IIIb(X)
% criteria2005IIIb evalates the criterium according the 
% BCI competition 2005, dataset IIIb. 
%
%  X = bci4eval(...); 
%  [crit] = criteria2005IIIb(X)
%
% see also: SUMSKIPNAN, PLOTA, BCI3EVAL, BCI4EVAL


%    $Revision: 1.1 $
%    $Id: criteria2005IIIb.m,v 1.1 2009-01-30 06:04:50 arno Exp $
%    Copyright (C) 2004,2005 by Alois Schloegl <a.schloegl@ieee.org>	
%    This is part of the BIOSIG-toolbox http://biosig.sf.net/

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


crit = NaN; 
%if isfield(X,'datatype') & strcmp(X.datatype,'BCI_TSD7')
if isfield(X,'T') & isfield(X,'I');
        [crit,t] = max(X.I./(X.T-3).*(X.T>3.5));
        t = X.T(t);
end;
