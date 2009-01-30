function INI = tvaar(signal,arg2,arg3,arg4)
% TVAAR wrapper around adaptive autoregressive estimator. 
%       X = tvaar(signal,p,UC)
%       X = tvaar(signal,X)
%  
% INPUT:
%   X.MOP=[d,p,q]
%       d = 1: with mean term; d=0: without mean term
%       p:    order of AutoRegressive part 
%       q:    order of Moving Average
%   X.UC  update coefficient
%   X.Mode = [amode,vmode]: choose estimation algorithm    
%   X.Z0  covariance of state estimates
%   X.z0  initial state vector 
%   X.W0  covariance of system noise 
%   X.V0  variance of observation noise
%
% OUTPUT:
%   X.Z0  average covariance of state estimates
%   X.z0  average state vector 
%   X.W0  covariance of system noise 
%   X.V0  variance of prediction error
%   X.AAR estimated AAR parameters
%   X.E   prediction error, residuum, 
%   X.PE  time-varying variance of residual process
% 
% REFERENCES: 
% [1] Schlögl A.(2000)
%   The electroencephalogram and the adaptive autoregressive model: theory and applications
%   Shaker Verlag, Aachen, Germany,(ISBN3-8265-7640-3). 
% [2] Schlögl A, Lee FY, Bischof H, Pfurtscheller G
%   Characterization of Four-Class Motor Imagery EEG Data for the BCI-Competition 2005.
%   Journal of neural engineering 2 (2005) 4, S. L14-L22

%  	$Id: tvaar.m,v 1.1 2009-01-30 06:04:44 arno Exp $
%	Copyright (c) 2003-2004,2007 by Alois Schloegl <a.schloegl@ieee.org>	
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


if isstruct(arg2) 
        INI = arg2;
else
        INI.MOP = arg2;
        INI.UC = arg3;
end;

if ~isfield(INI,'z0');
        INI.z0 = zeros(1,sum(INI.MOP));
        INI.Mode = [17,4];
end;
if ~isfield(INI,'Z0');
        INI.Z0 = eye(sum(INI.MOP));
end;
if ~isfield(INI,'V0');
        INI.V0 = 1;
end;
if ~isfield(INI,'W0');
        [z,E,ESU,REV,V,Z,SPUR] = amarma(signal, INI.Mode, INI.MOP, INI.UC, INI.z0, INI.Z0, INI.V0);        
else
        INI.Mode = [99,5];
        [z,E,ESU,REV,V,Z,SPUR] = amarma(signal, INI.Mode, INI.MOP, INI.UC, INI.z0, INI.Z0, INI.V0, INI.W0);
end;

ix  = ~isnan(E);
PE  = E; 
tmp = filter(INI.UC,[1,INI.UC-1],[mean(V);E(ix).*E(ix)]);
PE(find(ix)) = tmp(2:end);

tmp    = z(~isnan(E),:);
INI.C0 = covm(tmp,'E');
%INI.Z0 = covm(tmp,'D');
%INI.z0 = mean(tmp);
[INI.z0,sd,INI.Z0] = decovm(INI.C0);
INI.W0 = covm(diff(tmp,[],1));
INI.V0 = meansq(E);

INI.AAR= z;
INI.E  = E;
INI.PE = PE;

INI.S  = signal; 

INI.datatype = 'AMARMA';
