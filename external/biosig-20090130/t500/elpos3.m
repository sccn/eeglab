function [xyz,code]=elpos3(Label);
% ELPOS3 provides electrode positions in 3-D according to FEF 
%   
% [YX,code]=elpos3(Label);
%
% see also: ELPOS
%
% Reference(s): 
% [1] FEF: 
%   File Exchange Format for Vital Signs - Normative Annex A (Normative). 
%   The Medical Data Information Base (MDIB), Nomenclature, Data Dictionary and Codes
%   Table A.6.3: Nomenclature and Codes for Electrode Sites 
%   for Electroencephalography according to the International 10-20 system.
%   CEN/TC251/PT-40 (2001)
 

%	$Id: elpos3.m,v 1.1 2009-01-30 06:04:51 arno Exp $
%	Copyright (c) 1997,1998,2004,2007 by Alois Schloegl
%	a.schloegl@ieee.org	
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


tmp = leadidcodexyz('Fp1');   % make sure BIOSIG_GLOBAL.XYZ is loaded 
global BIOSIG_GLOBAL


if nargin<1,
        Label='';
end
[nr,nc]=size(Label);
code=zeros(nr,1);
xyz=ones(nr,3);

XYZ = BIOSIG_GLOBAL.XYZ;

for k=1:nr,
for l=1:length(BIOSIG_GLOBAL.Label),%size(Electrode.Theta,2), 
	if strcmp(upper(deblank(Label(k,:))),upper(BIOSIG_GLOBAL.Label{l}))
		code(k)=l;
		xyz = BIOSIG_GLOBAL.XYZ(l,:);
break;
	end;
end;
end;


K=code(code>0)';
plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'x',XYZ(K,1),XYZ(K,2),XYZ(K,3),'ro');
for k=1:size(XYZ,1),
	if all(~isnan(XYZ(k,:))),
		text(XYZ(k,1),XYZ(k,2),XYZ(k,3),BIOSIG_GLOBAL.Label{k});
	end;
end; 
set(gca,'xtick',0,'ytick',0)


