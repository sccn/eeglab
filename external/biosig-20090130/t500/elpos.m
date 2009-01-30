function [YX,code]=elpos(Label);
% ELPOS provides electrode positions in 2-D according to [1]
%   
% [YX,code]=elpos(Label);
%
% see also: ELPOS3
%
% Reference(s): 
% [1] FEF: 
%   File Exchange Format for Vital Signs - Normative Annex A (Normative). 
%   The Medical Data Information Base (MDIB), Nomenclature, Data Dictionary and Codes
%   Table A.6.3: Nomenclature and Codes for Electrode Sites 
%   for Electroencephalography according to the International 10-20 system.
%   CEN/TC251/PT-40 (2001)

%	$Id: elpos.m,v 1.1 2009-01-30 06:04:50 arno Exp $
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
XYZ = BIOSIG_GLOBAL.XYZ;
Phi   = atan(XYZ(:,1)./XYZ(:,2))*180/pi;
Theta = acos(XYZ(:,3)./sqrt(sum(XYZ.^2,2)));
xy  = Theta.*exp(i*angle(XYZ*[1;i;0]))*180/pi; 

for k=1:nr;
for l=1:length(BIOSIG_GLOBAL.Label),
	if strcmp(upper(deblank(Label(k,:))),upper(BIOSIG_GLOBAL.Label{l}))
		code(k)=l;
break;
	end;	
end;
end;

K=code(code>0)';

T=0:.001:2*pi;
R=180/pi*2;
plot(real(xy),imag(xy),'x',real(xy(K)),imag(xy(K)),'ro',R*sin(T),R*cos(T),'b-',-R+R/10*sin(-T/2),10*cos(-T/2),'b-',10*sin(T/2)+R,10*cos(T/2),'b-',[-10 0 10],[R R+10 R],'b-');
for k=1:size(XYZ,1), 
	if all(~isnan(XYZ(k,:))),
		text(real(xy(k)),imag(xy(k)),BIOSIG_GLOBAL.Label{k});
	end;
end;
set(gca,'xtick',0,'ytick',0,'xticklabel','','yticklabel','')


