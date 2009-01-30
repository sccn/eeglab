function topo2(V,HDR,maplimits)
% TOPO2 makes a 2-D topographic map
%
%  topo2(value,HDR [,limits]); 
%  
%
%  

%	$Id: topo2.m,v 1.1 2009-01-30 06:04:51 arno Exp $
%	Copyright (C) 2007,2008 by Alois Schloegl <a.schloegl@ieee.org>	
%       This is part of the BIOSIG-toolbox http://biosig.sf.net/
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



Label = []; 
if isstruct(HDR)
	if isfield(HDR,'ELEC')
		XYZ = HDR.ELEC.XYZ; 
		Label = HDR.Label; 
	elseif isfield(HDR,'XYZ')
		XYZ = HDR.XYZ; 
	end; 
elseif isnumeric(HDR),
	XYZ = HDR; 	
end; 	
size(XYZ),size(V),

ix = find(all(~isnan(XYZ),2) & ~isnan(V(:))); 

R   = sqrt(sum(XYZ(ix,:).^2,2));
th  = atan(XYZ(ix,2)./XYZ(ix,1)) + pi*(XYZ(ix,1)<0);
rd  = acos(XYZ(ix,3)./R);
th(isnan(th))=0;
V = V(ix); 
Label = Label(ix); 

if nargin<3,
	maplimits=[-1,1]*max(abs(V)); 
%elseif isnumeric(maplimits)
elseif ischar(maplimits)
 	if strcmp(maplimits,'absmax'); 
 		maplimits=[-1,1]*max(abs(V)); 
	elseif strcmp(maplimits,'maxmin'); 
		maplimits=[min(V),max(V)];
	end; 	 
end; 
 
 
r   = max(rd)*1.05;
x   = -rd.*sin(th-pi/2);
y   = rd.*cos(th-pi/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% interpolation         %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GRID_SCALE = 67; 
xi = linspace(-r,+r,GRID_SCALE); 
yi = linspace(-r,+r,GRID_SCALE); 

[Xi,Yi,Zi] = griddata(x,y,V,xi',yi,'invdist'); % interpolate data

Zi(Xi.^2 + Yi.^2 > r.^2) = NaN; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% make topographic plot %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%surf(Xi,Yi,Zi); 
caxis(maplimits);
surface(Xi,Yi,zeros(size(Zi)),Zi,'EdgeColor','none','FaceColor','flat'); 
hold on; 
%contour(Xi,Yi,Zi,6,'k');
t = 0:.01:2*pi;
plot(x,y,'k.',r*sin([-.1,0,.1]'),r*[cos([-.1,0,.1]')+[0,.1,0]'],'k-',r*sin(t),r*cos(t),'k-'); 
%for k=1:length(V),text(x(k),y(k),Label{k}),end;
set(gca,'visible','off','xlim',[-1,1]*r,'ylim',[-1,1.1]*r); 
hold off;

 
