function showfield_general(z,loc,pars); 
% usage showfield_general(z,loc); 
% displays fields/potentials specified z in channels at locations 
% specified in locs as a contour-plot
%
% pars is optional
%   pars.scale sets the scale of the color map. Either a 1x2 vector
%         corresponding to minimum and a maximum or just a number (say x) 
%          then the scale is from [-x x]. The default is 
%          scale=[ -max(abs(z)) max(abs(z)) ]

%  $Id: showfield_general.m,v 1.1 2009-01-30 06:04:51 arno Exp $ 
%  Copyright (C) 2003,2004 Guido Nolte
%  Adapted by Alois Schloegl 
%  This function is part of the BioSig project
%  http://biosig.sf.net/	

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the  License, or (at your option) any later version.
% 


figkont=1;


if nargin>2
    if isfield(pars,'scale')
        scal=pars.scale;
        if length(scal)==1
            scal=[-scal scal];
        end
    else
        scal=[-max(max(abs(z))),max(max(abs(z)))];
    end

    if isfield(pars,'figkont')
       figkont=pars.figkont;
    end
  else
     scal=[-max(max(abs(z))),max(max(abs(z)))];
  
end



[n,m]=size(loc);
if m==2;
  x=loc(:,1);
  y=loc(:,2);
else;
  x=loc(:,2);
  y=loc(:,3);
end


xlin = linspace(1.4*min(x),1.4*max(x),250);
ylin = linspace(1.4*min(y),1.4*max(y),250);
[X,Y] = meshgrid(xlin,ylin);
Z = griddata(x,y,z,X,Y,'invdist');
%Z = griddata(x,y,z,X,Y,'nearest');


  % Take data within head
  rmax=1.02*max(sqrt(x.^2+y.^2));
  mask = (sqrt(X.^2+Y.^2) <= rmax);
  ii = find(mask == 0);
  Z(ii) = NaN;
  
  
surface(X,Y,zeros(size(Z)),Z,'edgecolor','none');shading interp;
%caxis([ - max(abs(z)) max(abs(z))]);
caxis(scal);
hold on;
plot(x,y,'.k','markersize',2);


%meanx=mean(loc(:,2)*.85+.45);
%meany=mean(loc(:,3)*.85+.45);
scalx=1;
drawhead(0,.0,rmax,scalx);
set(gca,'visible','off');

%axis([-1.2*rmax 1.2*rmax -1.0*rmax 1.4*rmax]);
if figkont==1
axis([-1.2*rmax 1.2*rmax -1.0*rmax 1.4*rmax]);
h=colorbar;set(h,'fontweight','bold')
else
axis([-1.4*rmax 1.4*rmax -1.0*rmax 1.4*rmax]);
end

%plot(.985*rmax*sin((0:1000)/1000*2*pi), .985*rmax*cos((0:1000)/1000*2*pi),'linewidth',2,'color','k'); 
return; 
function drawhead(x,y,size,scalx);

cirx=(x+scalx*size*cos((1:1000)*2*pi/1000) )';ciry=(y+size*sin((1:1000)*2*pi/1000))';

plot(cirx,ciry,'k','linewidth',1);
hold on;

ndiff=20;
plot( [x  cirx(250-ndiff) ],[y+1.1*size ciry(250-ndiff)],'k','linewidth',1);
plot( [x  cirx(250+ndiff) ],[y+1.1*size ciry(250+ndiff)],'k','linewidth',1);


return;
