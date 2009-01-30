   function s2d=sensor3d2sensor2d(sensors,lr,ud);   
% SENSOR3D2SENSOR2D
%  --### docu missing ###	
   
%  $Id: sensor3d2sensor2d.m,v 1.1 2009-01-30 06:04:51 arno Exp $ 
%  Copyright (C) 2003,2004 Guido Nolte
%  Adapted by Alois Schloegl 
%  This function is part of the BioSig project
%  http://biosig.sf.net/	

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the  License, or (at your option) any later version.
% 


   [ns,ndum]=size(sensors);
   [center,radius]=sphfit(sensors);
   s=sensors-repmat(center,ns,1);
   cms=mean(sensors);
   zz=cms-center;
   zz=zz/norm(zz);
   snorms=sqrt(sum((s').^2)');
   s=s./repmat(snorms,1,3);
   thetas=acos( s*zz');
   ex=[1;0;0];ey=[0;1;0];
   yy=cross(zz',ex);yy=yy/norm(yy);
   xx=cross(yy,zz');xx=xx/norm(xx);
   xx=xx';
   yy=yy';
   
   xxx=s*xx';
   yyy=s*yy';
   
   phis=angle(xxx+sqrt(-1)*yyy);
   
   locs=[cos(phis).*thetas,sin(phis).*thetas];
   
   if nargin>1 
       if lr==1
           locs=[-locs(:,1),locs(:,2)];
       end
   end
   
   if nargin>2 
       if ud==1
           locs=[locs(:,1),-locs(:,2)];
       end
   end

   locs=locs/max(max(abs(locs)))/2;
   
   
   
   
   s2d=locs;
   
   return; 
   
   
   
   
   