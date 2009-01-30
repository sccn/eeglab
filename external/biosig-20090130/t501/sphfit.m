function [center,radius]=sphfit(vc)
% SPHFIT fits a sphere to a set of surface points
%
% input: 
% vc   nx3 matrix, where each row represents the location
%      of one surface point. vc can have more than 3 columns 
%      (e.g. orientations) - then only the first 3 columns are used
%
% center  1x3 vector denoting the center
% radius  scalar denoting the radius 

%  $Id: sphfit.m,v 1.1 2009-01-30 06:04:51 arno Exp $ 
%  Copyright (C) 2003,2004 Guido Nolte
%  Adapted by Alois Schloegl 
%  This function is part of the BioSig project
%  http://biosig.sf.net/	

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the  License, or (at your option) any later version.
% 

vc=vc(:,1:3);
[nvc,ndum]=size(vc);

center_0=mean(vc);
vcx=vc-repmat(center_0,nvc,1);
radius_0=mean(sqrt(vcx(:,1).^2+vcx(:,2).^2+vcx(:,3).^2));

alpha=1;
err_0=costfun(vc,center_0,radius_0);

for k=1:5;
    
    [center_new,radius_new]=lm1step(vc,center_0,radius_0,alpha);
    
    err_new=costfun(vc,center_new,radius_new);
    %disp([k,err_0,err_new,center_new,radius_new]);
    
    if err_new<err_0;
        center_0=center_new;
        radius_0=radius_new;
        err_0=err_new;
        alpha=alpha/5;
    else
        alpha=alpha*5;
    end
    
    radius=radius_0;
    center=center_0;
    
    
end

return;

function err=costfun(vc,center,radius);
[nvc,ndum]=size(vc);


 vcx=vc-repmat(center,nvc,1);
 
 err=sqrt(sqrt(mean( (vcx(:,1).^2+vcx(:,2).^2+vcx(:,3).^2-radius^2).^2)));
 
 return;
 
 function  [center_new,radius_new]=lm1step(vc,center,radius,alpha);
 
 [nvc,ndum]=size(vc);
 vcx=vc-repmat(center,nvc,1);
 f=vcx(:,1).^2+vcx(:,2).^2+vcx(:,3).^2-radius^2;
 
 L=2*[vcx,repmat(radius,nvc,1)];
 
 par_new=inv(L'*L+alpha*eye(4))*L'*f;
 
 center_new=center+par_new(1:3)';
 radius_new=radius+par_new(4);
 
 return;
 
 
 
 