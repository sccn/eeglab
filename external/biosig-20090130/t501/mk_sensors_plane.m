function loc_phys=mk_sensors_plane(sensors,pars);
% MK_SENSORS_PLANE
%  ### docu incomplete ###

%  $Id: mk_sensors_plane.m,v 1.1 2009-01-30 06:04:51 arno Exp $ 
%  Copyright (C) 2003,2004 Guido Nolte
%  Adapted by Alois Schloegl 
%  This function is part of the BioSig project
%  http://biosig.sf.net/	

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the  License, or (at your option) any later version.
% 

[ns_ori,ndim]=size(sensors);


if nargin>1
  if isfield(pars,'rot');
      rot=pars.rot;
  else
      rot=0;
  end
  if isfield(pars,'nin');
      nin=pars.nin;
  else
      nin=ns_ori;
  end
  if isfield(pars,'indices');
      indices=pars.indices;
  else
      indices=[1:ns_ori]';
  end
  if isfield(pars,'circle_shift');
      shift_cont=pars.circle_shift;
  else
      shift_cont=1;
  end

else
    rot=0;
    nin=ns_ori;
    indices=[1:ns_ori]';
    shift_cont=1;
end



if nin<ns_ori;
    [sensors_n,inds]=select_chans(sensors,nin);
else
    sensors_n=sensors;
end


if ndim==3
   s2d=sensor3d2sensor2d(sensors);
else
   s2d=sensors;
end


if nin<ns_ori;
    s2d_n=s2d(inds,:);
else
    sensors_n=sensors;
    s2d_n=s2d;
end

[ns,ndum]=size(sensors_n);

phi=rot*pi/180;
   
s2d=([[cos(phi),-sin(phi)];[sin(phi),cos(phi)]]*s2d')';
s2d_n=([[cos(phi),-sin(phi)];[sin(phi),cos(phi)]]*s2d_n')';
 
     
if shift_cont>0  
   loc_phys_sparse=locphys2locphys([[1:ns]',s2d_n]);
else
    loc_phys_sparse=[[1:ns]',s2d_n,s2d_n];
   figure;
   for i=1:ns;
       text(s2d_n(i,1),s2d_n(i,2),num2str(i));
   end
   axis([-.6 .6 -.6 .6]);  

end

   if nin<ns_ori
       loc_phys=[-(1:ns_ori)',s2d,s2d];
       loc_phys(inds,:)=[inds,loc_phys_sparse(:,2:5)];
   else
       loc_phys=loc_phys_sparse;
   end
   
loc_phys(:,1)=indices.*sign(loc_phys(:,1));

return;