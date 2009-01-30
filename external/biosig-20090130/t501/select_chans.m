function [locs_n,inds]=select_chans(locs,nin)
% usage [locs_n,inds]=select_chans(locs,nin)
% selects a sparse subset of locations, such 
% the selected locations are as equally distributed as 
% possible
% 
% input:
% locs: locations of sensors (2D or 3D), each row is the location of one sensor
% nin:  number of desired channels
%      
% output:
% locs_n: selcted locations
% inds: indices of selected locations

%  $Id: select_chans.m,v 1.1 2009-01-30 06:04:51 arno Exp $ 
%  Copyright (C) 2003,2004 Guido Nolte
%  Adapted by Alois Schloegl 
%  This function is part of the BioSig project
%  http://biosig.sf.net/	

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the  License, or (at your option) any later version.
% 


[n,m]=size(locs);

inds=[1];

for k=1:nin-1
    
    %loct=locs(inds);
    
    p=0;
    dis=zeros(n,1);
    for i=1:n
        disloc=zeros(k,1);
        for j=1:k
        disloc(j)=norm(locs(i,:)-locs(inds(j),:));
        end
        dis(i)=min(disloc);
    end
    
    
    [dmax,imax]=max(dis);
    
    inds=[inds;imax];
end
locs_n=locs(inds,:);
return;
    