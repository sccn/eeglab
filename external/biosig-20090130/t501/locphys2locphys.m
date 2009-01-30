function [loc_phys,names]=locphys2locphys(loc_phys_in);
% Purpose: create 'physical' locations from neoroscan-data
%
% usage: [loc_phys,names]=loc2locphys(fn,ind_phys,loctype)
%
% output: loc_phys  nX5 matrix; each row contains information 
%                   about a physical electrode in the form 
%                   [channelnumber x_coord y_coord x'_coord y'_coord]
%                   the primed coordinates  differ from the 
%                   non-primed coordinates if loctype is set to 1. 
%                   Then the primed coordinates are slightly shifted to 
%                    avoid overlapping subplots in plots with subplots 
%                    placed on electrode locations. 
%           names  list of channel names
%
% input: fn    filename of neurocan-file with electrode locations and names
%        ind_phys  row-vector of indices which are regarded as 'physical'
%                   i.e. electrodes on the scalp
%                  (e.g. if we have 64 channels and channels 31,32,61,62,63,64 
%                   do not refer to electrodes on the scalp, then 
%                   ind_phys=[1:30,33:60];
%        loctype    if set to 1:    calculate modified coordinates for plotting purposes 
%                   any other value: set the shifted locations to the original 
%
%  $Id: locphys2locphys.m,v 1.1 2009-01-30 06:04:51 arno Exp $ 
%  Copyright (C) 2003,2004 Guido Nolte
%  Adapted by Alois Schloegl 
%  This function is part of the BioSig project
%  http://biosig.sf.net/	

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the  License, or (at your option) any later version.
% 

ind_phys=loc_phys_in(:,1);
[nphys,m]=size(ind_phys);

xx=loc_phys_in(:,2);
yy=loc_phys_in(:,3);


xx(abs(xx)<.01)=0;
yy(abs(yy)<.01)=0;


figure
plot(xx,yy,'.');
hold on;
for i=1:nphys,
  text(xx(i),yy(i),num2str(i));
end

xp=xx;
yp=yy;
step=.0001;
mag=200;
%d=.13;
d=1.2*sqrt(pi/4/nphys);
  [xp,yp]=plotcoordc(xx,yy,d,mag,step);
  
  figure
  plot(xp,yp,'.');
  hold on;
  for i=1:nphys
    text(xp(i),yp(i),num2str(i));
  end



loc_phys=[ind_phys,xx,yy,xp,yp];


return; 
function [xp,yp]=plotcoordc(x,y,d,mag,step);

[n,m]=size(x);

dx=zeros(n,1);
dy=zeros(n,1);




for k=1:150
  [x,y]=onestep(x,y,d,mag,step);
  [minmin,minmax]=mindis(x,y);
%  display([k,minmin,minmax]);
  

  if round(k/1)*1==k
    %figure;
    %plot(x,y,'.');
  end
end
xp=x;yp=y;

return;

function [x,y]=onestep(x,y,d,mag,step);

 xn=x;yn=y;
 [n,m]=size(x);
 for i=1:n
     fallx=0;fally=0;
     r=norm([x(i),y(i)])+sqrt(eps);
     for j=1:n
         if j~=i
             [fx,fy]=force(x(i),x(j),y(i),y(j),d,mag);
             fallx=fallx+fx;fally=fally+fy;
         end
         
     end
     fallx=fallx-x(i)/(r);fally=fally-y(i)/(r);
     
     if abs(x(i))>.01
       xn(i)=x(i)+step*fallx;
     end
     if abs(y(i))>.01
       yn(i)=y(i)+step*fally;
     end
     
 end
             
 x=xn;y=yn;
          
 return;
 
 function [fx,fy]=force(x1,x2,y1,y2,d,mag);

 dis=max(abs(x1-x2),abs(y1-y2));
 dis=norm([x1-x2,y1-y2]);
 diffvec=[x2-x1;y2-y1];diffvec=diffvec/norm(diffvec);
 fx=0;fy=0;
 if dis<d
     fx=-mag*diffvec(1);
     fy=-mag*diffvec(2);
 end
     fx=-mag*exp(-dis^4/d^4)*diffvec(1);
     fy=-mag*exp(-dis^4/d^4)*diffvec(2);
     
 
 return; 
     
 
 
 
 
 
 
function [nover,deltax,deltay]=overlap(x,y,i,d,step);

[n,m]=size(x);
nover=0;deltax=0;deltay=0;

for j=1:n
if j~=i 

   dlx=abs(x(j)-x(i));
   dly=abs(y(j)-y(i));
   if  dlx<d && dly<d
       nover=1;
     if dlx<dly 
        deltax=deltax-step*(x(j)-x(i));
     else 
        deltay=deltay-step*(y(j)-y(i));
     end   
   else
       deltax=-step*sign(x(i));
       deltay=-step*sign(x(i));
        
   end
      
end
end


return; 

function [minmin,minmax]=mindis(x,y);

[n,m]=size(x);

minall=zeros(n,1);



for i=1:n
    md=1.e8;
    for j=1:n
        if j~=i
          %dis=max( abs(x(i)-x(j)),abs(y(i)-y(j)));
          dis=norm([x(i)-x(j),y(i)-y(j)]);
          if dis<md
             md=dis;
          end
        end     
    end
    minall(i)=md;
end



minmin=min(minall);
[minmax,imax]=max(minall);

return;

