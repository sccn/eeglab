% cart2topo() - convert xyz-cartesian channel coordinates 
%               to polar topoplot coordinates. Input data
%               are points on a sphere around a given
%               [x y z] center or around (0,0,0) {default}.
%
% Usage:        >> [th r] = cart2topo(xyz);   % 3-column data
%               >> [th r] = cart2topo(x,y,z); % separate x,y,z vectors
%
%               >> [th r] = cart2topo(x,y,z,center,squeeze,optim); 
%
% Optional inputs:
%
%    center  = [ X Y Z] known center different from [0 0 0]
%
%    For the following options put center = [] if center not known:
%
%    squeeze = plotting  squeeze-in factor (0[default]->1). -1 optimize
%              the squeeze facto automatically.
%    optim   = 1: find best fitting sphere center
%
% Example:      >> [th r] = cart2topo(xyz,[1 0 4]);
%
% Notes: topoplot() does not plot channels with radius>0.5
%        Shrink radii to within this range to interpolate 
%        all channels.
%        Radii are further squeezed in of factor squeeze
%
% Important: 
%   The completed chan.locs file must have four colums
%   [channums th r chanlabels] and the chanlabels must all be 4-char 
%   strings (with . for spaces). See >> topoplot example
%
% Authors: Scott Makeig&Luca Finelli, SCCN/INC/UCSD, La Jolla, 11/1999-03/2002 
%
% See also: topo2sph(), sph2topo()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 11/1999 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: not supported by cvs2svn $

% 3-16-00 improved help message -sm
% 1-25-02 put spherror subfunction inside cart2topo -ad
% 1-25-02 help pops-up if no arguments -ad
% 01-25-02 reformated help & license -ad 
% 02-13-02 now center fitting works, need spherror outside cart2topo -lf
% 02-14-02 radii are squeezed of squeeze in to fit topoplot circle -lf
% 03-31-02 center fitting is optional
% 04-01-02 automatic squeeze calculation -ad & sm
 
function [th,r] = cart2topo(x,y,z,center,squeeze,optim)

if nargin<3
    help cart2topo
    return;
end;

if nargin<6
  optim=0;
end

if nargin<5
  squeeze=0;
end

if exist('squeeze') & (squeeze>1)
  fprintf('Warning: Squeeze must be less than 1.\n');
  return
end

if size(x,2)==3 % separate 3-column data
   if nargin>1
     center = y;
   end
   z = x(:,3);
   y = x(:,2);
   x = x(:,1);
end

if exist('center') & (~isempty(center) & size(center,2) ~= 3)
  fprintf('Warning: Center must be [x y z].\n');
  return
end

if size(x,2)>1 % convert to columns
   x = x'
end

if nargin>1
  if size(y,2)>1
    x = x'
 end
end

if nargin>2
  if size(z,2)>1
   x = x';
  end
end

options = [];
if ~exist('center') | isempty(center)

% Find center
% ----------------------------------------------
  fprintf('Optimising center position...\n');
  center = fminsearch('spherror',[0 0 0],options,x,y,z);
  fprintf('Best center is [%g,%g,%g].\n',center(1),center(2),center(3));
end
x = x - center(1);  % center the data
y = y - center(2);
z = z - center(3);
radius = (sqrt(x.^2+y.^2+z.^2));   % assume xyz values are on a sphere
wobble = std(radius);              % test if xyz values are on a sphere
fprintf('Radius values: %g (mean) +/- %g (std)\n',mean(radius),wobble);

if  wobble/mean(radius) > 0.01 & optim==1
  kk=0;
  while wobble/mean(radius) > 0.01 & kk<5
    newcenter = fminsearch('spherror',center,options,x,y,z);
    nx = x - newcenter(1);  % re-center the data
    ny = y - newcenter(2);
    nz = z - newcenter(3);
    nradius = (sqrt(nx.^2+ny.^2+nz.^2));   % assume xyz values are on a sphere
    newobble = std(nradius);   
    if newobble<wobble
      center=newcenter;
      fprintf('Wobble too high (%3.2g%%)! Re-centering data on (%g,%g,%g)\n',...
                100*wobble/mean(radius),newcenter(1),newcenter(2),newcenter(3))
      x = nx;  % re-center the data
      y = ny;
      z = nz;
      radius=nradius;
      wobble=newobble;
      kk=kk+1;
    else
      kk=5;
    end
  end
   fprintf('Wobble (%3.2g%%) after centering data on (%g,%g,%g)\n',...
              100*wobble/mean(radius),center(1),center(2), ...
	   center(3))
%else
%  fprintf('Wobble (%3.2g%%) after centering data on (%g,%g,%g)\n',...
%              100*wobble/mean(radius),center(1),center(2),center(3))
end

x = x./radius; % make radius 1
y = y./radius;
z = z./radius;

r = x; th=x;

for n=1:size(x,1)
  if x(n)==0 & y(n)==0
    r(n) = 0;
  else
    r(n)  = pi/2-atan(z(n)./sqrt(x(n).^2+y(n).^2));
  end
end

r = r/pi;        % scale r to max 0.500
if squeeze < 0
    squeeze = 1 - 0.5/max(r); %(2*max(r)-1)/(2*rmax);
    fprintf('Electrodes will be squeezed to show all using squeezing act of %2.3g\n', squeeze);
end;
r = r-squeeze*r; % squeeze electrodes in squeeze*100% to have all inside

for n=1:size(x,1)
  if abs(y(n))<1e-6
    if x(n)>0
      th(n) = -90;
    else % x(n) <= 0
      th(n) = 90;
    end
  else
    th(n) = atan(x(n)./y(n))*180/pi+90;
    if y(n)<0
      th(n) = th(n)+180;
    end
  end
  if th(n)>180 
     th(n) = th(n)-360;
  end
  if th(n)<-180 
     th(n) = th(n)+360;
  end
end
