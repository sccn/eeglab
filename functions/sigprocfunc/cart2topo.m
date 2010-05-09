% cart2topo() - convert xyz-cartesian channel coordinates 
%               to polar topoplot() coordinates. Input data
%               are points on a sphere centered at (0,0,0)
%               or at optional input 'center'. This function
%               is now DEPRECATED! See Important warning below.
%
% Usage:  >> [th r] = cart2topo(xyz);   % 3-column [x y z] position data
%         >> [th r] = cart2topo(x,y,z); % separate x,y,z vectors
%         >> [th r x y z] = cart2topo(xyz,'key', 'val', ...); 
%         >> [th r x y z] = cart2topo(x,y,z,'key', 'val', ...); 
%
% Optional inputs:
%    'center'  = [X Y Z] use a known sphere center {Default: [0 0 0]}
%    'squeeze' = plotting squeeze factor (0[default]->1). -1 -> optimize
%                the squeeze factor automatically so that maximum radius is 0.5. 
%    'optim'   = [0|1] find the best-fitting sphere center minimizing the
%                standard deviation of the radius values.
%    'gui'     = ['on'|'off'] pops up a gui for optional arguments. 
%                Default is off.
%
% Example: >> [th r] = cart2topo(xyz,[1 0 4]);
%
% Notes: topoplot() does not plot channels with radius>0.5
%        Shrink radii to within this range to plot all channels.
%        [x y z] are returned after the optimization of the center
%        and optionally squeezing r towards it by factor 'squeeze'
%
% Important: 
%   DEPRECATED: cart2topo() should NOT be used if elevation angle is less than 0 
%   (for electrodes below zero plane) since then it returns INNACURATE results. 
%   SUBSTITUTE: Use cart2topo = cart2sph() -> sph2topo().
%
% Authors: Scott Makeig, Luca Finelli & Arnaud Delorme SCCN/INC/UCSD,
%          La Jolla, 11/1999-03/2002 
%
% See also: topo2sph(), sph2topo(), chancenter()

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

% 3-16-00 improved help message -sm
% 1-25-02 put spherror subfunction inside cart2topo -ad
% 1-25-02 help pops-up if no arguments -ad
% 01-25-02 reformated help & license -ad 
% 02-13-02 now center fitting works, need spherror outside cart2topo -lf
% 02-14-02 radii are squeezed of squeeze in to fit topoplot circle -lf
% 03-31-02 center fitting is optional
% 04-01-02 automatic squeeze calculation -ad & sm
 
function [th,r,xx,yy,zz] = cart2topo(x,varargin)

if nargin<1
    help cart2topo
    return;
end;
if nargin >= 2
	if ~ischar(varargin{1})
		y = varargin{1};
		z = varargin{2};
		varargin = varargin(3:end);
	end;
end;
if exist('y') ~= 1 
	if size(x,2)==3 % separate 3-column data
		z = x(:,3);
		y = x(:,2);
		x = x(:,1);
	else
		error('Insufficient data in first argument');
	end
end;

g = [];
if ~isempty(varargin)
    try, g = struct(varargin{:}); 
    catch, error('Argument error in the {''param'', value} sequence'); end; 
end;

try, g.optim;      catch, g.optim = 0; end;
try, g.squeeze;    catch, g.squeeze = 0; end;
try, g.center;     catch, g.center = [0 0 0]; end;
try, g.gui;        catch, g.gui = 'off'; end;

if g.squeeze>1
  fprintf('Warning: Squeeze must be less than 1.\n');
  return
end

if ~isempty(g.center) & size(g.center,2) ~= 3
  fprintf('Warning: Center must be [x y z].\n');
  return
end

% convert to columns
x = -x(:); % minus os for consistency between measures
y = -y(:);
z = z(:);

if any(z < 0)
    disp('WARNING: some electrodes lie below the z=0 plane, result may be innacurate')
    disp('         Instead use cart2sph() then sph2topo().')
end;
if strcmp(g.gui, 'on')
	[x y z newcenter] = chancenter(x, y, z, [], 1);
else 
	if g.optim == 1
		[x y z newcenter] = chancenter(x, y, z, []);
	else
		[x y z newcenter] = chancenter(x, y, z, g.center);
	end;
end;
radius = (sqrt(x.^2+y.^2+z.^2));   % assume xyz values are on a sphere
xx=-x; yy=-y; zz=z;
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
if g.squeeze < 0
    g.squeeze = 1 - 0.5/max(r); %(2*max(r)-1)/(2*rmax);
    fprintf('Electrodes will be squeezed together by %2.3g to show all\n', g.squeeze);
end;
r = r-g.squeeze*r; % squeeze electrodes in squeeze*100% to have all inside

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
