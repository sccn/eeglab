% chancenter() - recenter cartesian X,Y,Z channel coordinates
%
% Usage:  >> [x y z newcenter] = chancenter(x,y,z,center); 
%
% Optional inputs:
%    x,y,z     = 3D coordintates of the channels
%    center    = [X Y Z] known center different from [0 0 0]
%                [] will optimize the center location according
%                to the best sphere. Default is [0 0 0].
%
% Note: 4th input gui is obsolete. Use pop_chancenter instead.
%
% Authors: Arnaud Delorme, Luca Finelli & Scott Makeig SCCN/INC/UCSD,
%          La Jolla, 11/1999-03/2002 
%
% See also: spherror(), cart2topo()

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
% 1-25-02 put spherror subfunction inside chancenter -ad
% 1-25-02 help pops-up if no arguments -ad
% 01-25-02 reformated help & license -ad 
% 02-13-02 now center fitting works, need spherror outside chancenter -lf
% 02-14-02 radii are squeezed of squeeze in to fit topoplot circle -lf
% 03-31-02 center fitting is optional
% 04-01-02 automatic squeeze calculation -ad & sm
 
function [ x, y, z, newcenter, optim] = chancenter( x, y, z, center, gui)

optim = 0;

if nargin<4
    help chancenter
    return;
end;

if nargin > 4 & gui
    error('Chancenter: 4th input'' gui'' is obsolete. Use pop_chancenter instead');
else 
	if isempty(center)
		optim = 1;
		center = [0 0 0];
	end;
end;

options = {'MaxFunEvals',1000*length(center)};
x = x - center(1);  % center the data
y = y - center(2);
z = z - center(3);
radius = (sqrt(x.^2+y.^2+z.^2));   % assume xyz values are on a sphere
if ~isempty(radius)
     wobble = std(radius);              % test if xyz values are on a sphere
else wobble = [];
end;
fprintf('Radius values: %g (mean) +/- %g (std)\n',mean(radius),wobble);
newcenter = center;

if  wobble/mean(radius) > 0.01 & optim==1
	% Find center
	% ----------------------------------------------
	fprintf('Optimizing center position...\n');
	kk=0;
	while wobble/mean(radius) > 0.01 & kk<5
        try
    		newcenter = fminsearch('spherror',center,struct(options{:}),x,y,z);
        catch
            try
        		newcenter = fminsearch('spherror',center,options,x,y,z);
            catch
            	newcenter = fminsearch('spherror',center,[], [], x,y,z);
            end
        end
		nx = x - newcenter(1);  % re-center the data
		ny = y - newcenter(2);
		nz = z - newcenter(3);
		nradius = (sqrt(nx.^2+ny.^2+nz.^2));   % assume xyz values are on a sphere
		newobble = std(nradius);   
		if newobble<wobble
			center=newcenter;
			fprintf('Wobble too strong (%3.2g%%)! Re-centering data on (%g,%g,%g)\n',...
					100*wobble/mean(radius),newcenter(1),newcenter(2),newcenter(3))
			x = nx;  % re-center the data
			y = ny;
			z = nz;
			radius=nradius;
			wobble=newobble;
			kk=kk+1;
		else
			newcenter = center;
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
