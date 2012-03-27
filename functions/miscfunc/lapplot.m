% lapplot() -  Compute the discrete laplacian of EEG scalp distribution(s)
%                
% Usage:
%   >> laplace = lapplot(map,eloc_file,draw)
% 
% Inputs:
%    map        - Activity levels, size (nelectrodes,nmaps)
%    eloc_file	- Electrode location filename (.loc file) 
%                 For format, see  >> topoplot example 
%    draw       - If defined, draw the map(s) {default: no}
%
% Output:
%    laplace    - Laplacian map, size (nelectrodes,nmaps)
%
% Note: uses del2()
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 1998 
%
% See also: topoplot(), gradplot()

% Copyright (C) Scott Makeig, SCCN/INC/UCSD, La Jolla, 1998 
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

% 01-25-02 reformated help & license, added links -ad 

function [laplac] = lapplot(map,filename,draw)

if nargin < 2
	help lapplot;
	return;
end;

MAXCHANS = size(map,1);
GRID_SCALE = 2*MAXCHANS+5;
MAX_RADIUS = 0.5;

% ---------------------
% Read the channel file
% ---------------------
if isstr( filename )
	fid = fopen(filename); 
	locations = fscanf(fid,'%d %f %f %s',[7 MAXCHANS]);
	fclose(fid);
	locations = locations';
	Th = pi/180*locations(:,2);   % convert degrees to rads
	Rd = locations(:,3);
	ii = find(Rd <= MAX_RADIUS); % interpolate on-scalp channels only
	Th = Th(ii);
	Rd = Rd(ii);
	[x,y] = pol2cart(Th,Rd);
else
	x = real(filename);
	y = imag(filename);
end;	

% ---------------------------------------------------
% Locate nearest position of an electrode in the grid 
% ---------------------------------------------------
xi = linspace(-0.5,0.5,GRID_SCALE);   % x-axis description (row vector)
yi = linspace(-0.5,0.5,GRID_SCALE);   % y-axis description (row vector)
for i=1:MAXCHANS
   [useless_var horizidx(i)] = min(abs(y(i) - xi));    % find pointers to electrode
   [useless_var vertidx(i)] = min(abs(x(i) - yi));     % positions in Zi
end;
   
% -----------------
% Compute laplacian
% -----------------
for i=1:size(map,2) 
   	[Xi,Yi,Zi] = griddata(y,x,map(:,i),yi',xi, 'v4');   % interpolate data

   	laplac2D = del2(Zi);
	positions = horizidx + (vertidx-1)*GRID_SCALE;
	laplac(:,i) = laplac2D(positions(:));

	% ------------------
	% Draw laplacian map
	% ------------------
	if exist('draw');
        mask = (sqrt(Xi.^2+Yi.^2) <= MAX_RADIUS);
        laplac2D(find(mask==0)) = NaN;

		subplot(ceil(sqrt(size(map,2))), ceil(sqrt(size(map,2))), i);
		contour(laplac2D); 
		title( int2str(i) );

% %%% Draw Head %%%%
ax = axis; 
width = ax(2)-ax(1);
axis([ax(1)-width/3 ax(2)+width/3 ax(3)-width/3 ax(4)+width/3])
steps = 0:2*pi/100:2*pi;
basex = .18*MAX_RADIUS;  
tip = MAX_RADIUS*1.15; 
base = MAX_RADIUS-.004;
EarX = [.497 .510 .518 .5299 .5419 .54 .547 .532 .510 .489];
EarY = [.0555 .0775 .0783 .0746 .0555 -.0055 -.0932 -.1313 -.1384 -.1199];

HCOLOR = 'k';
HLINEWIDTH = 1.8;

% Plot Head, Ears, Nose
hold on
plot(1+width/2+cos(steps).*MAX_RADIUS*width,...
     1+width/2+sin(steps).*MAX_RADIUS*width,...
    'color',HCOLOR,'Linestyle','-','LineWidth',HLINEWIDTH); % head

plot(1+width/2+[.18*MAX_RADIUS*width;0;-.18*MAX_RADIUS*width],...
     1+width/2+[base;tip;base]*width,...
    'Color',HCOLOR,'LineWidth',HLINEWIDTH);                 % nose
   
plot(1+width/2+EarX*width,...
     1+width/2+EarY*width,...
           'color',HCOLOR,'LineWidth',HLINEWIDTH)           % l ear
plot(1+width/2-EarX*width,...
     1+width/2+EarY*width,...
           'color',HCOLOR,'LineWidth',HLINEWIDTH)           % r ear

hold off
axis off

	end;
end;                   

return;
