% del2map() - compute the discrete laplacian of an EEG distribution.
%
% Usage:
% >> [ laplac ] = del2map( map, filename, draw );
%
% Inputs:
%    map        - level of activity (size: nbelectrodes * nbChannel)
%    filename	- filename (.loc file) countaining the coordinates
%                 of the electrodes, or array countaining complex positions 		 
%    draw       - integer, if not nul draw the gradient (default:0)
%
% Output:
%    laplac     - laplacian map 
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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
% Revision 1.3  2003/10/30 22:02:16  arno
% allowing to read any channel location file
%
% Revision 1.2  2003/09/04 23:20:38  arno
% draw
%
% Revision 1.1  2002/04/05 17:39:45  jorn
% Initial revision
%
% 01-25-02 reformated help & license -ad 

function [ laplac, sumLaplac2D ] = del2map( map, filename, draw )

if nargin < 2
	help del2map;
	return;
end;

% process several maps
if size(map,2) > 1
    if size(map,1) > 1
        for index = 1:size(map,2)
            laplac(:,index) = del2map( map(:,index), filename);
        end;
        return;
    else
        map = map';
    end;
end;

MAXCHANS = size(map,1);
GRID_SCALE = 2*MAXCHANS+5;

% Read the channel file
% ---------------------
if isstr( filename ) | isstruct( filename )
    [tmp lb Th Rd] = readlocs(filename);
	Th = pi/180*Th;                               % convert degrees to rads
	[x,y] = pol2cart(Th,Rd);
else
	x = real(filename);
	y = imag(filename);		
    if exist('draw') == 1 & draw ~= 0
        line( [(x-0.01)' (x+0.01)']', [(y-0.01)' (y+0.01)']');
        line( [(x+0.01)' (x-0.01)']', [(y-0.01)' (y+0.01)']');
    end;
end;	

% locates nearest position of electrod in the grid 
% ------------------------------------------------
xi = linspace(-0.5,0.5,GRID_SCALE);   % x-axis description (row vector)
yi = linspace(-0.5,0.5,GRID_SCALE);   % y-axis description (row vector)
for i=1:MAXCHANS
   [useless_var horizidx(i)] = min(abs(x(i) - xi));    % find pointers to electrode
   [useless_var vertidx(i)] = min(abs(y(i) - yi));     % positions in Zi
end;
   
% Compute gradient
% ----------------
sumLaplac2D = zeros(GRID_SCALE, GRID_SCALE);
for i=1:size(map,2) 
   	[Xi,Yi,Zi] = griddata(y,x,map(:,i),yi',xi, 'invdist');   % interpolate data
   	laplac2D = del2(Zi);
	positions = horizidx + (vertidx-1)*GRID_SCALE;

	laplac(:,i) = laplac2D(positions(:));
	sumLaplac2D = sumLaplac2D + laplac2D;

	% Draw gradient
	% -------------
	if exist('draw') == 1 & draw ~= 0
		if size(map,2) > 1
			subplot(ceil(sqrt(size(map,2))), ceil(sqrt(size(map,2))), i);
		end;
		plot(y, x, 'x', 'Color', 'black', 'markersize', 5); hold on
		contour(Xi, Yi, laplac2D); hold off;
		%line( [(x-0.01)' (x+0.01)']', [(y-0.01)' (y+0.01)']', 'Color','black');
		%line( [(x+0.01)' (x-0.01)']', [(y-0.01)' (y+0.01)']', 'Color','black');
		title( int2str(i) );
	end;
end;                                                             %
plot(y, x, 'x', 'Color', 'black', 'markersize', 5); hold on
contour(Xi, Yi, sumLaplac2D); hold off;

return;

