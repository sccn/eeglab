% del2map() - compute the discrete laplacian of an EEG distribution.
%
% Usage:
% >> [ laplac ] = del2map( map, filename, draw );
%
% Inputs:
%    map        - level of activity (size: nbChannel)
%    filename	- filename (.loc file) countaining the coordinates
%                 of the electrodes, or array countaining complex positions 		 
%    draw       - integer, if not nul draw the gradient (default:0)
%
% Output:
%    laplac     - laplacian map. If the input values are in microV and the
%                 the sensors placement are in mm, the output values are
%                 returned in microV/mm^2. In order to use current density 
%                 units like milliamps/mm2, you would  need to know skin 
%                 conductance information, which as far as we know is not
%                 really known with enough accuracy to be worthwhile. 
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%         Thanks Ramesh Srinivasan and Tom Campbell for the discussion
%         on laplacian output units.

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [ laplac, sumLaplac2D ] = del2map( map, filename, draw )

if nargin < 2
	help del2map;
	return;
end

% process several maps
if size(map,2) > 1
    if size(map,1) > 1
        for index = 1:size(map,2)
            laplac(:,index) = del2map( map(:,index), filename);
        end
        return;
    else
        map = map';
    end
end

MAXCHANS = size(map,1);
GRID_SCALE = 2*MAXCHANS+5;

% Read the channel file
% ---------------------
if ischar( filename ) || isstruct( filename )
    [tmp lb Th Rd] = readlocs(filename);
	Th = pi/180*Th;                               % convert degrees to rads
	[x,y] = pol2cart(Th,Rd);
else
	x = real(filename);
	y = imag(filename);		
    if exist('draw') == 1 && draw ~= 0
        line( [(x-0.01)' (x+0.01)']', [(y-0.01)' (y+0.01)']');
        line( [(x+0.01)' (x-0.01)']', [(y-0.01)' (y+0.01)']');
    end
end;	

% locates nearest position of electrod in the grid 
% ------------------------------------------------
xi = linspace(-0.5,0.5,GRID_SCALE);   % x-axis description (row vector)
yi = linspace(-0.5,0.5,GRID_SCALE);   % y-axis description (row vector)
for i=1:MAXCHANS
   [useless_var horizidx(i)] = min(abs(x(i) - xi));    % find pointers to electrode
   [useless_var vertidx(i)] = min(abs(y(i) - yi));     % positions in Zi
end
   
% Compute gradient
% ----------------
sumLaplac2D = zeros(GRID_SCALE, GRID_SCALE);
for i=1:size(map,2) 
   	[Xi,Yi,Zi] = griddata(y,x,map(:,i),yi',xi, 'v4');   % interpolate data
   	laplac2D = del2(Zi);
	positions = horizidx + (vertidx-1)*GRID_SCALE;

	laplac(:,i) = laplac2D(positions(:));
	sumLaplac2D = sumLaplac2D + laplac2D;

	% Draw gradient
	% -------------
	if exist('draw') == 1 && draw ~= 0
		if size(map,2) > 1
			subplot(ceil(sqrt(size(map,2))), ceil(sqrt(size(map,2))), i);
		end
		plot(y, x, 'x', 'Color', 'black', 'markersize', 5); hold on
		contour(Xi, Yi, laplac2D); hold off;
		%line( [(x-0.01)' (x+0.01)']', [(y-0.01)' (y+0.01)']', 'Color','black');
		%line( [(x+0.01)' (x-0.01)']', [(y-0.01)' (y+0.01)']', 'Color','black');
		title( int2str(i) );
	end
end;                                                             %

return;

