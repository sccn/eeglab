% gradmap() - compute the gradient of an EEG spatial distribution.
%
% Usage:
%    >> [gradX, gradY ] = gradmap( map, filename, draw )
%
% Inputs:
%    map      - level of activity (size: nbelectrodes x nbChannel)
%    filename - filename (.loc file) countaining the coordinates
%               of the electrodes, or array countaining complex positions 		 
%               Can also be a EEGLAB channel structure.
%    draw     - integer, if not nul draw the gradient (default:0)
%
% Output:
%    gradX    - gradient over X 
%    gradY    - gradient over Y (use cart2pol to get polar coordinates) 
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001

%  This section bastardizes topoplot.m in order to get gradient maps
%  for all of the component maps brought back from ClusMapSpec.m
%  This is done to improve the clustering results.

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

function [gradx, grady] = gradmap( map, filename, draw ) 

if nargin < 2
	help gradmap;
	return;
end

MAXCHANS = size(map,1);
GRID_SCALE = 2*MAXCHANS+5;

% Read the channel file
% ---------------------
if isnumeric( filename )
	x = real(filename);
	y = imag(filename);
else
    tmploc = readlocs( filename);
    fe = find(cellfun('isempty', { tmploc.theta }));
    tmploc(fe) = [];
    map(fe,:) = [];
    [x,y] = pol2cart( [ tmploc.theta ]/180*pi , [ tmploc.radius ] );
end;	

% locates nearest position of electrod in the grid 
% ------------------------------------------------
xi = linspace(-0.5,0.5,GRID_SCALE);   % x-axis description (row vector)
yi = linspace(-0.5,0.5,GRID_SCALE);   % y-axis description (row vector)
horizidx=zeros(1, MAXCHANS); %preallocation
vertidx=zeros(1, MAXCHANS); % preallocation
for i=1:MAXCHANS
   [useless_var horizidx(i)] = min(abs(y(i) - xi));    % find pointers to electrode
   [useless_var vertidx(i)] = min(abs(x(i) - yi));     % positions in Zi
end
   
draw = 1;

% Compute gradient
% ----------------
for i=1:size(map,2) 
   	[Xi,Yi,Zi] = griddata(y,x,map(:,i),yi',xi, 'v4');   % interpolate data
   	[FX,FY] = gradient(Zi);
	positions = horizidx + (vertidx-1)*GRID_SCALE;

	gradx(:,i) = FX(positions(:));
	grady(:,i) = FY(positions(:));

	% Draw gradient
	% -------------
	if exist('draw')
        if exist('imresize')
            subplot(ceil(sqrt(size(map,2))), ceil(sqrt(size(map,2))), i);
            contour(imresize(Zi,0.5)); hold on
            quiver(imresize(FX, 0.5), imresize(FY, 0.5)); hold off
        else
            disp('Warning: cannot plot because imresize function is missing (image processing toolbox)');
        end
	end
end;                                                       

return;

