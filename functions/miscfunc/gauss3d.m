% gauss3d() - generate a 3-dimensional gaussian matrix
%
% Usage:
%   >> [ gaussmatrix ] = gauss2d( nX, nY, nZ);
%   >> [ gaussmatrix ] = gauss2d( nX, nY, nZ, ...
%                              sigmaX, sigmaY, sigmaZ, ...
%                              centerX, centerY, centerZ, mask)
%
% Example:
%   >> gauss3d(3,3,3); % generate a 3x3x3 gaussian matrix
%
% Inputs:
%   nX    - number of values in first dimension
%   nY    - number of values in second dimension
%   nZ    - number of values in third dimension
%   sigmaX  - width of standard deviation in first dim (default: nX/5)
%   sigmaY  - width of standard deviation in second dim (default: nY/5)
%   sigmaZ  - width of standard deviation in third dim (default: nZ/5)
%   centerX - location of center (default: nX/2)
%   centerY - location of center (default: nY/2)
%   centerZ - location of center (default: nZ/2)
%   mask    - (0->1) percentage of low values in the matrix to mask 
%             with zeros (default: 0 or none)
%
% Ouput:
%   gaussmatrix - 3-D gaussian matrix
%
% Author: Arnaud Delorme, 2009

% Copyright (C) 2009 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function mat = gauss3d( sizeX, sizeY, sizeZ, sigmaX, sigmaY, sigmaZ, meanX, meanY, meanZ, cut);

if nargin < 2
	help gauss2d
	return; 
end;
if nargin < 4
	sigmaX = sizeX/5;
end;
if nargin < 5
	sigmaY = sizeY/5;
end;
if nargin < 6
	sigmaZ = sizeZ/5;
end;
if nargin < 7
	meanX = (sizeX+1)/2;
end;
if nargin < 8
	meanY = (sizeY+1)/2;
end;
if nargin < 9
	meanZ = (sizeZ+1)/2;
end;
if nargin < 10
	cut = 0;
end;

[X,Y,Z] = ndgrid(1:sizeX,1:sizeY,1:sizeZ);

mat = exp(-0.5*(  ((X-meanX)/sigmaX).*((X-meanX)/sigmaX)...
				 +((Y-meanY)/sigmaY).*((Y-meanY)/sigmaY)... 
				 +((Z-meanZ)/sigmaZ).*((Z-meanZ)/sigmaZ)))... 
            			/((sigmaX*sigmaY*sigmaZ)^(0.5)*pi); 

if cut > 0
	maximun = max(mat(:))*cut;
	I = find(mat < maximun);
	mat(I) = 0;
end;

return;

