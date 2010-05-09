% gauss2d() - generate a 2-dimensional gaussian matrix
%
% Usage:
%   >> [ gaussmatrix ] = gauss2d( rows, columns, ...
%                              sigmaR, sigmaC, peakR, peakC, mask)
%
% Example:
%   >> imagesc(gauss2d(50, 50)); % image a size (50,50) 2-D gaussian matrix
%
% Inputs:
%   rows    - number of rows in matrix
%   columns - number of columns in matrix
%   sigmaR  - width of standard deviation in rows (default: rows/5)
%   sigmaC  - width of standard deviation in columns (default: columns/5)
%   peakR   - location of the peak in each row (default: rows/2)
%   peakC   - location of the peak in each column (default: columns/2)
%   mask    - (0->1) portion of the matrix to mask with zeros (default: 0)
%
% Ouput:
%   gaussmatrix - 2-D gaussian matrix
%
% Author: Arnaud Delorme, CNL/Salk Institute, 2001

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

function mat = gauss2d( sizeX, sizeY, sigmaX, sigmaY, meanX, meanY, cut);

if nargin < 2
	help gauss2d
	return; 
end;
if nargin < 3
	sigmaX = sizeX/5;
end;
if nargin < 4
	sigmaY = sizeY/5;
end;
if nargin < 5
	meanX = (sizeX+1)/2;
end;
if nargin < 6
	meanY = (sizeY+1)/2;
end;
if nargin < 7
	cut = 0;
end;

X = linspace(1, sizeX, sizeX)'* ones(1,sizeY);
Y = ones(1,sizeX)'   		  * linspace(1, sizeY, sizeY);
%[-sizeX/2:sizeX/2]'*ones(1,sizeX+1);
%Y = ones(1,sizeY+1)'   *[-sizeY/2:sizeY/2];

mat = exp(-0.5*(  ((X-meanX)/sigmaX).*((X-meanX)/sigmaX)...
				+((Y-meanY)/sigmaY).*((Y-meanY)/sigmaY)))... 
            			/((sigmaX*sigmaY)^(0.5)*pi); 

if cut > 0
	maximun = max(max(mat))*cut;
	I = find(mat < maximun);
	mat(I) = 0;
end;

return;

