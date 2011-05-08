% laplac2d() - generate a 2 dimensional gaussian matrice
%
% Usage :
%    >> [ gaussmatrix ] = laplac2d( rows, columns, sigma, ...
%                                       meanR, meanC, cut)
%
% Example :
%   >> laplac2d( 5, 5)
%
% Inputs:
%   rows    - number of rows 
%   columns - number of columns 
%   sigma   - standart deviation (default: rows/5)
%   meanR   - mean for rows (default: center of the row)
%   meanC   - mean for columns (default: center of the column)
%   cut	    - percentage (0->1) of the maximum value for removing 
%             values in the matrix (default: 0) 
%
% Note: this function implements a simple laplacian for exploratory
%       research. For a more rigorous validated approach use the freely 
%       available Current Source Density Matlab toolbox.
%
% See also: eeg_laplac()
%
% Author: Arnaud Delorme, CNL, Salk Institute, 2001

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

function mat = laplac2d( sizeX, sizeY, sigma, meanX, meanY, cut);

if nargin < 2
	help laplac2d
	return; 
end;
if nargin < 3
	sigma = sizeX/5;
end;
if nargin < 4
	meanX = (sizeX+1)/2;
end;
if nargin < 5
	meanY = (sizeY+1)/2;
end;
if nargin < 6
	cut = 0;
end;

X = linspace(1, sizeX, sizeX)'* ones(1,sizeY);
Y = ones(1,sizeX)'   		  * linspace(1, sizeY, sizeY);
%[-sizeX/2:sizeX/2]'*ones(1,sizeX+1);
%Y = ones(1,sizeY+1)'   *[-sizeY/2:sizeY/2];

r2 = (X-meanX).*(X-meanX) + (Y-meanY).*(Y-meanY);
sigma2 = sigma*sigma;

mat = - exp(-0.5*r2/sigma2) .* ((r2 - sigma2)/(sigma2*sigma2)); 
% zeros crossing at r = -/+ sigma;
% mat = r2;
if cut > 0
	maximun = max(max(mat))*cut;
	I = find(mat < maximun);
	mat(I) = 0;
end;

return;
