% laplac2d() - generate a 2 dimensional laplacian matrice
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
%   sigma   - standard deviation (default: rows/5)
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

function mat = laplac2d( sizeX, sizeY, sigma, meanX, meanY, cut);

if nargin < 2
	help laplac2d
	return; 
end
if nargin < 3
	sigma = sizeX/5;
end
if nargin < 4
	meanX = (sizeX+1)/2;
end
if nargin < 5
	meanY = (sizeY+1)/2;
end
if nargin < 6
	cut = 0;
end

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
	maximum = max(max(mat))*cut;
	I = find(mat < maximum);
	mat(I) = 0;
end

return;
