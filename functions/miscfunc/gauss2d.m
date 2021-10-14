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
% Output:
%   gaussmatrix - 2-D gaussian matrix
%
% Author: Arnaud Delorme, CNL/Salk Institute, 2001

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

function mat = gauss2d( sizeX, sizeY, sigmaX, sigmaY, meanX, meanY, cut);

if nargin < 2
	help gauss2d
	return; 
end
if nargin < 3
	sigmaX = sizeX/5;
end
if nargin < 4
	sigmaY = sizeY/5;
end
if nargin < 5
	meanX = (sizeX+1)/2;
end
if nargin < 6
	meanY = (sizeY+1)/2;
end
if nargin < 7
	cut = 0;
end

X = linspace(1, sizeX, sizeX)'* ones(1,sizeY);
Y = ones(1,sizeX)'   		  * linspace(1, sizeY, sizeY);
%[-sizeX/2:sizeX/2]'*ones(1,sizeX+1);
%Y = ones(1,sizeY+1)'   *[-sizeY/2:sizeY/2];

mat = exp(-0.5*(  ((X-meanX)/sigmaX).*((X-meanX)/sigmaX)...
				+((Y-meanY)/sigmaY).*((Y-meanY)/sigmaY)))... 
            			/((sigmaX*sigmaY)^(0.5)*pi); 

if cut > 0
	maximum = max(max(mat))*cut;
	I = find(mat < maximum);
	mat(I) = 0;
end

return;

