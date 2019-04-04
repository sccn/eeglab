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

function mat = gauss3d( sizeX, sizeY, sizeZ, sigmaX, sigmaY, sigmaZ, meanX, meanY, meanZ, cut);

if nargin < 2
	help gauss2d
	return; 
end
if nargin < 4
	sigmaX = sizeX/5;
end
if nargin < 5
	sigmaY = sizeY/5;
end
if nargin < 6
	sigmaZ = sizeZ/5;
end
if nargin < 7
	meanX = (sizeX+1)/2;
end
if nargin < 8
	meanY = (sizeY+1)/2;
end
if nargin < 9
	meanZ = (sizeZ+1)/2;
end
if nargin < 10
	cut = 0;
end

[X,Y,Z] = ndgrid(1:sizeX,1:sizeY,1:sizeZ);

mat = exp(-0.5*(  ((X-meanX)/sigmaX).*((X-meanX)/sigmaX)...
				 +((Y-meanY)/sigmaY).*((Y-meanY)/sigmaY)... 
				 +((Z-meanZ)/sigmaZ).*((Z-meanZ)/sigmaZ)))... 
            			/((sigmaX*sigmaY*sigmaZ)^(0.5)*pi); 

if cut > 0
	maximun = max(mat(:))*cut;
	I = find(mat < maximun);
	mat(I) = 0;
end

return;

