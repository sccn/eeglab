% gabor2d() - generate a two-dimensional gabor matrice.
%
% Usage:
%   >> [ matrix ] = gabor2d(rows, columns);
%   >> [ matrix ] = gabor2d( rows, columns, freq, ...
%                             angle, sigmaR, sigmaC, meanR, meanC, dephase, cut)
% Example :
%	>> imagesc(gabor2d( 50, 50))
%
% Inputs:
%   rows        - number of rows 
%   columns     - number of columns 
%   freq        - frequency of the sinusoidal function in degrees (default: 360/rows)
%   angle       - angle of rotation of the resulting 2-D array in
%                 degrees of angle {default: 0}.
%   sigmaR      - standard deviation for rows {default: rows/5}
%   sigmaC      - standard deviation for columns {default: columns/5}
%   meanR       - mean for rows {default: center of the row}
%   meanC       - mean for columns {default: center of the column}
%   dephase     - phase offset in  degrees {default: 0}. A complex Gabor wavelet 
%                 can be build using gabor2dd(...., 0) + i*gabor2d(...., 90), 
%                 0 and 90 being the phase offset of the real and imaginary parts
%   cut	        - percentage (0->1) of maximum value below which to remove values 
%                 from the matrix {default: 0}
% Output:
%   matrix - output gabor matrix
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001

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

function mat = gabor2d( sizeX, sizeY, freq, angle, sigmaX, sigmaY, meanX, ...
meanY, dephase, cut);

if nargin < 2
	help gabor2d
	return; 
end
if nargin < 3
	freq = 360/sizeX;
end
if nargin < 4
	angle = 0;
end
if nargin < 5
	sigmaX = sizeX/5;
end
if nargin < 6
	sigmaY = sizeY/5;
end
if nargin < 7
	meanX = (sizeX+1)/2;
end
if nargin < 8
	meanY = (sizeY+1)/2;
end
if nargin < 9
	dephase = 0;
end
if nargin < 10
	cut = 0;
end
freq = freq/180*pi;

X = linspace(1, sizeX, sizeX)'* ones(1,sizeY);
Y = ones(1,sizeX)'   		  * linspace(1, sizeY, sizeY);
%[-sizeX/2:sizeX/2]'*ones(1,sizeX+1);
%Y = ones(1,sizeY+1)'   *[-sizeY/2:sizeY/2];

rotatedmat = ((X-meanX)+i*(Y-meanY)) * exp(i*angle/180*pi);
mat = sin(real(rotatedmat)*freq + dephase/180*pi).*exp(-0.5*(  ((X-meanX)/sigmaX).*((X-meanX)/sigmaX)...
				+((Y-meanY)/sigmaY).*((Y-meanY)/sigmaY)))... 
            			/((sigmaX*sigmaY)^(0.5)*pi); 

if cut > 0
	maximum = max(max(mat))*cut;
	I = find(mat < maximum);
	mat(I) = 0;
end

return;

% other solution
% --------------

for X = 1:sizeX
    for Y = 1:sizeY
        mat(X,Y) = sin(real((X-meanX+j*(Y-meanY))*exp(i*angle/180*pi))*freq + dephase/180*pi) ...
            .*exp(-0.5*(  ((X-meanX)/sigmaX).*((X-meanX)/sigmaX)...
                          +((Y-meanY)/sigmaY).*((Y-meanY)/sigmaY)))... 
            			/((sigmaX*sigmaY)^(0.5)*pi); 
    end
end

return;
