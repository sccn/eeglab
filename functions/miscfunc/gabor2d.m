% gabor2d() - generate a 2 dimensional gabor matrice.
%
% Usage:
%   >> [ matrix ] = gabor2d( rows, columns, freq, ...
%                 angle, sigmaR, sigmaC, meanR, meanC, dephase, cut)
%
% Example :
%	>> imagesc(gabor2d( 50, 50))
%
% Inputs:
%   rows        - number of rows 
%   columns     - number of columns 
%   freq        - frequency of sinus function in degree (default: 360/rows)
%   angle       - angle for rotation of the resulting 2D array. In
%                 degree of angle. Default is 0.
%   sigmaR      - standart deviation for rows (default: rows/5)
%   sigmaC      - standart deviation for columns (default: columns/5)
%   meanR       - mean for rows (default: center of the row)
%   meanC       - mean for columns (default: center of the column)
%   dephase     - dephase in  degrees (default:0). A Gabor wavelet can be
%                 build using gabor2dd(...., 0) + i*gabor2d(...., 90), 0
%                 and 90 being the dephasage
%   cut	        - percentage (0->1) of the maximum value for removing 
%                 values in the matrix (default: 0)
%
% Ouput:
%   matrix - output gabor matrix
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
% Revision 1.1  2002/04/05 17:39:45  jorn
% Initial revision
%
% 01-25-02 reformated help & license -ad 

function mat = gabor2d( sizeX, sizeY, freq, angle, sigmaX, sigmaY, meanX, ...
meanY, dephase, cut);

if nargin < 2
	help gabor2d
	return; 
end;
if nargin < 3
	freq = 360/sizeX;
end;
if nargin < 4
	angle = 0;
end;
if nargin < 5
	sigmaX = sizeX/5;
end;
if nargin < 6
	sigmaY = sizeY/5;
end;
if nargin < 7
	meanX = (sizeX+1)/2;
end;
if nargin < 8
	meanY = (sizeY+1)/2;
end;
if nargin < 9
	dephase = 0;
end;
if nargin < 10
	cut = 0;
end;
freq = freq/180*pi;

X = linspace(1, sizeX, sizeX)'* ones(1,sizeY);
Y = ones(1,sizeX)'   		  * linspace(1, sizeY, sizeY);
%[-sizeX/2:sizeX/2]'*ones(1,sizeX+1);
%Y = ones(1,sizeY+1)'   *[-sizeY/2:sizeY/2];

rotatedmat = ((X-meanX)+i*(Y-meanY)) * exp(i*angle/180*pi);
mat = sin(real(rotatedmat)*freq + dephase/180*pi).*exp(-0.5*(  ((X-meanX)/sigmaX).*((X-meanX)/sigmaX)...
				+((Y-meanY)/sigmaY).*((Y-meanY)/sigmaY)))... 
            			/((sigmaX*sigmaY)^(0.5)*pi); 

return;

for X = 1:sizeX
    for Y = 1:sizeY
        mat(X,Y) = sin(real((X-meanX+j*(Y-meanY))*exp(i*angle/180*pi))*freq + dephase/180*pi) ...
            .*exp(-0.5*(  ((X-meanX)/sigmaX).*((X-meanX)/sigmaX)...
                          +((Y-meanY)/sigmaY).*((Y-meanY)/sigmaY)))... 
            			/((sigmaX*sigmaY)^(0.5)*pi); 
    end;
end;


if cut > 0
	maximun = max(max(mat))*cut;
	I = find(mat < maximun);
	mat(I) = 0;
end;

return;
