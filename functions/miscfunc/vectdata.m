% vectdata() - vector data interpolation with optional moving 
%              average.
%
% Usage:
%   >> [interparray timesout] = vectdata( array, timesin, 'key', 'val', ... );
%
% Inputs:
%   array      - 1-D or 2-D float array. If 2-D, the second dimension
%                only is interpolated.
%   timesin    - [float vector] time point indices. Same dimension as
%                the interpolated dimension in array.
%
% Optional inputs
%   'timesout' - [float vector] time point indices for interpolating
%                data.
%   'method'   - method for interpolation
%        'linear'    -> Triangle-based linear interpolation (default).
%        'cubic'     -> Triangle-based cubic interpolation.
%        'nearest'   -> Nearest neighbor interpolation.
%        'v4'        -> MATLAB 4 griddata method.
%   'average'  - [real] moving average in the dimension of timesin
%                note that extreme values might be inacurate. Default
%                none or [].
%
% Outputs:
%   interparray - interpolated array
%   timesout    - output time points
%
% Author: Arnaud Delorme, CNL / Salk Institute, 20 Oct 2002
%
% See also: griddata()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [interparray, timesout] = vectdata( array, timevect, varargin );

if nargin < 3
    help vectdata;
    return;
end;

g = finputcheck( varargin, { 'timesout'   'real'  []           [];
                             'average'    'real'  []           [];
                             'gauss'      'real'  []           [];
                             'method'     'string' { 'linear' 'cubic' 'nearest' 'v4' } 'linear'});
if isstr(g), error(g); end;

if size(array,2) == 1
    array = transpose(array);
end;

if ~isempty(g.average)
    if any( (timevect(2:end) -timevect(1:end-1)) - mean(timevect)) % not uniform values
        fprintf('Data has to be interpolated uniformly for moving average\n');
        minspace = min(timevect(2:end) -timevect(1:end-1));
        newtimevect = linspace(timevect(1), timevect(end), ceil((timevect(end)-timevect(1))/minspace)); 
        array = interpolate( array, timevect, newtimevect, g.method);
        timevect = newtimevect;
    end;
    g.average = round(g.average/(timevect(2)-timevect(1)));
    fprintf('Moving average updated to %d\n', g.average);
    array = convolve(array, ones(1, g.average));
end;

interparray = interpolate( array, timevect, g.timesout, g.method);
timesout = g.timesout;

% interpolation function
% ----------------------
function [interparray] = interpolate( array, timesin, timesout, method);
    interparray = zeros(size(array,1), length(timesout));
    for index = 1:size(array,1)
        tmpa = [array(1,:) ; array(1,:)];
        
        [Xi,Yi,Zi] = griddata(timesin, [1 2]', tmpa, timesout, [1 2]', method); % Interpolate data
        interparray(index,:) = Zi(1,:);
    end;
