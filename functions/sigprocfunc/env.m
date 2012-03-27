% env() - return envelope of rows of a data matrix, or optionally
%         of the data interpolated to a different sampling rate.
% Usage:
%   >> envdata = env(data);
%   >> envdata = env(data, timelimits, timearray);
%
% Inputs:
%   data       - (nchannels,ntimepoints) data array
%
% Optional Inputs:
%   timelimits - (start_time, end_time) timelimits (default: none required)
%   timearray  - Optional times array to interpolate the data (default: none)
%
% Outputs:
%   envdata    - A (2,timepoints) array containing the "envelope" of 
%                a multichannel data set = the maximum and minimum values,
%                across all the channels, at each time point. That is,
%                   >> envdata = [max(data');min(data')];
%
% Author: Scott Makeig & Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: ENVTOPO

% Scott Makeig & Arnaud Delorme - CNL / Salk Institute, La Jolla 8/8/97
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

% 2001 - extrapolation -ad
% 01-25-02 reformated help & license -ad 

function envdata = env(data, timelimits, timearray )

maxdata = max(data);
mindata = min(data);

% extrapolate these values if necessary
% -------------------------------------
if nargin > 2
        timelimits = timelimits(:)';  % make row vector
        if size(timelimits,2)~=2 | size(timelimits,2)~=2
           error('timelimits array must be a [start_time, end_time] vector')
        end
	X = linspace(timelimits(1),timelimits(2),length(maxdata));   % x-axis description (row vector)
	Y = ones(1,size(X,2));
        if size(timearray,1)>1 & size(timearray,2)>1
           error('timearray must be a vector')
        end
	Xi = timearray(:)';   % make a row vector
	Yi = ones(1,length(timearray));

    try
        [tmp1,tmp2,Zi] = griddata(Y, X, maxdata, Yi, Xi, 'v4');   % interpolate data
        maxdata = Zi;
        [tmp1,tmp2,Zi] = griddata(Y, X, mindata, Yi, Xi, 'v4');   % interpolate data
        mindata = Zi;
    catch
        disp('Warning, "v4" interpolation failed, using linear interpolation instead (probably running Octave)');
        [tmp1,tmp2,Zi] = griddata(Y, X, maxdata, Yi, Xi, 'linear');   % interpolate data
        maxdata = Zi;
        [tmp1,tmp2,Zi] = griddata(Y, X, mindata, Yi, Xi, 'linear');   % interpolate data
        mindata = Zi;
    end
end;	

envdata = [maxdata;mindata];
return;
