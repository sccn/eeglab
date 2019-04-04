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
% See also: envtopo()

% Copyright (C) Scott Makeig & Arnaud Delorme - CNL / Salk Institute, La Jolla 8/8/97
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

% 2001 - extrapolation -ad
% 01-25-02 reformated help & license -ad 

function envdata = env(data, timelimits, timearray )

maxdata = max(data);
mindata = min(data);

% extrapolate these values if necessary
% -------------------------------------
if nargin > 2
        timelimits = timelimits(:)';  % make row vector
        if size(timelimits,2)~=2 || size(timelimits,2)~=2
           error('timelimits array must be a [start_time, end_time] vector')
        end
	X = linspace(timelimits(1),timelimits(2),length(maxdata));   % x-axis description (row vector)
	Y = ones(1,size(X,2));
        if size(timearray,1)>1 && size(timearray,2)>1
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
