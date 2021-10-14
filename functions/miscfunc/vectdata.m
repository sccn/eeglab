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
%                note that extreme values might be inaccurate (see 'borders'). 
%                Default none or [].
%   'avgtype'  - ['const'|'gauss'] use a const value when averaging (array of 
%                ones) or a gaussian window. Default is 'const'.
%   'border'   - ['on'|'off'] correct border effect when smoothing.
%                default is 'off'.
%
% Outputs:
%   interparray - interpolated array
%   timesout    - output time points
%
% Author: Arnaud Delorme, CNL / Salk Institute, 20 Oct 2002
%
% See also: griddata()

% Copyright (C) 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [interparray, timesout] = vectdata( array, timevect, varargin );

if nargin < 3
    help vectdata;
    return;
end

g = finputcheck( varargin, { 'timesout'   'real'  []           [];
                             'average'    'real'  []           [];
                             'gauss'      'real'  []           [];
                             'border'     'string' { 'on','off' } 'off';
                             'avgtype'    'string' { 'const','gauss' } 'const';
                             'method'     'string' { 'linear','cubic','nearest','v4' } 'linear'});
if ischar(g), error(g); end

if size(array,2) == 1
    array = transpose(array);
end

if ~isempty(g.average)
    timediff = timevect(2:end) -timevect(1:end-1);
    if any( (timediff - mean(timediff)) > 1e-8 ) % not uniform values
        fprintf('Data has to be interpolated uniformly for moving average\n');
        minspace = mean(timediff);
        newtimevect = linspace(timevect(1), timevect(end), ceil((timevect(end)-timevect(1))/minspace)); 
        array = interpolate( array, timevect, newtimevect, g.method);
        timevect = newtimevect;
    end
    oldavg = g.average;
    g.average = round(g.average/(timevect(2)-timevect(1)));
    if oldavg ~= g.average
        fprintf('Moving average updated from %3.2f to %3.2f (=%d points)\n', ...
                oldavg, g.average*(timevect(2)-timevect(1)), g.average);
    end
    if strcmpi(g.border, 'on')
        if strcmpi(g.avgtype, 'const')
            array = convolve(array, ones(1, g.average));
        else
            convolution = gauss2d(1,g.average,1,round(0.15*g.average));
            array = convolve(array, convolution);
        end
    else
        if strcmpi(g.avgtype, 'const')
            array = conv2(array, ones(1, g.average)/g.average, 'same');
        else
            convolution = gauss2d(1,g.average,1,round(0.15*g.average));
            array = conv2(array, convolution/sum(convolution), 'same');
        end
    end
end

interparray = interpolate( array, timevect, g.timesout, g.method);
timesout = g.timesout;

% interpolation function
% ----------------------
function [interparray] = interpolate( array, timesin, timesout, method);
    interparray = zeros(size(array,1), length(timesout));
    for index = 1:size(array,1)
        tmpa = [array(index,:) ; array(index,:)];
        
        [Xi,Yi,Zi] = griddata(timesin, [1 2]', tmpa, timesout, [1 2]', method); % Interpolate data
        interparray(index,:) = Zi(1,:);
    end
