% env() - return envelope of rows of a data matrix.
%
% Usage:
%   >> envdata = env( data, timelimits, timearray);
%
% Inputs:
%   data       - nbchannel x points data
%   timelimits - timelimits (default: none)
%   timearray  - time array to extrapolate data (default: none)
%
% Outputs:
%   envdata    - The "envelope" of a multichannel data set is the maximum
%                and minimum value at each time point, i.e.
%                envdata = [rowmax;rowmin];
%
% Author: Scott Makeig & Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: ENVTOPO

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% 2001 - extrapolation -ad
% 01-25-02 reformated help & license -ad 

function envdata = env(data, timelimits, timearray )

maxdata = max(data);
mindata = min(data);

% extrapolate these values if necessary
% -------------------------------------
if nargin > 2
	X = linspace(timelimits(1),timelimits(2),length(maxdata));   % x-axis description (row vector)
	Y = ones(1,size(X,2));
	Xi = [timearray];
	Yi = ones(1,size(timearray,2));

    warning off;
	[tmp1,tmp2,Zi] = griddata(Y, X, maxdata, Yi, Xi, 'invdist');   % interpolate data
	maxdata = Zi;
	[tmp1,tmp2,Zi] = griddata(Y, X, mindata, Yi, Xi, 'invdist');   % interpolate data
	mindata = Zi;
    warning on;
end;	

envdata = [maxdata;mindata];
return;
