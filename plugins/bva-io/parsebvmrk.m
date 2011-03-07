% parsebvmrk() - convert Brain Vision Data Exchange format marker
%                configuration structure to EEGLAB event structure
%
% Usage:
%   >> EVENT = parsebvmrk(MRK);
%
% Inputs:
%   MRK   - marker configuration structure
%
% Outputs:
%   EVENT - EEGLAB event structure
%
% Author: Andreas Widmann, University of Leipzig, 2007

% Copyright (C) 2007 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
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

% $Id: parsebvmrk.m 37 2007-06-26 12:56:17Z andreaswidmann $

function EVENT = parsebvmrk(MRK)

for idx = 1:size(MRK.markerinfos, 1)
    [mrkType mrkDesc EVENT(idx).latency EVENT(idx).duration  EVENT(idx).channel EVENT(idx).bvtime] = ...
        strread(MRK.markerinfos{idx, 1}, '%s%s%f%d%d%d', 'delimiter', ',');
    EVENT(idx).bvmknum = MRK.markerinfos{idx, 2};

    if strcmpi(mrkType, 'New Segment') || strcmpi(mrkType, 'DC Correction')
        EVENT(idx).type = 'boundary';
    else
        EVENT(idx).type = char(mrkDesc);
    end

    EVENT(idx).code = char(mrkType);
    EVENT(idx).urevent = idx;
end
