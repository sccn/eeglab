% fieldtripchan2eeglab() - convert Fieldtrip channel location structure
%                          to EEGLAB channel location structure
%
% Usage:
%   >> chanlocs = fieldtripchan2eeglab( fieldlocs );
%
% Inputs:
%   fieldlocs - Fieldtrip channel structure. See help readlocs()
%
% Outputs:
%   chanlocs  - EEGLAB channel location structure.
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2006-
%
% See also: readlocs()

% Copyright (C) 2003 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function chanlocs = fieldtripchan2eeglab( loc );
    
    if nargin < 1
        help fieldtripchan2eeglab;
        return;
    end;
    
    chanlocs = struct('labels', loc.label', 'X', mattocell(loc.pnt(:,1)'), ...
                                            'Y', mattocell(loc.pnt(:,2)'), ...
                                            'Z', mattocell(loc.pnt(:,3)'));
    chanlocs = convertlocs(chanlocs, 'cart2all');
