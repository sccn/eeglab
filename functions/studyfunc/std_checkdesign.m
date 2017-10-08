% std_checkdesign() - Check if STUDY design is compatible with EEGLAB
%                     2-categorical factor statistics
% Usage:
%   >> [STUDY] = std_checkdesign(STUDY, designind);
%
% Inputs:
%   STUDY      - EEGLAB STUDY set
%   designind  - [integer] design index
%
% Authors: Arnaud Delorme

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, 2017, hilit@sccn.ucsd.edu
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

function retVal = std_checkdesign(STUDY, designind)

if nargin < 2
    help std_checkdesign;
    return;
end

retVal = 1;
if length(STUDY.design(designind).variable) > 2 || ...
        ( ~isempty(STUDY.design(designind).variable) && ...
    any(strmatch('continuous', {STUDY.design(designind).variable.vartype })))
    warndlg2( [ 'These plotting functions cannot process this design.' 10 ...
        'This design either has more than 2 categorical variables' 10 ...
        'or it uses one continuous variable. To process this' 10 ...
        'design, use the new LIMO EEG tools.' ]);
    retVal = 0;
end

