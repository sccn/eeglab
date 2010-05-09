% std_selectdesign() - select an existing STUDY design. 
%                      Use std_makedesign() to add a new STUDY.design.
%
% Usage:
%   >> [STUDY] = std_selectdesign(STUDY, ALLEEG, designind);
%
% Inputs:
%   STUDY     - EEGLAB STUDY structure
%   STUDY     - EEGLAB ALLEEG structure
%   designind - desired (existing) design index
%
% Outputs:
%   STUDY     - EEGLAB STUDY structure with currentdesign set to the input design index. 
%
% Author: Arnaud Delorme, Institute for Neural Computation, UCSD, 2010-

% Copyright (C) Arnaud Delorme, arno@salk.edu
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

function STUDY = std_selectdesign(STUDY, ALLEEG, designind);

if nargin < 3
    help std_selectdesign;
    return;
end;

if designind < 1 || designind > length(STUDY.design) || isempty(STUDY.design(designind).name)
    disp('Cannot select an empty STUDY.design');
    return;
end;
STUDY.currentdesign = designind;
STUDY = std_rmalldatafields( STUDY );

% remake setinds and allinds
% --------------------------
STUDY = std_changroup(STUDY, ALLEEG, [], 'on'); % with interpolation HAVE TO FIX THAT

% update the component indices
% ----------------------------
STUDY.cluster(1).setinds = {};
STUDY.cluster(1).allinds = {};
for index = 1:length(STUDY.cluster)
    STUDY.cluster(index) = std_setcomps2cell(STUDY, index);
end;
