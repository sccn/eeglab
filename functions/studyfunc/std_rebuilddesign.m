% std_rebuilddesign - reduild design structure when datasets have been 
%                     removed or added.
% Usage:
%    STUDY = std_rebuilddesign(STUDY, ALLEEG);
%    STUDY = std_rebuilddesign(STUDY, ALLEEG, designind);
%
% Inputs:
%   STUDY     - EEGLAB STUDY set
%   ALLEEG    - vector of the EEG datasets included in the STUDY structure 
%
% Optional inputs:
%   designind - [integer>0] indices (number) of the design to rebuild.
%               Default is all.
% Ouput:
%   STUDY     - updated EEGLAB STUDY set
%
% Author: Arnaud Delorme, Institute for Neural Computation UCSD, 2010-

% Copyright (C) Arnaud Delorme, arno@ucsd.edu
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

function STUDY = std_rebuilddesign(STUDY, ALLEEG, designind);

if nargin < 2
    help std_rebuilddesign;
    return;
end;

if nargin < 3
    designind = 1:length(STUDY.design);
end;

disp('Rebuilding designs...');
[indvars indvarvals] = std_getindvar(STUDY);
for indDesign = designind
    
    % find out if some independent variables or independent 
    % variable values have been removed
    tmpdesign = STUDY.design(indDesign);
    for iVar = length(tmpdesign.variable):-1:1
        indVar = strmatch(tmpdesign.variable(iVar).label, indvars, 'exact');
        if isempty(indVar)
            fprintf('Design %d: removing independent variable "%s"\n', indDesign, tmpdesign.variable(iVar).label);
            tmpdesign.variable(iVar) = []; % remove independent var
        else
            tmpvals = myintersect( tmpdesign.variable(iVar).value, indvarvals{indVar} );
            if isempty(tmpvals)
                 fprintf('Design %d: removing independent variable "%s"\n', indDesign, tmpdesign.variable(iVar).label);
                 tmpdesign.variable(iVar) = []; % remove independent var
            else tmpdesign.variable(iVar).value = tmpvals;
            end;
        end;
    end;
    
    STUDY = std_makedesign(STUDY, ALLEEG, indDesign, tmpdesign);
end;

STUDY = std_selectdesign(STUDY, ALLEEG, STUDY.currentdesign);

% take the intersection for independent variables
% -----------------------------------------------
function a = myintersect(a,b);

    if isempty(b) || isempty(a), a = {}; return; end;
    
    for index = length(a):-1:1
        if isstr(a{index})
            res = intersect_bc(a(index), b);
            if isempty(res)
                 a{index} = [];
            else a(index) = res;
            end;
        elseif iscell(a{index})
            a{index} = intersect_bc(a{index}, b);
        elseif isnumeric(a{index})
            a{index} = intersect_bc(a{index}, [ b{:} ]);
        end;
        if isempty(a{index}), a(index) = []; end;
    end;
        