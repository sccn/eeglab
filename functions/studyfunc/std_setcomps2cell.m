% std_setcomps2cell - convert .sets and .comps to cell array. The .sets and
%                     .comps format is useful for GUI but the cell array
%                     format is used for plotting and statistics.
%            
% Usage:
%   [ struct setinds allinds ] = std_setcomps2cell(STUDY, clustind);
%   [ struct setinds allinds ] = std_setcomps2cell(STUDY, sets, comps);
%
% Author: Arnaud Delorme, CERCO/CNRS, UCSD, 2009-

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

% $Log: std_setcomps2cell.m,v $
% Revision 1.1  2009/10/20 02:28:35  arno
% Updated conversion between sets and indices formats
%

function [ tmpstruct setinds allinds ] = std_setcomps2cell(STUDY, sets, comps)

if nargin < 3
    tmpstruct = STUDY.cluster(sets);
    sets  = tmpstruct.sets;
    comps = tmpstruct.comps; % old format
else
    tmpstruct     = [];
end;
comps   = repmat(comps, [size(sets,1) 1]);
oldsets = sets;
sets    = reshape(sets , 1, size(sets ,1)*size(sets ,2));
comps   = reshape(comps, 1, size(comps,1)*size(comps,2));

% get indices for all groups and conditions
% -----------------------------------------
setinfo       = STUDY.design(STUDY.currentdesign).setinfo;
allconditions = STUDY.design(STUDY.currentdesign).condition;
allgroups     = STUDY.design(STUDY.currentdesign).group;
nc = max(length(allconditions),1);
ng = max(length(allgroups),    1);
allinds = cell( nc, ng );
setinds = cell( nc, ng );

for index = 1:length(setinfo)
    condind = strmatch( setinfo(index).condition, allconditions, 'exact');
    grpind  = strmatch( setinfo(index).group    , allgroups    , 'exact');

    if isempty(allconditions), condind = 1; end;
    if isempty(allgroups),     grpind  = 1; end;

    % get the position in sets where the dataset is
    % if several datasets check that they all have the same
    % ICA and component index
    % -----------------------
    datind  = setinfo(index).setindex;
    ind     = find(datind(1) == sets);
    if ~isempty(ind) && length(datind) > 1
        [ind1 ind2] = find(datind(1) == oldsets);
        columnica   = oldsets(ind2(1),:);
        if ~all(ismember(datind, columnica));
            warning(' ***** Some datasets in the STUDY design have uniform ICA distributions');
        end;
    end;
        
    allinds{ condind, grpind } = [ allinds{ condind, grpind } comps(ind) ];
    setinds{ condind, grpind } = [ setinds{ condind, grpind } repmat(index, [1 length(ind)]) ];
end;
tmpstruct.allinds = allinds;
tmpstruct.setinds = setinds;
