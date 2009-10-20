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

% $Log: not supported by cvs2svn $

function [ tmpstruct setinds allinds ] = std_setcomps2cell(STUDY, sets, comps)

if nargin < 3
    tmpstruct = STUDY.cluster(sets);
    alldatasets = tmpstruct.sets;
    allchanorcomp = repmat(tmpstruct.comps, [size(tmpstruct.sets,1) 1]); % old format
else
    tmpstruct     = [];
    alldatasets   = sets;
    allchanorcomp = comps;
    if length(comps) < length(sets(:))
        allchanorcomp = repmat(comps, [size(sets,1) 1]); % old format
    end;
end;

alldatasets   = alldatasets(:)';
allchanorcomp = allchanorcomp(:)';

% get indices for all groups and conditions
% -----------------------------------------
nc = max(length(STUDY.condition),1);
ng = max(length(STUDY.group),1);
allinds = cell( nc, ng );
setinds = cell( nc, ng );
for indtmp = 1:length(alldatasets)
    if ~isnan(alldatasets(indtmp))
        index = alldatasets(indtmp);
        condind = strmatch( STUDY.datasetinfo(index).condition, STUDY.condition, 'exact'); if isempty(condind), condind = 1; end;
        grpind  = strmatch( STUDY.datasetinfo(index).group    , STUDY.group    , 'exact'); if isempty(grpind) , grpind  = 1; end;
        indcellarray = length(allinds{condind, grpind})+1;
    end
    
    % load data
    % ---------
    tmpind = allchanorcomp(indtmp);
    if ~isnan(tmpind)
        allinds{ condind, grpind}(indcellarray) = tmpind;
        setinds{ condind, grpind}(indcellarray) = index;
    end;
end;
tmpstruct.allinds = allinds;
tmpstruct.setinds = setinds;