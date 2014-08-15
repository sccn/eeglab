% std_setcomps2cell - convert .sets and .comps to cell array. The .sets and
%                     .comps format is useful for GUI but the cell array
%                     format is used for plotting and statistics.
%            
% Usage:
%   [ struct setinds allinds ] = std_setcomps2cell(STUDY, clustind);
%   [ struct setinds allinds ] = std_setcomps2cell(STUDY, sets, comps);
%   [ struct setinds allinds measurecell] = std_setcomps2cell(STUDY, sets, comps, measure);
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

function [ tmpstruct setinds allinds measurecell] = std_setcomps2cell(STUDY, sets, comps, measure, generateerror)

if nargin < 5
    generateerror = 0;
end;
if nargin == 4 && length(measure) == 1
    generateerror = measure;
    measure       = [];
end;
if nargin < 3
    tmpstruct = STUDY.cluster(sets);
    sets  = tmpstruct.sets;
    comps = tmpstruct.comps; % old format
else
    tmpstruct     = [];
end;
if nargin < 4 || isempty(measure)
    measure = comps;
end;
measure = repmat(measure, [size(sets,1) 1]);
comps   = repmat(comps  , [size(sets,1) 1]);
oldsets = sets;
sets    = reshape(sets , 1, size(sets ,1)*size(sets ,2));
measure = reshape(measure, 1, size(measure,1)*size(measure,2));
comps   = reshape(comps  , 1, size(comps  ,1)*size(comps  ,2));

% get indices for all groups and conditions
% -----------------------------------------
setinfo       = STUDY.design(STUDY.currentdesign).cell;
allconditions = STUDY.design(STUDY.currentdesign).variable(1).value;
allgroups     = STUDY.design(STUDY.currentdesign).variable(2).value;
nc = max(length(allconditions),1);
ng = max(length(allgroups),    1);
allinds     = cell( nc, ng );
setinds     = cell( nc, ng );
measurecell = cell( nc, ng );

for index = 1:length(setinfo)
    % get index of independent variables
    % ----------------------------------
    condind = std_indvarmatch( setinfo(index).value{1}, allconditions);
    grpind  = std_indvarmatch( setinfo(index).value{2}, allgroups    );
    if isempty(allconditions), condind = 1; end;
    if isempty(allgroups),     grpind  = 1; end;

    % get the position in sets where the dataset is
    % if several datasets check that they all have the same
    % ICA and component index
    % -----------------------
    datind  = setinfo(index).dataset;
    ind     = find(datind(1) == sets);
    if ~isempty(ind) && length(datind) > 1
        [ind1 ind2] = find(datind(1) == oldsets);
        columnica   = oldsets(:,ind2(1));
        if ~all(ismember(datind, columnica));
            disp('Warning: STUDY design combines datasets with different ICA - use ICA only for artifact rejection');
        end;
    end;
        
    measurecell{ condind, grpind } = [ measurecell{ condind, grpind } measure(ind) ];
    allinds{     condind, grpind } = [ allinds{     condind, grpind } comps(  ind) ];
    setinds{     condind, grpind } = [ setinds{     condind, grpind } repmat(index, [1 length(ind)]) ];
end;
tmpstruct.allinds = allinds;
tmpstruct.setinds = setinds;

if generateerror && isempty(setinds{1})
    error( [ 'Some datasets not included in preclustering' 10 ... 
             'because of partial STUDY design. You need to' 10 ...
             'use a STUDY design that includes all datasets.' ]);
end;

