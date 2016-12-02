% std_combtrialinfo() - Combine and integrate all the info into the
%                       field trialinfo
%
% Usage:
%   >>   trialinfo = std_combtrialinfo(datasetinfo, [1 2]);
%   >>   trialinfo = std_combtrialinfo(datasetinfo, 'S01');
%
% Inputs:
%      datasetinfo         -  Datasetinfo structure from STUDY or
%                             datafiles
%      inds                - Either subject name  or indices of subject
%                            in datasetinfo structure.
%      trials              - All EEG dataset trial numbers
%
% Outputs:
%    trialinfo  - Updated trialinfo structure
%
% See also:
%   std_plotinfocluster
%
% Author: Arnaud Delorme
%         Ramon Martinez-Cancino
%
% Copyright (C) 2015  Arnaud Delorme, Ramon Martinez-Cancino,  INC, SCCN
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

function trialinfo = std_combtrialinfo(datasetinfo, inds, trials)

if nargin < 1
    help std_combtrialinfo;
    return;
end;

% Inds or subjectname
if isnumeric(inds)
    if ~isempty(find(~ismember(inds(:), [datasetinfo.index]))) %#ok<EFIND>
         error(['std_combtrialinfo error: Indices '' ' num2str(inds) ' '' aout of range']);
    end
elseif ischar(inds)
    if ~ismember(inds, {datasetinfo.subject})
        error(['std_combtrialinfo error: Subject name '' ' inds ' '' is invalid']);
    else
        % Looking for indices of subject 'ind'
        inds = find(cellfun(@(x) strcmp(x,inds), {datasetinfo.subject}));
    end
    
end

if ~isfield(datasetinfo, 'trialinfo')
    trialinfo = struct([]);
    nvals = [ 1 cumsum( [trials( [ datasetinfo(inds).index ]) ])+ 1 ];
else
    % check if duration field is present
    for iDat = inds(:)'
        if ~isfield(datasetinfo(iDat).trialinfo, 'duration') datasetinfo(iDat).trialinfo(1).duration = []; end;
    end;
    try
        trialinfo = [ datasetinfo(inds(:)').trialinfo ];
    catch
        datasetinfo(inds(:)').trialinfo
        error( [ 'Trial information field differ for subject ' datasetinfo(inds(1)).subject ' (see above)' ]);
    end;
    nvals     = [ 1 cumsum(cellfun(@length, { datasetinfo(inds).trialinfo }))+1 ];
end;

for iDat = 1:length(inds)
    [trialinfo(nvals(iDat):nvals(iDat+1)-1).condition] = deal( datasetinfo(inds(iDat)).condition );
    [trialinfo(nvals(iDat):nvals(iDat+1)-1).group    ] = deal( datasetinfo(inds(iDat)).group     );
    [trialinfo(nvals(iDat):nvals(iDat+1)-1).session  ] = deal( datasetinfo(inds(iDat)).session   );
end;