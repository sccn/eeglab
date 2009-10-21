%  std_createclust()  - dreate a new empty cluster.  After creation, components 
%                       may be (re)assigned to it using std_movecomp().
% Usage:
%                    >> [STUDY] = std_createclust(STUDY, ALLEEG, name);
% Inputs:
%   STUDY    - STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG   - vector of EEG datasets included in the STUDY, typically created 
%              using load_ALLEEG().
% Optional inputs:
%   name     - ['string'] name of the new cluster {default: 'Cls #', where 
%              '#' is the next available cluster number}
% Outputs:
%   STUDY    - the input STUDY set structure modified with the new cluster.
%
% Example: 
%            >> [STUDY] = std_createclust(STUDY, ALLEEG, 'eye_movements');
%            % Create a new cluster named 'eye_movements'.
%
%  See also  pop_clustedit(), std_movecomp()
%
% Authors:  Hilit Serby, Arnaud Delorme, Scott Makeig, SCCN, INC, UCSD, June, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, June 07, 2005, hilit@sccn.ucsd.edu
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

function [STUDY, clusters] = std_createclust(STUDY, ALLEEG, varargin)

if nargin< 2
    help std_createclust;
    return;
end;

% decoding options
% ----------------
options = {};
if length(varargin) < 2
    if isempty(varargin) | strcmpi(varargin,'')
        options = 'Cls';
    else
        options = { 'name' varargin{1} };
    end
else
    options =  varargin;
end;
opt = finputcheck(options, { 'name'       'string'   []  '';
                             'subjects'   'cell'     []  {};
                             'datasets'   {'cell' 'integer'}  { [] [] } [];
                             'components' {'cell' 'integer'}  { [] [] } [] }, 'std_createclust');
if isstr(opt), error(opt); end;
    
% Find out the highst cluster id number (in cluster name), to find
% next available cluster index
max_id = 0;
if ~isfield(STUDY, 'cluster'), STUDY.cluster = []; end;
for k = 1:length(STUDY.cluster)
    ti = strfind(STUDY.cluster(k).name, ' ');
    clus_id = STUDY.cluster(k).name(ti(end) + 1:end);
    max_id = max(max_id, str2num(clus_id));
end
max_id = max_id + 1;
opt.name = sprintf('%s %d', opt.name, max_id);
clustind = length(STUDY.cluster)+1;
% Initialize the new cluster fields.
if length(STUDY.cluster) > 0
    STUDY.cluster(clustind).parent{1} = STUDY.cluster(1).name;
    if ~iscell(STUDY.cluster(1).child)
         STUDY.cluster(1).child = { opt.name };
    else STUDY.cluster(1).child = { STUDY.cluster(1).child{:} opt.name };
    end;
else  
    STUDY.cluster(clustind).parent{1} = 'manual'; % update parent cluster if exists.
end;
STUDY.cluster(clustind).name = opt.name;
STUDY.cluster(clustind).child = [];
STUDY.cluster(clustind).comps = [];
STUDY.cluster(clustind).sets = [];
STUDY.cluster(clustind).algorithm = [];
STUDY.cluster(clustind).centroid = [];
STUDY.cluster(clustind).preclust.preclustparams = [];
STUDY.cluster(clustind).preclust.preclustdata = [];

if (~isempty(opt.datasets) | ~isempty(opt.subjects)) & ~isempty(opt.components)
    
    % convert subjects to dataset indices
    % -----------------------------------
    if ~isempty(opt.subjects)
        if length(opt.subjects) ~= length(opt.components)
            error('If subjects are specified, the length of the cell array must be the same as for the components');
        end;
        alls = { ALLEEG.subject };
        for index = 1:length(opt.subjects)
            tmpinds = strmatch(opt.subjects{index}, alls, 'exact');
            if isempty(tmpinds)
                error('Cannot find subject');
            end;
            opt.datasets(1:length(tmpinds),index) = tmpinds;
        end;
        opt.datasets(opt.datasets(:) == 0) = NaN;
    end;
    
    % deal with cell array inputs
    % ---------------------------
    if iscell(opt.components)
        newcomps = [];
        newdats  = [];
        for ind1 = 1:length(opt.components)
            for ind2 = 1:length(opt.components{ind1})
                if iscell(opt.datasets)
                     newdats  = [ newdats  opt.datasets{ind1}' ];
                else newdats  = [ newdats  opt.datasets(:,ind1) ];
                end;
                newcomps = [ newcomps opt.components{ind1}(ind2) ];
            end;
        end;
        opt.datasets   = newdats;
        opt.components = newcomps;
    end;
    
    % create .sets, .comps, .setinds, .allinds fields
    % -----------------------------------------------
    [tmp setinds allinds] = std_setcomps2cell( STUDY, opt.datasets, opt.components);
    STUDY.cluster(clustind).setinds = setinds;
    STUDY.cluster(clustind).allinds = allinds;
    STUDY.cluster(clustind) = std_cell2setcomps( STUDY, ALLEEG, clustind); 
    STUDY.cluster(clustind) = std_setcomps2cell( STUDY, clustind);
    %[ STUDY.cluster(finalinds(ind)) setinds allinds ] =
        %std_setcomps2cell(STUDY, finalinds(ind));
end;
