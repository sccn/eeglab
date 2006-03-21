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

if isempty(varargin) | strcmpi(varargin,'')
    name = 'Cls';
else
    name = varargin{1};
end
% Find out the highst cluster id number (in cluster name), to find
% next available cluster index
max_id = 0;
for k = 1:length(STUDY.cluster)
    ti = strfind(STUDY.cluster(k).name, ' ');
    clus_id = STUDY.cluster(k).name(ti(end) + 1:end);
    max_id = max(max_id, str2num(clus_id));
end
max_id = max_id + 1;
name = sprintf('%s %d', name, max_id);
STUDY.cluster(end+1).name = name;
% Initialize the new cluster fields.
STUDY.cluster(end).parent{1} = 'manual'; % update parent cluster if exists.
STUDY.cluster(end).child = [];
STUDY.cluster(end).comps = [];
STUDY.cluster(end).sets = [];
STUDY.cluster(end).algorithm = [];
STUDY.cluster(end).centroid = [];
STUDY.cluster(end).preclust.preclustparams = [];
STUDY.cluster(end).preclust.preclustdata = [];



