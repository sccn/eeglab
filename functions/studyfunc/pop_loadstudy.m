% pop_loadstudy() - load an existing EEGLAB STUDY set of EEG datasets plus 
%                   its corresponding ALLEEG structure. Calls load_alleeg().
% Usage:
%   >> [STUDY ALLEEG] = pop_loadstudy; % pop up a window to collect filename
%   >> [STUDY ALLEEG] = pop_loadstudy( 'key', 'val', ...); % no pop-up
%
% Optional inputs:
%   'filename' - [string] filename of the STUDY set file to load.
%   'filepath' - [string] filepath of the STUDY set file to load.
%
% Outputs:
%   STUDY      - the requested STUDY set structure.
%   ALLEEG     - the corresponding ALLEEG structure containing 
%                the (loaded) STUDY EEG datasets.    
%
% See also: load_alleeg(), pop_savestudy(), pop_createstudy()
%
% Authors: Hilit Serby SCCN, INC, UCSD, September 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, Spetember 2005, hilit@sccn.ucsd.edu
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

% Coding notes: Useful information on functions and global variables used.


function [STUDY, ALLEEG, com] = pop_loadstudy(varargin)

STUDY  = [];
ALLEEG = [];
com = '';
if isempty(varargin)
    [filename, filepath] = uigetfile2('*.study', 'Load a STUDY -- pop_loadstudy()'); 
    if filename(1) == 0, return; end;
    if ~strncmp(filename(end-5:end), '.study',6)
        if isempty(strfind(filename,'.'))
            filename = [filename '.study'];
        else
            filename = [filename(1:strfind(filename,'.')-1) '.study'];
        end
    end
else
    for k = 1:2:length(varargin)
        switch varargin{k}
            case 'filename'
                filename = varargin{k+1};
            case 'filepath'
                filepath = varargin{k+1};
        end
    end
end

if (~isempty(filename)) & (~isempty(filepath))
    STUDYfile = fullfile(filepath,filename);
    try 
        load('-mat', STUDYfile);
    catch
        error(['pop_loadstudy(): STUDY set file ''STUDYfile'' not loaded -- check filename and path']);
    end
    [filepath filename ext] = fileparts(STUDYfile);
    STUDY.filename = [filename ext];
    STUDY.filepath = filepath;
else
    error(['pop_loadstudy(): No STUDY set file provided.']);
end
STUDY.saved = 'yes';
  
ALLEEG = load_alleeg(STUDY);

% Update the pointers from STUDY to the ALLEEG datasets
for k = 1:length(STUDY.datasetinfo)
    STUDY.datasetinfo(k).index = k;
end
       
com = sprintf('[STUDY ALLEEG] = pop_loadstudy(''filename'', ''%s'', ''filepath'', ''%s'');', filename, filepath);