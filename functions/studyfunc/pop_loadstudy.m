% pop_loadstudy() - load an existing STUDY set and its corresponding ALLEEG structure
%
% Usage:
%   >> [STUDY ALLEEG] = pop_loadstudy; % uses an interactive pop-up window
%                                      % to select the STUDY set to load
%   >> [STUDY ALLEEG] = pop_loadstudy( 'key', 'val', ...); % no pop-up, uses the input
%                                      % parameters to select the STUDY set to load   
%                                              
%
% Optional inputs:
%   'filename' - [string] filename of the STUDY set file to load.
%   'filepath' - [string] filepath of the STUDY set file to load.
%
% Outputs:
% STUDY - the requested STUDY set structure.
% ALLEEG - the corresponding ALLEEG structure containing the 
%           datasets in the STUDY.    
%
% see also: load_ALLEEG, pop_savestudy, pop_createstudy
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


function [STUDY, ALLEEG] = pop_loadstudy(varargin)

if isempty(varargin)
    [filename, filepath] = uigetfile2('*.study', 'Load a STUDY -- pop_loadstudy()'); 
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
        eval(['load ' STUDYfile ' -mat']);
    catch
        error(['pop_loadstudy: the Study file ' STUDYfile ' could not be loaded, check file name and file path']);
    end
else
    error(['pop_loadstudy: No Study file to load was provided.']);
end
  
ALLEEG = load_ALLEEG(STUDY);

    
       