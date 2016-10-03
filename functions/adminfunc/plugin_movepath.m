%plugin_movepath()- Given a path to a plugin folder, this function will
%                     put the plugin at the bottom of the path.
%
% Usage:
%         plugin_movepath('x','begin');           % Put plugin 'x' at the top of the path
%         plugin_movepath('x','end');             % Put plugin 'x' at the bottom of the path
%         plugin_movepath('x','begin','warns',1); % Put plugin 'x' at the top of the path and show warnings
%                                                                  
%  Inputs:
%         foldername   - [string] Plugin name or part of it
%         pluginpos    - {'begin','end'} Position to move the plugin in the path.
%                      To the top ('begin') or to the bottom ('end').
% Optional inputs:
%        warns        -[0,1] Allow isplay [1] or do not display [0] warnings.
%                      Warnings are restored at the end of the process. Default[0] 
%  Outputs:
%        oldpath  - Original MATLAB path before entering this function
%        newpath  - MATLAB path after being modified by this function
% 
% Author: Ramon Martinez-Cancino and Arnaud Delorme, SCCN, 2016
%         
% Copyright (C) 2016  Ramon Martinez-Cancino,INC, SCCN
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

function[oldpath, newpath] = plugin_movepath(foldername,pluginpos, varargin)
oldpath = []; newpath = [];

try
    options = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else g = []; end;
catch
    disp('plugin_movepath() error: calling convention {''key'', value, ... } error'); return;
end;

try g.warns;   catch, g.warns  = 0;  end; % NO warnings by default

% Checking entries
if sum(strcmp(pluginpos,{'begin','end'})) ~= 1
    fprintf(2,'plugin_movepath error: Invalid function argument\n');
    return;
end
    
% Look in the plugin list for the plugin name provided
global PLUGINLIST
hitindx = find(~cellfun(@isempty,strfind(lower({PLUGINLIST.plugin}),lower(foldername))));
if isempty(hitindx)
    fprintf(2,'plugin_movepath error: Unidentified plugin folder\n')
    return;
else
    eeglabfolder = fileparts(which('eeglab.m'));
    pluginfolder = fullfile(eeglabfolder,'plugins',PLUGINLIST(hitindx).foldername );
end

% Backing up old path
if ismatlab
    oldpath = matlabpath;
else
    oldpath = path;
end;

% Retreiving path
comp = computer;
if strcmpi(comp(1:2), 'PC')
    newpathtest = [ pluginfolder ';' ];
else
    newpathtest = [ pluginfolder ':' ];
end;
ind = strfind(oldpath, newpathtest);

% Checking out if the work is already done
if strcmp(pluginpos,'begin') && ind == 1, return; end;
if strcmp(pluginpos,'end')   && strcmp(oldpath(end-length(pluginfolder):end),pluginfolder), return; end;
    
% Shooting down warnings
if ~g.warns
    tmpwarn = warning;
    warning off;
end

% Remove folder and subfolders from path
rmpath(genpath(pluginfolder));

% Add folder to the path again in the requested position
 if strcmp(pluginpos,'begin')
     addpath(genpath(pluginfolder),'-begin');
     fprintf(1,['EEGLAB warning: to avoid name conflict ' PLUGINLIST(hitindx).foldername ' functions relocated at the end of the path\n']);
 elseif strcmp(pluginpos,'end')
     addpath(genpath(pluginfolder),'-end');
     fprintf(1,['EEGLAB warning: to avoid name conflict ' PLUGINLIST(hitindx).foldername ' functions relocated at the end of the path\n']);
 end
 
 % Restoring warnings
 if ~g.warns
     warning(tmpwarn); 
 end
 
 newpath = path; % Retreiving new path

 % required here because path not added yet
% to the admin folder
function res = ismatlab;

v = version;
if v(1) > '4'
    res = 1;
else
    res = 0;
end;