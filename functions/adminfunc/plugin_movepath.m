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
%        warns        -[0,1] Allow display [1] or do not display [0] warnings.
%                      Warnings are restored at the end of the process. Default[0] 
%  Outputs:
%        oldpath  - Original MATLAB path before entering this function
%        newpath  - MATLAB path after being modified by this function
% 
% Author: Ramon Martinez-Cancino and Arnaud Delorme, SCCN, 2016
%         
% Copyright (C) 2016  Ramon Martinez-Cancino,INC, SCCN
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function[oldpath, newpath] = plugin_movepath(foldername,pluginpos, varargin)
oldpath = []; newpath = [];

try
    options = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else g = []; end
catch
    disp('plugin_movepath() error: calling convention {''key'', value, ... } error'); return;
end

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
end

% Retrieving path
comp = computer;
if strcmpi(comp(1:2), 'PC')
    newpathtest = [ pluginfolder ';' ];
else
    newpathtest = [ pluginfolder ':' ];
end
ind = strfind(oldpath, newpathtest);

% Checking out if the work is already done
if strcmp(pluginpos,'begin') && ind == 1, return; end
if strcmp(pluginpos,'end')   && strcmp(oldpath(end-length(pluginfolder):end),pluginfolder), return; end
    
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
 
 newpath = path; % Retrieving new path

 % required here because path not added yet
% to the admin folder
function res = ismatlab;

v = version;
if v(1) > '4'
    res = 1;
else
    res = 0;
end
