% Plugin support function to remove plugin

% Copyright (C) 2012- Arnaud Delorme
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

function plugin_remove(foldername)

    % get plugin path
    % ---------------
    fullpluginfolder = fullfile(fileparts(which('eeglab.m')), 'plugins', foldername);
    if ~exist(fullpluginfolder)
        error([ 'Could not find folder ' foldername ' in plugins folder' ]);
    end

    disp([ 'Removing plugin folder ' foldername ]);
    try
        rmdir(fullpluginfolder, 's');
    catch
        eeglab_error;
    end
