% uigetfile2() - same as uigetfile but remember folder location.
%
% Usage: >> uigetfile2(...)
%
% Inputs: Same as uigetfile
%
% Author: Arnaud Delorme & Hilit Serby, Scott Makeig, SCCN, UCSD, 2004
%         thanks to inputs from Bas Kortmann
%
% Copyright (C) Arnaud Delorme & Hilit Serby, Scott Makeig, SCCN, UCSD, 2004
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

function varargout = uigetfile2(varargin);
    
    if nargin < 1
        help uigetfile2;
        return;
    end;
    
    % remember old folder
    %% search for the file which contains the latest used directory (== mat file)
    % -------------------
    olddir = pwd;
    if exist(fullfile(getenv('TEMP'),'eeglab.cfg'))
        load(fullfile(getenv('TEMP'),'eeglab.cfg'),'Path','-mat');
        s = ['cd([''' Path '''])'];
        if exist(Path) == 7, eval(s); end;
    end;

    %% show the open dialog and save the latest directory to the file
    % ---------------------------------------------------------------
    [varargout{1} varargout{2}] = uigetfile(varargin{:});
    if varargout{1} ~= 0
        Path = varargout{2};
        save(fullfile(getenv('TEMP'),'eeglab.cfg'),'Path','-mat');
    end;
    cd(olddir)
    