% uiputfile2() - same as uigputfile but remember folder location.
%
% Usage: >> uiputfile2(...)
%
% Inputs: Same as uiputfile
%
% Author: Arnaud Delorme & Hilit Serby, Scott Makeig, SCCN, UCSD, 2004
%         Thanks to input from Bas Kortmann
%
% Copyright (C) Arnaud Delorme & Hilit Serby, Scott Makeig, SCCN, UCSD, 2004

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


function varargout = uiputfile2(varargin);
    
    if nargin < 1
        help uiputfile2;
        return;
    end;
    
    % remember old folder
    %% search for the (mat) file which contains the latest used directory 
    % -------------------
    olddir = pwd;
    try,
        eeglab_options;
        if option_rememberfolder
            tmp_fld = getenv('TEMP');
            if isempty(tmp_fld) & isunix
                if exist('/tmp') == 7
                    tmp_fld = '/tmp';
                end;
            end;
            if exist(fullfile(tmp_fld,'eeglab.cfg'))
                load(fullfile(tmp_fld,'eeglab.cfg'),'Path','-mat');
                s = ['cd([''' Path '''])'];
                if exist(Path) == 7, eval(s); end;
            end;
        end;
    catch, end;

    %% Show the open dialog and save the latest directory to the file
    % ---------------------------------------------------------------
    [varargout{1} varargout{2}] = uiputfile(varargin{:});
    try,
        if option_rememberfolder
            if varargout{1} ~= 0
                Path = varargout{2};
                try, save(fullfile(tmp_fld,'eeglab.cfg'),'Path','-mat','-V6'); % Matlab 7
                catch, 
                    try,  save(fullfile(tmp_fld,'eeglab.cfg'),'Path','-mat');
                    catch, error('uigetfile2: save error, out of space or file permission problem');
                    end
                end
                if isunix
                    eval(['cd ' tmp_fld]);
                    system('chmod 777 eeglab.cfg');
                end
            end;
        end;
    catch, end;
    cd(olddir) 
    

    
