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

% $Log: not supported by cvs2svn $
%
% Revision 1.1  2004/12/03 16:00:00  hilit
% Initial revision
%

function varargout = uiputfile2(varargin);
    
    if nargin < 1
        help uiputfile2;
        return;
    end;
    
    % remember old folder
    %% search for the (mat) file which contains the latest used directory 
    % -------------------
    olddir = pwd;
    if exist(fullfile(getenv('TEMP'),'eeglab.cfg'))
        load(fullfile(getenv('TEMP'),'eeglab.cfg'),'Path','-mat');
        s = ['cd([''' Path '''])'];
        if exist(Path) == 7, eval(s); end;
    end;

    %% Show the open dialog and save the latest directory to the file
    % ---------------------------------------------------------------
    [varargout{1} varargout{2}] = uiputfile(varargin{:});
    if varargout{1} ~= 0
        Path = varargout{2};
        cd(olddir);
        try, save(fullfile(getenv('TEMP'),'eeglab.cfg'),'Path','-mat','-V6'); % Matlab 7
        catch, 
            try,  save(fullfile(getenv('TEMP'),'eeglab.cfg'),'Path','-mat');
          catch, error('uiputfile2: save error, out of space or file permission problem');
            end
        end
        if isunix
            system('chmod 777 eeglab.cfg');
        end
    else
       cd(olddir) 
    end;
     
