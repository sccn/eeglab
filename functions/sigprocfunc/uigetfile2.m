% uigetfile2() - same as uigetfile but remember folder location.
%
% Usage: >> uigetfile2(...)
%
% Inputs: Same as uigetfile
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
% Revision 1.5  2004/12/01 23:18:14  hilit
% now changed correctly
%
% Revision 1.4  2004/12/01 23:15:31  hilit
% file can be now read from Matlab6 even if created in Matlab 7
%
% Revision 1.3  2004/11/29 23:12:43  hilit
% verifing that everyone have perimissions to eeglab.cfg
%
% Revision 1.2  2004/11/23 19:27:14  scott
% edit help
%
% Revision 1.1  2004/09/09 22:56:44  arno
% Initial revision
%

function varargout = uigetfile2(varargin);
    
    if nargin < 1
        help uigetfile2;
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
    [varargout{1} varargout{2}] = uigetfile(varargin{:});
    if varargout{1} ~= 0
        Path = varargout{2};
        cd(olddir);
        try, save(fullfile(getenv('TEMP'),'eeglab.cfg'),'Path','-mat','-V6'); % Matlab 7
        catch, 
            try,  save(fullfile(getenv('TEMP'),'eeglab.cfg'),'Path','-mat');
          catch, error('uigetfile2: save error, out of space or file permission problem');
            end
        end
        if isunix
            system('chmod 777 eeglab.cfg');
        end
    else
       cd(olddir) 
    end;
     
