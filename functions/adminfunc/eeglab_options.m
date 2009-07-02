% eeglab_options() - handle EEGLAB options. This script (not function)
%                    set the various options in the eeg_options() file.
%
% Usage:
%   eeglab_options;
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2006-

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, 2006-
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
% Revision 1.3  2008/04/16 17:43:42  arno
% nothing
%
% Revision 1.2  2006/02/03 18:53:30  arno
% more detailed message
%
% Revision 1.1  2006/01/31 20:07:08  arno
% Initial revision
%

% load local file
% ---------------
try,
    clear eeg_options;
    
    W_MAIN = findobj('tag', 'EEGLAB');
    if ~isempty(W_MAIN)
        tmpuserdata = get(W_MAIN, 'userdata');
        tmp_opt_path = tmpuserdata{3}; % this contain the default path to the option file
    else
        tmp_opt_path = '';
    end;
    
    tmp_opt_path2 = which('eeg_options');
    tmp_opt_path2 = fileparts( tmp_opt_path2 );
    if ~isempty(tmp_opt_path) & ~strcmpi(tmp_opt_path2, tmp_opt_path)
        if exist(fullfile(tmp_opt_path, 'eeg_options.m')) == 2
            fprintf('Warning: you should delete the eeg_option.m file in folder %s\n', tmp_opt_path2);
            fprintf('         using instead the eeg_option.m file accessible at EEGLAB startup in %s\n', tmp_opt_path);
            fprintf('         To use the first option file, restart EEGLAB and reload datasets\n');
            addpath(tmp_opt_path);
        else
            fprintf('IMPORTANT WARNING: the current eeg_option.m file in %s HAS BEEN DELETED\n', tmp_opt_path);
            fprintf('                   EEGLAB SHOULD BE RESTARTED TO ENSURE STABILITY\n', tmp_opt_path);
        end;
    end;
    
    eeg_optionsbackup;
    eeg_options;
    
    option_savematlab = ~option_savetwofiles;
    
catch 
    disp('Warning: could not access the local eeg_options file');
end;
