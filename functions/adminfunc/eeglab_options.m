% eeglab_options() - handle EEGLAB options. This script (not function)
%                    set the various options in the eeg_options() file.
%
% Usage:
%   eeglab_options;
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2006-

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

% load local file
% ---------------
try,    
    %clear eeg_options; % note: we instead clear this function handle in pop_editoptions()
    
    if iseeglabdeployed
        com1 = readtxtfile(fullfile(eeglabexefolder, 'eeg_optionsbackup.txt'));
        com2 = readtxtfile(fullfile(eeglabexefolder, 'eeg_options.txt'));
        eval( com1 );
        eval( com2 );
    else
        eeg_optionsbackup;
        W_MAIN = findobj('tag', 'EEGLAB');
        if ~isempty(W_MAIN)
            tmpuserdata = get(W_MAIN, 'userdata');
            tmp_opt_path = tmpuserdata{3}; % this contain the default path to the option file
            tmpp = pwd;
            
            tmpp = fileparts(which('eeg_options.m'));
            curpathconflict = 0;
            if ~strcmpi(tmpp, tmp_opt_path)
                if ~isempty(findstr(tmp_opt_path, path)), rmpath(tmp_opt_path); end;
                if strcmpi(pwd, tmpp)
                    curpathconflict = 1;
                    fprintf('Warning: conflicting eeg_options.m file in current path (ignored)\n');
                else
                    fprintf('Warning: adding path for eeg_options.m, %s\n', tmp_opt_path);
                end;
                addpath(tmp_opt_path);
                clear functions;
            end;
            if curpathconflict
                cd(tmp_opt_path);
                eeg_options;
                cd(tmpp);
            else
                eeg_options;
            end;
        else
            eeg_options;
        end;
    end;
    
    option_savematlab = ~option_savetwofiles;
    
catch 
    lasterr
    disp('Warning: could not access the local eeg_options file');
end;
