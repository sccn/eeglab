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
homefolder = '';
try,    
    %clear eeg_options; % note: we instead clear this function handle in pop_editoptions()
    
    if iseeglabdeployed
        com1 = readtxtfile(fullfile(eeglabexefolder, 'eeg_optionsbackup.txt'));
        com2 = readtxtfile(fullfile(eeglabexefolder, 'eeg_options.txt'));
        eval( com1 );
        eval( com2 );
    else
        eeg_optionsbackup;
        icadefs;
        
        % folder for eeg_options file (also update the pop_editoptions)
        if ~isempty(EEGOPTION_PATH)
             homefolder = EEGOPTION_PATH;
        elseif ispc
             if ~exist('evalc'), eval('evalc = @(x)(eval(x));'); end;
             homefolder = deblank(evalc('!echo %USERPROFILE%'));
        else homefolder = '~';
        end;
        
        option_file = fullfile(homefolder, 'eeg_options.m');
        oldp = pwd;
        try
            if ~isempty(dir(option_file))
                cd(homefolder);
            else
                tmpp2 = fileparts(which('eeglab_options.m'));
                cd(tmpp2);
            end;
        catch, end;
        eeg_options; % default one with EEGLAB
        cd(oldp);
    end;
    option_savematlab = ~option_savetwofiles;
    
catch 
    lasterr
    disp('Warning: could not access the local eeg_options file');
end;
