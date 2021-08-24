echo off;

% eeglab_options() - handle EEGLAB options. This script (not function)
%                    set the various options in the eeg_options() file.
%
% Usage:
%   eeglab_options;
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2006-

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, 2006-
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

% load local file
% ---------------
homefolder = '';
try 
    %clear eeg_options; % note: we instead clear this function handle in pop_editoptions()
    
    eeg_optionsbackup;
    if isdeployed || (exist('ismcc') && ismcc)
        fileName = which('eeg_options.txt');
        
        com2 = readtxtfile(fileName);
        eval( com2 );
    else
        icadefs;
        
        % folder for eeg_options file (also update the pop_editoptions)
        if ~isempty(EEGOPTION_PATH) % in icadefs above
             homefolder = EEGOPTION_PATH;
        elseif ispc
%              if ~exist('evalc'), eval('evalc = @(x)(eval(x));'); end
%              homefolder = deblank(evalc('!echo %USERPROFILE%'));
            homefolder = getenv('USERPROFILE');
        else homefolder = '~';
        end
        
        option_file = fullfile(homefolder, 'eeg_options.m');
        oldp = pwd;
        try
            if ~isempty(dir(option_file))
                cd(homefolder);
            else
                tmpp2 = fileparts(which('eeglab_options.m'));
                cd(tmpp2);
            end
        catch, end
        echo off;
        eeg_options; % default one with EEGLAB
        cd(oldp);
    end
    option_savematlab = ~option_savetwofiles;
    
    if option_donotusetoolboxes
        disp('Not using signal processing toolbox, if you experience problem, reset your Matlab path to default')
    end
catch 
    lasterr
    disp('Warning: could not access the local eeg_options file');
end
