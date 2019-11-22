% pop_writeeeg - write EEGLAB dataset to disk in EDF/GDF or BDF format
%
% pop_writeeeg( EEG ) % pops up a window
% pop_writeeeg( EEG, filename, 'key', 'val' )
%
% Inputs:
%  EEG        - EEGLAB dataset
%  filename   - [string] filename
%
% Optional keys (same as writeeeg):
%  'TYPE'         - ['GDF'|'EDF'|'BDF'|'CFWB'|'CNT'] file format for writing
%                   default is 'EDF'.
%  See writeeeg for more information
%
% Author: Arnaud Delorme, SCCN, UCSD/CERCO, 2009
%         Based on BIOSIG, sopen and swrite

% Copyright (C) 22 March 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [command] = pop_writeeeg(EEG, filename, varargin)
command = '';

% enforce use the str2double of Biosig
if ~plugin_askinstall('Biosig', 'sopen'), return; end
biosigpathfirst

if nargin < 2
    if EEG.trials > 1
        res = questdlg2( [ 'This dataset contains data epochs.' 10 'Do you want to export the concatenated' 10 'data epochs?' ], '', 'No', 'Yes', 'Yes');
        if strcmpi(res, 'No')
            return;
        end
    end
    
	% ask user
	[filename, filepath] = uiputfile('*.*', 'Enter a file name -- pop_writeeeg()'); 
	if filename == 0, return; end
	filename = fullfile(filepath,filename);
    
    % file format
    % -----------
    fileformats = { 'EDF' 'GDF' 'BDF' };
    uilist = { { 'style' 'text' 'String' 'File format' } ...
               { 'style' 'listbox' 'string' strvcat(fileformats) 'value' 1 } };
    geom = [1 1];
    result = inputgui( 'geometry', geom, 'uilist', uilist, 'helpcom', 'pophelp(''pop_writeeeg'')', ...
                     'title', 'Write data using BIOSIG -- pop_writeeeg()', 'geomvert', [1 2.5]);
    if isempty(result), return; end

    if result{1} == 3
        disp('WARNING: there is a potential issue BDF file header, see https://sccn.ucsd.edu/bugzilla/show_bug.cgi?id=1020');
    end
    
    options = { 'TYPE' fileformats{result{1}} };
else
    options = varargin;
end

warning('off', 'MATLAB:intConvertNonIntVal');
if ~isempty(EEG.chanlocs)
    tmpchanlocs = EEG.chanlocs;
    writeeeg(filename, EEG.data(:,:), EEG.srate, 'label', { tmpchanlocs.labels }, 'EVENT', EEG.event, options{:});
else
    writeeeg(filename, EEG.data(:,:), EEG.srate, 'EVENT', EEG.event, options{:});
end
warning('on', 'MATLAB:intConvertNonIntVal');

command = sprintf('pop_writeeeg(EEG, ''%s'', %s);', filename, vararg2str(options)); 
biosigpathlast

