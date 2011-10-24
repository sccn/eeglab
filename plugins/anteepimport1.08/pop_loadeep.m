% pop_loadeep() - Load an EEProbe continuous file (*.cnt).
%                 (pop out window if no arguments)
%
% Usage:
%   >> [EEG] = pop_loadeep;
%   >> [EEG] = pop_loadeep( filename, 'key', 'val', ...);
%
% Graphic interface:
%
%   "Time interval in seconds" - [edit box] specify time interval [min max]
%                                to import portion of data. Command line equivalent
%                                in loadeep: 'time1' and 'time2'
%   "Import triggers "         - [checkbox] set this option to import triggers from the
%                                trigger file (*.trg). Command line equivalent 'triggerfile'.
% Inputs:
%   filename                   - file name
%
% Optional inputs:
%   'triggerfile'               -'on' or 'off' (default = 'off')
%   Same as loadeep() function.
%
% Outputs:
%   [EEG]                       - EEGLAB data structure
%
% Note:
% This script is based on pop_loadcnt.m to make it compatible and easy to use in
% EEGLab.
%
% Author: Maarten-Jan Hoeve, ANT Software, Enschede, The Netherlands, 8 October 2003
%
% See also: eeglab(), loadeep()
%

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2003 Maarten-Jan Hoeve, m.hoeve@ieee.org
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

% $Log: pop_loadeep.m,v $
% Revision 1.6  2008-06-20 10:36:05  mvelde
% fixed reading of triggers for epoch selection
%
% Revision 1.5  2006-09-25 14:04:03  mvelde
% updated for EEGLAB 5.03
%
% Revision 1.4  2005/07/21 10:12:51  mvelde
% updated information
%
% Revision 1.3  2005/06/24 13:46:46  mvelde
% fixed typo in matlab syntax
%
% Revision 1.2  2005/06/08 08:16:37  mvelde
% converted files to unix format
%
% Revision 1.1  2004/11/26 13:22:22  jwiskerke
% Added eeglab to cvs.
%
% Revision 1.3  2003/10/24 13:34:41  Maarten-Jan Hoeve
% Added GNU Licence and updated revision history
%
% Revision 1.1.1.2  2003/10/17 09:55:20  mvelde
% updated: consistent copyrights, arguments/data labels, fixed some typos
%
% Advanced Neuro Technology (ANT) BV, The Netherlands, www.ant-neuro.com / info@ant-neuro.com
%

function [EEG, command]=pop_loadeep(filename, varargin);

command = '';
filepath = '';
EEG=[];

if nargin < 1

	% ask user
	[filename, filepath] = uigetfile('*.CNT;*.cnt', 'Choose an EEProbe continuous file -- pop_loadeep()');
    drawnow;
	if filename == 0 return; end;

	% popup window parameters
	% -----------------------
    uigeom     = { [1 0.5] [1.09 0.13 0.4]};
    uilist   = { { 'style' 'text' 'string' 'Time interval in s (i.e. [0 100]; not compatible with importing triggers):' } ...
                 { 'style' 'edit' 'string' '' } ...
                 { 'style' 'text' 'string' 'Check to import triggers from EEProbe trigger file (*.trg)' } ...
                 { 'style' 'checkbox' 'string' '' } {} };

	result = inputgui(uigeom, uilist, 'pophelp(''pop_loadeep'')', 'Load an EEProbe dataset');
	if length( result ) == 0 return; end;

	% decode parameters
	% -----------------
    options = [];
    if ~isempty(result{1}),
        timer =  eval( [ '[' result{1} ']' ]);
        options = [ options ', ''time1'', ' num2str(timer(1)) ', ''time2'', ' num2str(timer(2)) ];
    end;
    if result{2}, options = [ options ', ''triggerfile'', ''on''' ]; end;
else
	options = vararg2str(varargin);
end;

% load datas
% ----------
EEG = eeg_emptyset;
if exist('filepath')
	fullFileName = sprintf('%s%s', filepath, filename);
else
	fullFileName = filename;
end;
if nargin > 0
	if ~isempty(varargin)
		r = loadeep( fullFileName, varargin{:});
	else
		r = loadeep( fullFileName);
	end;
else
	eval( [ 'r = loadeep( fullFileName ' options ');' ]);
end;

EEG.data            = r.dat;
EEG.comments        = [ 'Original file: ' fullfile(filepath, filename) ];
EEG.setname         = 'EEProbe continuous data';
EEG.nbchan          = r.nchannels;
EEG.xmin            = (r.sample1-1)/r.rate;
EEG.srate           = r.rate;
EEG.pnts            = r.nsmpl;
EEG.chanlocs        = r.chanlocs;
EEG = eeg_checkset(EEG);

% import events
if ~isempty(findstr('triggerfile', lower(options)))
    if strcmp(r.triggerfile,'on')
        [datdir,name,ext]=fileparts(fullFileName);
        tfilename=fullfile(datdir,[name '.trg']);
        tfexist=exist(tfilename);
        if tfexist > 0
            disp(['Loading file ' tfilename ' ...']);
            % read all triggers
            trg = read_eep_trg(tfilename);
            % -------------
            % add only triggers in loaded epoch
            j = 1;
            for i = 1:length(trg)
                % adjust latency (#samples) for srate and offset 'timer(1)'
                trg_latency = ( (trg(i).time/1000.0) - EEG.xmin ) * EEG.srate + 1;
                if trg_latency >= 0.5 && trg_latency < EEG.pnts*EEG.trials
                    EEG.event(j).type = trg(i).code;
                    EEG.event(j).latency = trg_latency;
                    j = j + 1;
                end;
            end;
            if j < length(trg)
                fprintf('pop_loadeep warning: %d/%d events had out-of-bounds latencies and were removed\n', ...
                    length(trg) - length(EEG.event), length(trg));
            end;
        else
            disp(['ERROR Trigger file: ' tfilename ' does not exist !!!!'])
        end;
    end;
end;

EEG = eeg_checkset(EEG);

if length(options) > 2
    command = sprintf('EEG = pop_loadeep(''%s'' %s);',fullFileName, options);
else
    command = sprintf('EEG = pop_loadeep(''%s'');',fullFileName);
end;
return;

