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
    uilist   = { { 'style' 'text' 'string' 'Time interval in seconds (i.e. [0 100]; default all):' } ...
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
%EEG = eeg_checkset(EEG);

if ~isempty(findstr('triggerfile', lower(options)))
    if strcmp(r.triggerfile,'on')
        [datdir,name,ext]=fileparts(fullFileName);
        tfilename=fullfile(datdir,[name '.trg']);
        tfexist=exist(tfilename);
        if tfexist > 0
            disp(['Loading file ' tfilename ' ...']);
            EEG = pop_importevent( EEG,  'append', 'no', 'event',tfilename, 'fields',{'latency', 'byte', 'type'},...
                'skipline',1, 'timeunit',1, 'align',NaN);
            EEG = pop_editeventfield( EEG,'byte', [], 'init_index', [], 'init_time', []);
        else
            disp(['ERROR Trigger file: ' tfilename ' does not exist !!!!'])
        end
        
        % make latencies integer (as they should be)
        % ------------------------------------------
        for i=1:length(EEG.event)
            EEG.event(i).latency = round(EEG.event(i).latency);
        end;
    end
end

%EEG = eeg_checkset(EEG);

if length(options) > 2
    command = sprintf('EEG = pop_loadeep(''%s'' %s);',fullFileName, options); 
else
    command = sprintf('EEG = pop_loadeep(''%s'');',fullFileName); 
end;
return;

