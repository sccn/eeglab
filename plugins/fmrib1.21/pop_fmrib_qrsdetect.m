% pop_fmrib_qrsdetect() -  GUI for fmrib_qrsdetect:
%   Detect QRS events from ECG channel and add them to
%   EEGLAB data event structure.
%
% Usage:
%    >> [EEGOUT, COM] = pop_fmrib_qrsdetec(EEG) % for gui
% or >> [EEGOUT, COM] = pop_fmrib_qrsdetec(EEG,ecgchan,'ename','yes') 
% for command line. 'no' or 'yes' can be used in the last argument:
%
% ecgchan: is the ecg channel number
% 'ename': event name to call the found peaks
% 'yes' or 'no': do you want to delete ecg channel?
%
% Inputs: EEG data structure
%
% GUI:
%   Enter ECG Channel Number: ECG Data channel to extract QRS peaks from.
%   What would you like to name the QRS events?: Event type to assign to 
%       the QRS peaks detected.
%   Delete ECG channel when finished (yes/no)?: I think this is obvious!
%  
%
% Author: Rami Niazy, FMRIB Centre, University of Oxford.
%
%
%
% Copyright (c) 2004 University of Oxford.
%


% Copyright (C) 2004 University of Oxford
% Author:   Rami K. Niazy, FMRIB Centre
%           rami@fmrib.ox.ac.uk
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

% Last update:
% Oct 26, 2004

function [EEGOUT, command] = pop_fmrib_qrsdetect(EEG,varargin); 

EEGOUT = [];
command = '';

% check inputs
%--------------------
nargchk(1,4,nargin);

if nargin>1 & nargin<4
    error('incorrect input arguments','Please see help.');
end

if isempty(EEG.data)
    errordlg2(strvcat('Error: no data'), 'Error');
    error('pop_fmrib_qrsdetect error: no data','fmrib_qrsdetect() error!'); return;
end


% check for DSP Toolbox
%----------------------
if ~exist('firls')
   error('QRS event detection requires the signal processing toolbox.','fmrib_qrsdetect() error!');
end

%Launch user input GUI and get inputs
%------------------------------------
if nargin ==1
	prompt={'Enter ECG channel number:','What would you like to name the QRS events?','Delete ECG channel when finished (yes/no)?'};
	def={'','qrs','no'};
	dlgtitle='pop_fmrib_qrsdetect()';
	helpfunc='fmrib_qrsdetect';
	answer=inputdlg2(prompt,dlgtitle,1,def,helpfunc);
    if ~isempty(answer)
		ecgchan=str2num(char(answer(1)));
		qrsevent=char(answer(2));
		ecgdel=char(answer(3));
    else
        return;
    end
else
    ecgchan=varargin{1};
    qrsevent=varargin{2};
    ecgdel=varargin{3};
end

%Detect QRS
%----------
fprintf('\nFinding QRS Peaks...\n');
Peaks=fmrib_qrsdetect(EEG,ecgchan);

%Add Peaks to events
%-------------------
fprintf('Writing QRS events to event structure...\n');
for index = 1:length(Peaks)
	EEG.event(end+1).type  = qrsevent;
	EEG.event(end).latency = Peaks(index);
    if EEG.trials > 1 | isfield(EEG.event, 'epoch');
        EEG.event(end).epoch = 1+floor((EEG.event(end).latency-1) / EEG.pnts);
    end
end

if EEG.trials > 1
    EEG = pop_editeventvals( EEG, 'sort', {  'epoch' 0 'latency', [0] } );
else
    EEG = pop_editeventvals( EEG, 'sort', { 'latency', [0] } );
end

if isfield(EEG.event, 'urevent')
    EEG.event = rmfield(EEG.event, 'urevent');
end

EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = eeg_checkset(EEG, 'makeur');

% delete ECG channel
% ---------------
if strcmp(lower(ecgdel), 'yes')
	EEG = pop_select(EEG, 'nochannel', ecgchan);
end

%return command
%--------------
EEGOUT=eeg_checkset(EEG);
command = sprintf('EEG=pop_fmrib_qrsdetect(EEG,%d,''%s'',''%s'')',ecgchan,qrsevent,ecgdel);
return;
