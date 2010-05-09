% pop_readsegegi() - load a segmented EGI EEG file. Pop up query 
%                    window if no arguments.
% Usage:
%   >> EEG = pop_readsegegi;             % a window pops up
%   >> EEG = pop_readsegegi( filename ); % no pop-up window
%
% Inputs:
%   filename       - first EGI file name
% 
% Outputs:
%   EEG            - EEGLAB data structure
%
% Author: Arnaud Delorme, CNL / Salk Institute, 10 April 2003
%
% See also: eeglab(), readegi(), readegihdr()

% Copyright (C) 12 Nov 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [EEG, command] = pop_readegi(filename); 
    
EEG = [];
command = '';
if nargin < 1 
	% ask user
	[filename, filepath] = uigetfile('*.RAW;*.raw', 'Choose first EGI RAW file -- pop_readsegegi()'); 
    drawnow;
	if filename == 0 return; end;
	filename = fullfile(filepath, filename);
end;

% load datas
% ----------
EEG = eeg_emptyset;

tailname = filename(end-3:end);
basename = filename(1:end-7);
index = 1;
cont = 1;
Eventdata = [];

disp('Removing trailing character of selected file to find base file name');
fprintf('Base file name is: %s\n', basename);
orifilename = [ basename sprintf('%3.3d', index) tailname ];
if ~exist(orifilename)
    disp ([ 'First file of series ''' orifilename ''' not found' ] );
    error([ 'First file of series ''' orifilename ''' not found' ] );
end;

while cont
    tmpfilename = [ basename sprintf('%3.3d', index) tailname ];
    try,
        disp(['Importing ' tmpfilename ]);
        [Head tmpdata tmpevent] = readegi( tmpfilename );
        EEG.data  = [ EEG.data  tmpdata ];
        Eventdata = [ Eventdata tmpevent ];
        index = index + 1;
    catch,
        cont = 0;
    end;
end;

% add one channel with the event data
% -----------------------------------
if ~isempty(Eventdata) & size(Eventdata,2) == size(EEG.data,2)
    EEG.data(end+1:end+size(Eventdata,1),:) = Eventdata;
end;
EEG.comments        = [ 'Original files: ' orifilename ' to ' tmpfilename ];
EEG.filepath        = '';
EEG.setname 		= 'EGI file';
EEG.nbchan          = size(EEG.data,1);
EEG.srate           = Head.samp_rate;
EEG.trials          = Head.segments;
EEG.pnts            = Head.segsamps;
EEG.xmin            = 0; 

% importing the events
% --------------------
if ~isempty(Eventdata)
    orinbchans = EEG.nbchan;
    for index = size(Eventdata,1):-1:1
        EEG = pop_chanevent( EEG, orinbchans-size(Eventdata,1)+index, 'edge', 'leading', ...
                             'delevent', 'off', 'typename', Head.eventcode(index,:), ...
                             'nbtype', 1, 'delchan', 'on');
    end;
end;

% importing channel locations
% ---------------------------
if all(EEG.data(end,1:10) == 0)
    disp('Deleting empty data reference channel (reference channel location is retained)');
    EEG.data(end,:)   = [];
    EEG.nbchan        = size(EEG.data,1);
    EEG = eeg_checkset(EEG);
end;
EEG = readegilocs(EEG);

EEG = eeg_checkset(EEG);
command = sprintf('EEG = pop_readsegegi(''%s'');', filename); 

return;
