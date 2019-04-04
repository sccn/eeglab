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

function [EEG, command] = pop_readegi(filename); 
    
EEG = [];
command = '';
if nargin < 1 
	% ask user
	[filename, filepath] = uigetfile('*.RAW;*.raw', 'Choose first EGI RAW file -- pop_readsegegi()'); 
    drawnow;
	if filename == 0 return; end
	filename = fullfile(filepath, filename);
end

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
end

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
    end
end

% add one channel with the event data
% -----------------------------------
if ~isempty(Eventdata) && size(Eventdata,2) == size(EEG.data,2)
    EEG.data(end+1:end+size(Eventdata,1),:) = Eventdata;
end
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
    end
end

% importing channel locations
% ---------------------------
if all(EEG.data(end,1:10) == 0)
    disp('Deleting empty data reference channel (reference channel location is retained)');
    EEG.data(end,:)   = [];
    EEG.nbchan        = size(EEG.data,1);
    EEG = eeg_checkset(EEG);
end
EEG = readegilocs(EEG);

EEG = eeg_checkset(EEG);
command = sprintf('EEG = pop_readsegegi(''%s'');', filename); 

return;
