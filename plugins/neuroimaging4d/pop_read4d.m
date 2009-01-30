% pop_read4d() - load a 4d MEG file (pop out window if no arguments).
%
% Usage:
%   >> EEG = pop_read4d;             % a window pops up
%   >> EEG = pop_read4d( filename );
%
% Inputs:
%   filename       - 4D set file name
% 
% Outputs:
%   EEG            - EEGLAB data structure
%
% Author: Christian Wienbruch / Konstanz 04 2004
%
% See also: eeglab(), read4d(), read4dhdr()

% $Log: not supported by cvs2svn $
% Revision 1.1  2005/06/08 13:41:21  chris
% Initial revision
%
% Revision 1.1  2004/11/15 19:36:53  chris
% Initial revision
%

function [EEG, command] = pop_read4d(filename); 
    
EEG = [];
command = '';
if nargin < 1 
	% ask user
	[setname, filepath] = uigetfile('*.m4d', 'Choose an 4D pdf file -- pop_read4d()'); 
    drawnow;
	if setname == 0 return; end;
	filename = [filepath setname];
end;

% load datas
% ----------
EEG = eeg_emptyset;

[Head EEG.data Eventdata] = read4d( filename );
if ~isempty(Eventdata) & size(Eventdata,2) == size(EEG.data,2)
    EEG.data(end+1:end+size(Eventdata,1),:) = Eventdata;
end;

EEG.comments        = [ 'Original file: ' filename ];
EEG.setname         = '4D file';
EEG.nbchan          = size(EEG.data,1);
EEG.srate           = Head.samp_rate;
EEG.trials          = Head.segments;
EEG.pnts            = Head.segsamps;
EEG.xmin            = Head.FirstLatency;

EEG = eeg_checkset(EEG);
command = sprintf('EEG = pop_read4d(''%s'');', filename); 

return;
