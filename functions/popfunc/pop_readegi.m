% pop_readegi() - load a EGI EEG file (pop out window if no arguments).
%
% Usage:
%   >> EEG = pop_readegi;             % a window pops up
%   >> EEG = pop_readegi( filename );
%
% Inputs:
%   filename       - EGI file name
%   datachunks     - desired frame numbers (see readegi() help)
%                    option available from the command line only
% 
% Outputs:
%   EEG            - EEGLAB data structure
%
% Author: Arnaud Delorme, CNL / Salk Institute, 12 Nov 2002
%
% See also: eeglab(), readegi(), readegihdr()

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.12  2004/11/10 02:43:06  arno
% add new pop up gui for epoch file
%
% Revision 1.11  2004/11/10 02:10:40  arno
% reading chunks of data from the command line
%
% Revision 1.10  2003/09/22 23:42:19  arno
% debuging urevent
%
% Revision 1.9  2003/07/11 21:44:00  arno
% removing warning for urevent
%
% Revision 1.8  2003/04/10 17:58:36  arno
% filter for file read
%
% Revision 1.7  2002/12/06 03:18:39  arno
% same
%
% Revision 1.6  2002/12/06 03:07:19  arno
% debuging channel import
%
% Revision 1.5  2002/12/06 02:50:05  arno
% use leading edge
%
% Revision 1.4  2002/12/06 02:44:39  arno
% adding event import
%
% Revision 1.3  2002/12/05 02:50:25  arno
% debugging event reading
%
% Revision 1.2  2002/11/14 23:35:36  arno
% header
%
% Revision 1.1  2002/11/13 02:34:22  arno
% Initial revision
%

function [EEG, command] = pop_readegi(filename, datachunks); 
    
EEG = [];
command = '';
disp('Warning: this function can only import continuous files or');
disp('         epoch files with only one length for data epochs');
if nargin < 1 
	% ask user
	[filename, filepath] = uigetfile('*.RAW;*.raw', 'Choose an EGI RAW file -- pop_readegi()'); 
    drawnow;
	if filename == 0 return; end;
	filename = [filepath filename];

    fid = fopen(filename, 'rb', 'b');
    if fid == -1, error('Cannot open file'); end;
    head = readegihdr(fid);
    fclose(fid);
    
    if head.segments ~= 0
        promptstr    = { sprintf('Segment/frame number (default:1:%d)', head.segments) };
        inistr       = { '' };
        result       = inputdlg2( promptstr, 'Import EGI file -- pop_readegi()', 1,  inistr, 'pop_readegi');
        if length(result) == 0 return; end;
        datachunks   = eval( [ '['  result{1} ']' ] );
    else
        datachunks   = [];
        disp('Only one segment, cannot read portion of the file');
    end;
end;

% load datas
% ----------
EEG = eeg_emptyset;
if exist('datachunks')
    [Head EEG.data Eventdata] = readegi( filename, datachunks);
else
    [Head EEG.data Eventdata] = readegi( filename);
end;
if ~isempty(Eventdata) & size(Eventdata,2) == size(EEG.data,2)
    EEG.data(end+1:end+size(Eventdata,1),:) = Eventdata;
end;
EEG.filename        = filename;
EEG.filepath        = '';
EEG.setname 		= 'EGI file';
EEG.nbchan          = size(EEG.data,1);
EEG.srate           = Head.samp_rate;
EEG.trials          = Head.segments;
EEG.pnts            = Head.segsamps;
EEG.xmin            = 0; 

% importing the events
% --------------------
EEG = eeg_checkset(EEG);
if ~isempty(Eventdata)
    orinbchans = EEG.nbchan;
    for index = size(Eventdata,1):-1:1
        EEG = pop_chanevent( EEG, orinbchans-size(Eventdata,1)+index, 'edge', 'leading', ...
                             'delevent', 'off', 'typename', Head.eventcode(index,:), ...
                             'nbtype', 1, 'delchan', 'on');
         Head.eventcode(end,:) = [];
    end;

    % renaming event codes
    % --------------------
    try,
        alltypes = { EEG.event.type };
        if isstr(alltypes{1})
            indepoc = strmatch('epoc', lower(alltypes), 'exact');
            indtim  = strmatch('tim0', lower(alltypes), 'exact');
        
            % if epoc but no tim0 then epoc represent pauses in recording
            if isempty(indtim) & ~isempty(indepoc)
                for index = indepoc
                    EEG.event(index).type = 'boundary';
                end;
            end;
            % other wise if both non-empty data epochs
            if ~isempty(indtim) & ~isempty(indepoc)
                if rem(size(EEG.data,2) / (length(indepoc)+1),1) == 0
                    EEG.event(index) = []; % remove epoch events
                    EEG.trials       = length(indepoc)+1;
                else
                    disp('Warning: data epochs detected but wrong data size');
                end;
            end; 
        end;
     catch, disp('Warning: event renaming failed'); end;
end;


EEG = eeg_checkset(EEG, 'makeur');
EEG = eeg_checkset(EEG, 'eventconsistency');
if nargin < 1 
    command = sprintf('EEG = pop_readegi(''%s'', %s);', filename, vararg2str({datachunks}) ); 
end;

return;
