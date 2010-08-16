% pop_loadeep_avg() - Load an EEProbe average file (*.avr). 
%                     (pop out window if no arguments)
%
% Usage:
%   >> [EEG] = pop_loadeep_avg;
%   >> [EEG] = pop_loadeep_avg( filename);
%
% Inputs:
%   filename                   - file name
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
% See also: eeglab(), loadeep_avg()
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

% $Log: pop_loadeep_avg.m,v $
% Revision 1.4  2006-09-25 14:04:03  mvelde
% updated for EEGLAB 5.03
%
% Revision 1.3  2005/07/21 10:12:51  mvelde
% updated information
%
% Revision 1.2  2005/06/08 08:16:37  mvelde
% converted files to unix format
%
% Revision 1.1  2004/11/26 13:22:22  jwiskerke
% Added eeglab to cvs.
%
% Revision 1.2  2003/10/24 13:34:41  Maarten-Jan Hoeve
% Added GNU Licence and updated revision history
%
% Revision 1.1.1.2  2003/10/17 09:55:20  mvelde
% updated: consistent copyrights, arguments/data labels, fixed some typos
%
% Advanced Neuro Technology (ANT) BV, The Netherlands, www.ant-neuro.com / info@ant-neuro.com
%

function [EEG, command]=pop_loadeep_avg(filename); 

command = '';
EEG=[];

if nargin < 1 

	% ask user
	[filename, filepath] = uigetfile('*.AVR;*.avr', 'Choose an EEProbe average file -- pop_loadeep_avg()'); 
    drawnow;
	if filename == 0 return; end;
end;

% load datas
% ----------
EEG = eeg_emptyset;
if exist('filepath')
	fullFileName = sprintf('%s%s', filepath, filename);
else
	fullFileName = filename;
end;	

r = loadeep_avg(fullFileName);

EEG.data            = r.dat;
EEG.comments        = [ 'Original file: ' fullfile(filename, filepath) ];
EEG.setname 		= 'Averaged ANT EEP file';
EEG.nbchan          = r.nchannels; 
%EEG.xmin            = (r.xmin-1)/r.rate;
%EEG.xmax           = (r.xmax-1)/r.rate;
EEG.srate           = r.rate;
EEG.pnts            = r.nsmpl;
EEG.chanlocs        = r.chanlocs;
EEG = eeg_checkset(EEG);

command = sprintf('EEG = pop_loadeep_avg(''%s'');',fullFileName); 

return;
