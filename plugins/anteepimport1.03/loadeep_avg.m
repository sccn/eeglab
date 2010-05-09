% loadeep_avg() - Load EEProbe average data (*.avr).
%
% Usage:
%   >> [EEG] = loadeep_avg(file);
%
% Inputs:
%            filename    - name of EEProbe averaged file, including extension (*.avr)
%
% Outputs:
%   [EEG]                - data structure holding EEG information and other
%                          relevant information.
%
% Author: Maarten-Jan Hoeve, ANT Software, Enschede, The Netherlands, 8 October 2003
%
% See also: eeglab(), pop_loadeep_avg()
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

function r=loadeep_avg(file)

if nargin < 1
	help loadeep;
	return;
end;	

% defaults
if ~any(file=='.'), file=[file '.avr']; end
[datdir,name,ext]=fileparts(file);
disp(['Loading file ' file ' ...'])

eeg=read_eep_avr(file);

r.filename=file;
r.totsamples=eeg.npnt;
r.nchannels=eeg.nchan;
r.rate=eeg.rate;

% Create struct for holding channel labels
for i=1:r.nchannels
    chanlocs(i).labels=char(eeg.label(i));
    chanlocs(i).theta=0;
    chanlocs(i).radius=0;
    chanlocs(i).X=0;
    chanlocs(i).Y=0;
    chanlocs(i).Z=0;
    chanlocs(i).sph_theta=0;
    chanlocs(i).sph_phi=0;
    chanlocs(i).sph_radius=0;    
end
r.chanlocs=chanlocs;

r.nsmpl=eeg.npnt;
r.time=eeg.time;
r.dat=eeg.data;
r.nsweeps=eeg.nsweeps;     %number of trials averaged
r.xmin=eeg.xmin; 
r.xmax=eeg.xmax;
r.variance=eeg.variance;   % variance (nchan x npnt)
r.condlab=eeg.condlab;     % string with condition label
r.condcol=eeg.condcol;     % string with color code for condition
r.trialc=eeg.trialc;       % total number of trial in original data
r.rejtrialc=eeg.rejtrialc; % number of rejected trials
