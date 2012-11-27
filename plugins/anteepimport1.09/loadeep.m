% loadeep() - Load EEProbe continuous data (*.cnt).
%
% Usage:
%   >> [eep] = loadeep(file, varargin);
%
% Inputs:
%            filename    - name of ANT EEP file, including extension (*.cnt)
%
%   Optional inputs:
%           'time1'      - start time in seconds.
%           'sample1'    - start at sample1, default 0.
%                          (Overrules time1) 
%           'time2'      - end time of data to be  read in seconds, default = whole file.
%           'sample2'    - end sample of data to be read, default = whole file.
%                          (Overrules time2)
% Outputs:
%   [eep]                - data structure holding continous EEG information and other
%                          relevant information.
%
% Author: Maarten-Jan Hoeve, ANT Software, Enschede, The Netherlands, 8 October 2003
%
% See also: eeglab(), pop_loadeep()
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

% $Log: loadeep.m,v $
% Revision 1.3  2005/07/21 10:12:51  mvelde
% updated information
%
% Revision 1.2  2005/06/08 08:16:37  mvelde
% converted files to unix format
%
% Revision 1.1  2004/11/26 13:22:22  jwiskerke
% Added eeglab to cvs.
%
% Revision 1.4  2003/10/24 13:34:40  Maarten-Jan Hoeve
% Added GNU Licence and updated revision history
%
% Revision 1.1.2  2003/10/17 09:55:20  mvelde
% updated: consistent copyrights, arguments/data labels, fixed some typos
%
% Advanced Neuro Technology (ANT) BV, The Netherlands, www.ant-neuro.com / info@ant-neuro.com
%

function r=loadeep(file, varargin)

if nargin < 1
	help loadeep;
	return;
end;	

% defaults
if ~any(file=='.'), file=[file '.cnt']; end
[datdir,name,ext]=fileparts(file);
disp(['Loading file ' file ' ...'])

% read short piece of data to get sampling rate, channels etc
eeg = read_eep_cnt(file,1,2);

if ~isempty(varargin)
	 r=struct(varargin{:});
else r = []; 
end;

% add defaults
try, r.time1;       catch, r.time1=0; end
try, r.sample1;     catch, r.sample1=[]; end
try, r.time2;       catch, r.time2=[]; end
try, r.sample2;     catch, r.sample2=[]; end

r.filename=file;
r.totsamples=eeg.nsample;
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

if r.sample1
    s1=r.sample1;
else
    s1=((r.time1 * r.rate)+1);
    r.sample1=s1;
end
if r.sample1 > 1
    fprintf('WARNING: event are not imported correctly when selecting a data portion\n');
end;

if r.sample2
    s2=r.sample2;
else
    if r.time2
        s2=((r.time2 * r.rate)+1);
        if (s2 > r.totsamples)|(s2 < s1); s2=r.totsamples; end
    else
        s2=r.totsamples;
    end
end
        
eeg=read_eep_cnt(file,s1,s2);
r.nsmpl=eeg.npnt;
r.time=eeg.time;
r.dat=eeg.data;
