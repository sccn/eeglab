function [eeg,t] = read_eep_trial(filename, triggernumber, interval);
%
% READ_EEP_TRIAL reads a data from an EEProbe *.cnt file
%
% [eeg,t] = read_eep_trial(filename, triggernumber, interval);
%
% interval = [ -2, 5] reads a window of -2 to 5 seconds around
% the given trigger number
%
% Script returns eeg data structure: it contains the data, labels etc
% for the given interval
%
% t is the time in milliseconds of the trial (at the trigger, i.e., stimulus-time)
%
% Author: Michiel van Burik, ANT Software, Enschede, The Netherlands, 8 October 2003
%
% See also READ_EEP_TRG, READ_EEP_REJ, READ_EEP_AVR
%

% Copyright (C) Michiel van Burik, ANT Software BV
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

% $Log: read_eep_trial.m,v $
% Revision 1.2  2005/06/08 08:16:37  mvelde
% converted files to unix format
%
% Revision 1.1  2004/11/26 13:17:02  jwiskerke
% Added m-files without binary code in maple distribution.
%
% Revision 1.3  2003/10/24 13:38:50  Maarten-Jan Hoeve
% Added GNU Licence and updated revision history
%
% Revision 1.1.2  2003/10/17 09:55:20  mvelde
% updated: consistent copyrights, arguments/data labels, fixed some typos
%
% ANT Software BV, The Netherlands, www.ant-software.nl / info@ant-software.nl
%

% read short piece of data to get sampling rate, channels etc
eeg = read_eep_cnt([filename '.cnt'],100,101);
samplerate = eeg.rate;
eeglabels = eeg.label;

% read trigger information from external trigger file
trg = read_eep_trg([filename '.trg']);
if triggernumber > length(trg)
	error('Invalid trigger number, trigger does not exist!'); return
else
   trg = trg(triggernumber);
end

% compute interval to extract from file in milliseconds
sample1 = trg.time + interval(1)*1000;
sample2 = trg.time + interval(2)*1000;

% convert interval to samples
sample1 = sample1 /1000 * samplerate + 1;
sample2 = sample2 /1000 * samplerate + 1;

% read data from file
eeg = read_eep_cnt([filename '.cnt'],sample1,sample2);
t = trg.time;

