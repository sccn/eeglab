function [eeg] = read_eep_cnt(fn, varargin);

% READ_EEP_CNT reads continuous EEG data from an EEProbe *.cnt file
% and returns a structure containing the header and data information.
%
% eeg = read_eep_cnt(filename, sample1, sample2)
%
% where sample1 and sample2 are the begin and end sample of the data
% to be read.
%
% eeg.label    ... labels of EEG channels
% eeg.rate     ... sampling rate
% eeg.npnt     ... number of sample in data segment
% eeg.nchan    ... number of channels
% eeg.nsample  
% eeg.time     ... array [1 x npnt] of time points (ms)
% eeg.data     ... array [nchan x npnt] containing eeg data (uV) 
%
% Author: Robert Oostenveld, Aalborg University, Denmark, 11 March 2003
%
% See also READ_EEP_TRG, READ_EEP_REJ, READ_EEP_AVR
%

% Copyright (C) 2002, Robert Oostenveld
%                     Aalborg University, Denmark
%                     http://www.smi.auc.dk/~roberto/
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

error('could not locate mex file');
