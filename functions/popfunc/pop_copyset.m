% pop_copyset() - Copy the current EEG dataset into another dataset.
%
% Usage:
%   >> ALLEEG = pop_copyset(ALLEEG, index1); % pop-up
%   >> [ ALLEEG EEG CURRENTSET ] = pop_copyset(ALLEEG, index1, index2 );
%
% Inputs:
%   ALLEEG     - array of dataset structure
%   index1     - input dataset number
%   index2     - index of dataset to copy into
%
% Inputs:
%   ALLEEG     - array of dataset structures
%   EEG        - new copied structure
%   CURRENTSET - index of the new dataset
%
% Note: this function performs ALLEEG(index2) = ALLEEG(index1);
%       with dataset checks
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeg_store(), pop_delset(), eeglab() 

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% 01-25-02 reformated help & license -ad 

function [ALLEEG, EEG, CURRENTSET, com] = pop_copyset(ALLEEG, set_in, set_out);

com = '';
if nargin < 2
	help pop_copyset;
	return;
end
if isempty(ALLEEG)
	error(['Pop_copyset error: cannot copy' 10 'single dataset mode']);
end;	
if set_in == 0
    error('Pop_copyset error: cannot copy dataset'); return;
end
if isempty(ALLEEG(set_in).data)
    error('Pop_copyset error: cannot copy empty dataset'); return;
end

if nargin < 3
	% which set to save
	% -----------------
	promptstr    = { 'Index of the new dataset:'};
	inistr       = { int2str(set_in+1) };
	result       = inputdlg2( promptstr, 'Copy dataset -- pop_copyset()', 1,  inistr, 'pop_copyset');
	if size( result ) == 0, EEG = []; CURRENTSET = 0; return; end
	set_out   	 = eval( result{1} );
end
ALLEEG = eeg_store(ALLEEG, eeg_retrieve(ALLEEG, set_in), set_out);
EEG    = eeg_retrieve(ALLEEG, set_out);
CURRENTSET = set_out;

com = sprintf('[ALLEEG EEG CURRENTSET] = pop_copyset( ALLEEG, %d, %d);', set_in, set_out);
return;
