% readegilocs() - look up locations for EGI EEG dataset.
%
% Usage:
%   >> EEG = readegilocs(EEG);
%   >> EEG = readegilocs(EEG, fileloc);
%
% Inputs:
%   EEG        - EEGLAB data structure
%
% Optional input:
%   fileloc    - EGI channel location file
%
% Outputs:
%   EEG        - EEGLAB data structure with channel location structure
%
% Author: Arnaud Delorme, SCCN/UCSD
%
% See also: pop_readegi(), readegi()

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

function EEG = readegilocs(EEG, fileloc);

if nargin < 1
    help readegilocs;
    return;
end

% importing channel locations
% ---------------------------
found = 1;
if nargin < 2
    switch EEG.nbchan
        case { 32 33 }, fileloc = 'GSN-HydroCel-32.sfp';
        case { 64 65 }, fileloc = 'GSN65v2_0.sfp';
        case { 128 129 }, fileloc = 'GSN129.sfp';
        case { 256 257 }, fileloc = 'GSN-HydroCel-257.sfp';
        otherwise, found = 0;
    end
end
if found
    fprintf('EGI channel location automatically detected %s ********* WARNING please check that this the proper file\n', fileloc);
    if EEG.nbchan == 64 || EEG.nbchan == 65 || EEG.nbchan == 256 || EEG.nbchan == 257
        fprintf( [ 'Warning: this function assumes you have a 64-channel system Version 2\n' ...
                   '         if this is not the case, update the channel location with the proper file' ]);
    end
    % remove the last channel for 33 channels

    peeglab = fileparts(which('eeglab.m'));
    fileloc = fullfile(peeglab, 'sample_locs', fileloc);
    locs = readlocs(fileloc);
    locs(1).type = 'FID';
    locs(2).type = 'FID';
    locs(3).type = 'FID';
    locs(end).type = 'REF';
    if EEG.nbchan == 256 || EEG.nbchan == 257 
        if EEG.nbchan == 256
            chaninfo.nodatchans = locs([end]);
            locs([end]) = [];
        end
    elseif mod(EEG.nbchan,2) == 0, 
        chaninfo.nodatchans = locs([1 2 3 end]);
        locs([1 2 3 end]) = [];
    else
        chaninfo.nodatchans = locs([1 2 3]);
        locs([1 2 3]) = [];
    end; % remove reference
    chaninfo.filename = fileloc;
    EEG.chanlocs   = locs;
    EEG.urchanlocs = locs;
    EEG.chaninfo   = chaninfo;
end
