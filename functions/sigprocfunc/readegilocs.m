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

function EEG = readegilocs(EEG, fileloc);

if nargin < 1
    help readegilocs;
    return;
end;

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
    end;
end;
if found
    fprintf('EGI channel location automatically detected %s ********* WARNING please check that this the proper file\n', fileloc);
    if EEG.nbchan == 64 || EEG.nbchan == 65 || EEG.nbchan == 256 || EEG.nbchan == 257
        fprintf( [ 'Warning: this function assumes you have a 64-channel system Version 2\n' ...
                   '         if this is not the case, update the channel location with the proper file' ]);
    end;
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
        end;
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
end;
