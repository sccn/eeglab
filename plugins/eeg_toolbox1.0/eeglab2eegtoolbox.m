% eeglab2eegtoolbox() - convert EEGLAB structure to EEG toolbox
%                       structure (for ERP peak picking mainly).
%
% Usage:
%   >> p = eeglab2eegtoolbox( EEG );
%
% Inputs:
%   EEG - EEG structure
%
% Outputs:
%   p   - EEG toolbox structure
%
% Author: Arnaud Delorme, SCCN, INC, UCSD
%
% See also: eeglab()

% Copyright (C) 2006 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function p = eeglab2eegtoolbox( EEG );
    
    if nargin < 1
        help eeglab2eegtoolbox;
        return;
    end;
    
    % default structure
    % -----------------
    p = eeg_toolbox_defaults;
    
    % exporting importing electrodes
    % ------------------------------
    fid = fopen('tmp.txt', 'w');
    % fiducials should also be added here
    for index = 1:length(EEG.chanlocs)
        tmp = EEG.chanlocs(index);
        fprintf(fid, '%s\t%d\t%f\t%f\t%f\n', tmp.labels, 69, tmp.X, tmp.Y, tmp.Z);
    end;
    fprintf(fid, 'Centroid\t%d\t%f\t%f\t%f\n\n', 99, 0, 0, 0);
    fclose(fid);
    
    p.elec.path = pwd;
    p.elec.file = 'tmp.txt';
    p.elec.n    = length(EEG.chanlocs)+1;
    p = elec_open(p);
    delete(fullfile(pwd, 'tmp.txt'));
    
    % importing data
    % --------------
    p.volt.data = double(mean(EEG.data,3)');
    p.volt.timeArray = EEG.times';
    p.volt.points    = EEG.pnts;
    p.volt.sampleHz  = EEG.srate;
    p.volt.sampleTime = EEG.srate;
    p.volt.sampleMsec = 1000/EEG.srate;
    p.volt.channels   = length(EEG.chanlocs);
    p.volt.epochStart = EEG.xmin*1000;
    p.volt.epochEnd   = EEG.xmax*1000;
    
