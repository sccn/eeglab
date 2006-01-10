% eeglab2fieldtrip() - do this ...
%
% Usage:    >> data = eeglab2fieldtrip( EEG, fieldbox, transform );
%
% Inputs:
%   EEG       - [struct] EEGLAB structure
%   fieldbox  - ['preprocessing'|'freqanalysis'|'timelockanalysis'|'companalysis']
%   transform - ['none'|'dipfit'] transform channel locations for DIPFIT
%               using the transformation matrix in the field
%               'coord_transform' of the dipfit substructure of the EEG
%               structure.
% Outputs:
%   data    - FIELDTRIP structure
%
% Author: Robert Oostenveld, F.C. Donders Centre, May, 2004.
%         Arnaud Delorme, SCCN, INC, UCSD
%
% See also: 

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2004 Robert Oostenveld, F.C. Donders Centre, roberto@smi.auc.dk
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

function data = eeglab2fieldtrip(EEG, fieldbox, transform)

if nargin < 2
  help eeglab2fieldtrip
  return;
end;

% start with an empty data object 
data = [];

% add the objects that are common to all fieldboxes
data.label = { EEG.chanlocs(1:EEG.nbchan).labels };
data.fsample = EEG.srate;

% get the electrode positions from the EEG structure: in principle, the number of 
% channels can be more or less than the number of channel locations, i.e. not 
% every channel has a position, or the potential was not measured on every
% position. This is not supported by EEGLAB, but it is supported by FIELDTRIP.

data.elec.pnt   = zeros(length( EEG.chanlocs ), 3);
for ind = 1:length( EEG.chanlocs )
    data.elec.label{ind} = EEG.chanlocs(ind).labels;
    if ~isempty(EEG.chanlocs(ind).X)
        data.elec.pnt(ind,1) = EEG.chanlocs(ind).X;
        data.elec.pnt(ind,2) = EEG.chanlocs(ind).Y;
        data.elec.pnt(ind,3) = EEG.chanlocs(ind).Z;
    else
        data.elec.pnt(ind,:) = [0 0 0];
    end;
end;

if nargin > 2
    if strcmpi(transform, 'dipfit') 
        if ~isempty(EEG.dipfit.coord_transform)
            disp('Transforming electrode coordinates to match head model');
            transfmat = traditional(EEG.dipfit.coord_transform);
            data.elec.pnt = transfmat * [ data.elec.pnt ones(size(data.elec.pnt,1),1) ]';
            data.elec.pnt = data.elec.pnt(1:3,:)';
        else
            disp('Warning: no transformation of electrode coordinates to match head model');
        end;
    end;
end;
        
switch fieldbox
  case 'preprocessing'
    for index = 1:EEG.trials
      data.trial{index}  = EEG.data(:,:,index);
      data.offset(index) = EEG.xmin*EEG.srate+1;                   % should be checked in FIELDTRIP
      data.time{index}   = linspace(EEG.xmin, EEG.xmax, EEG.pnts); % should be checked in FIELDTRIP
    end;
    
  case 'timelockanalysis'
    data.avg  = mean(EEG.data, 3);   
    data.var  = std(EEG.data, [], 3).^2;   
    data.time = linspace(EEG.xmin, EEG.xmax, EEG.pnts); % should be checked in FIELDTRIP
    
  case 'componentanalysis'
    for index = 1:EEG.trials
      % the trials correspond to the raw data trials, except that they
      % contain the component activations
      data.trial{index}  = EEG.icaact(:,:,index);
      data.offset(index) = EEG.xmin*EEG.srate+1;                   % should be checked in FIELDTRIP
      data.time{index}   = linspace(EEG.xmin, EEG.xmax, EEG.pnts); % should be checked in FIELDTRIP
    end;
    for comp = 1:size(EEG.icawinv,2)
      % the labels correspond to the component activations that are stored in data.trial
      data.label{comp} = sprintf('ica_%03d', comp);
    end
    % get the spatial distribution and electrode positions
    data.topolabel = { EEG.chanlocs(1:EEG.nbchan).labels };
    data.topo      = EEG.icawinv;
    
  case 'freqanalysis'
    error('freqanalysis fieldbox not implemented yet')
    
  otherwise
    error('unsupported fieldbox') 
end

try
  % get the full name of the function
  data.cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  data.cfg.version.name = st(i);
end

% add the version details of this function call to the configuration
data.cfg.version.id   = '$Id: eeglab2fieldtrip.m,v 1.3 2006-01-10 00:45:01 arno Exp $';

return
