% gui_erp_plot_eeglab - convert EEGLAB structure to EEG toolbox
%                       structure and plot ERP for ERP peak picking.
% Usage:
%   >> p = gui_erp_plot_eeglab( EEG );
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

function p = gui_erp_plot_eeglab( EEG );
    
    if nargin < 1
        help eeglab2eegtoolbox;
        return;
    end;
        
    p = eeglab2eegtoolbox( EEG );
    
    % print number-label correspondance
    % ---------------------------------
    disp('Channel indices are indicated in the top right corner of the GUI')
    disp('The correspondance between channel indices and labels is:')
    TMP = EEG.chanlocs;
    
    for index = 1:length(EEG.chanlocs)
        fprintf('%d: %s\n', index, EEG.chanlocs(index).labels);
    end;
    
    p = gui_erp_plot(p);
