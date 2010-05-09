% dipfit_erpeeg - fit multiple component dipoles using DIPFIT 
%
% Usage:
%         >> [ dipole model EEG] = dipfit_erpeeg(data, chanlocs, 'key', 'val', ...);
%
% Inputs:
%  data      - input data [channel x point]. One dipole per point is
%              returned.
%  chanlocs  - channel location structure (returned by readlocs()).
%
% Optional inputs:
%  'settings'  - [cell array] dipfit settings (arguments to the 
%                pop_dipfit_settings() function). Default is none.
%  'dipoles'   - [1|2] use either 1 dipole or 2 dipoles contrain in
%                symetry. Default is 1.
%  'dipplot'   - ['on'|'off'] plot dipoles. Default is 'off'.
%  'plotopt'   - [cell array] dipplot() 'key', 'val' options. Default is
%                'normlen', 'on', 'image', 'fullmri'
%
% Outputs:
%  dipole      - dipole structure ('posxyz' field is the position; 'momxyz'
%                field is the moment and 'rv' the residual variance)
%  model       - structure containing model information ('vol.r' field is 
%                radius, 'vol.c' conductances, 'vol.o' the 3-D origin and
%                'chansel', the selected channels).
%  EEG         - faked EEG structure containing erp activation at the place
%                of ICA components but allowing to plot ERP dipoles.
%
% Note: residual variance is set to NaN if Dipfit does not converge
%  
% Author: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, Nov. 2003

% Copyright (C) 10/2003 Arnaud Delorme, SCCN/INC/UCSD, arno@sccn.ucsd.edu
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

function [dipoles, model, EEG] = dipfit_erpeeg(DATA, chanlocs, varargin);
    
    if nargin < 1
        help dipfit_erpeeg;
        return;
    end;
    
    ncomps = size(DATA,2);
    if size(DATA,1) ~= length(chanlocs)
        error('# of row in ''DATA'' must equal # of channels in ''chanlocs''');
    end;
        
    % faking an EEG dataset
    % ---------------------
    EEG          = eeg_emptyset;
    EEG.data     = rand(size(DATA,1), 1000);
    EEG.nbchan   = size(DATA,1);
    EEG.pnts     = 1000;
    EEG.trials   = 1;
    EEG.chanlocs = chanlocs;
    EEG.icawinv    = [ DATA DATA ];
    EEG.icaweights = zeros(size([ DATA DATA ]))';
    EEG.icasphere  = zeros(size(DATA,1), size(DATA,1));
    %EEG            = eeg_checkset(EEG);
    EEG.icaact     = EEG.icaweights*EEG.icasphere*EEG.data(:,:);
    EEG.icaact     = reshape( EEG.icaact, size(EEG.icaact,1), size(EEG.data,2), size(EEG.data,3));
    
    % uses mutlifit to fit dipoles
    % ----------------------------
    EEG            = pop_multifit(EEG, [1:ncomps], varargin{:});
    
    % process outputs
    % ---------------
    dipoles = EEG.dipfit.model;
    if isfield(dipoles, 'active')
        dipoles = rmfield(dipoles, 'active');
    end;
    if isfield(dipoles, 'select')
        dipoles = rmfield(dipoles, 'select');
    end;
    model   = EEG.dipfit;
    if isfield(model, 'model')
        model   = rmfield(model, 'model');
    end;
    return;
    
