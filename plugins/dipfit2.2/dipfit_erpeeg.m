% dipfit_erpeeg - fit multiple component dipoles using DIPFIT 
%
% Usage:
%         >> [ dipole model ] = dipfit_erpeeg(data, chanlocs, 'key', 'val', ...);
%
% Inputs:
%  data      - input data [channel x point]. One dipole per point is
%              returned.
%  chanlocs  - channel location structure (returned by readlocs()).
%
% Optional inputs:
%  'settings'  - [cell array] dipfit settings. Default is none.
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
%
% Note: residual variance is set to NaN if Dipfit does not converge
%  
% Author: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, Nov. 2003

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.1  2003/11/06 18:28:54  arno
% Initial revision
%

function [dipoles, model] = dipfit_erpeeg(DATA, chanlocs, varargin);
    
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
    EEG.data     = zeros(size(DATA,1), 1000);
    EEG.nbchan   = size(DATA,1);
    EEG.pnts     = 1000;
    EEG.trials   = 1;
    EEG.chanlocs = chanlocs;
    EEG.icawinv    = DATA;
    EEG.icaweights = zeros(size(DATA))';
    EEG.icasphere  = zeros(size(DATA,1), size(DATA,1));
    
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
    