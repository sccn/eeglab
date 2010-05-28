% erpregout() - regress out the ERP from the data
%
% Usage:
%   newdata = erpregout(data);
%   [newdata erp factors] = erpregout(data, tlim, reglim);
%
% Inputs:
%   data    - [float] 2-D data (times x trials) or 3-D data
%             (channels x times x trials).
%
% Optional inputs:
%   tlim    - [min max] time limits in ms.
%   reglim  - [min max] regression time window in ms (by default
%             the whole time period is used
% Outputs:
%   newdata - data with ERP regressed out
%   erp     - data ERP
%   factors - factors used for regressing out the ERP (size is the same
%             as the number of trials or (channels x trials)
%
% Note: it is better to regress out the ERP about 4 times (launch the
%       function 4 times in a row) to really be able to regress out the
%       ERP and have a residual ERP close to 0.
%
% Author: Arnaud Delorme, Salk, SCCN, UCSD, CA, April 29, 2004

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2004 Arnaud Delorme
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

function [data, erp, factors] = erpregout(data, tlim, reglim);

    if nargin < 1
        help erpregout;
        return;
    end;
    if nargin < 2
        tlim = [0 1];
    end;
    if nargin < 3
        reglim = tlim;
    end;
    if ndims(data) == 2
        data = reshape(data, 1, size(data,1), size(data,2));
        redim = 1;
    else
        redim = 0;
    end;
    
    % find closest points
    % -------------------
    timevect = linspace(tlim(1), tlim(2), size(data,2));
    [tmp begpoint] = min( abs(timevect-reglim(1)) );    
    [tmp endpoint] = min( abs(timevect-reglim(2)) );
    erp = mean(data, 3);
    
    % regressing out erp in channels and trials
    % -----------------------------------------
    for chan = 1:size(data,1)
        fprintf('Channel %d (trials out of %d):', chan, size(data,3));
        for trial = 1:size(data,3)
            if ~mod(trial, 10) , fprintf('%d ', trial); end;
            if ~mod(trial, 200), fprintf('\n', trial); end;
            [factors(chan, trial) tmpf exitflag] = fminbnd('erpregoutfunc', 0, 10, [], ...
                                   data(chan, begpoint:endpoint, trial), erp(chan, begpoint:endpoint));
            data(chan,:,trial) = data(chan,:,trial) - factors(chan, trial)*erp(chan, :);
        end;
        fprintf('\n');
    end;
    
    if redim
        data = squeeze(data);
    end;