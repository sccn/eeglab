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

function [data, erp, factors] = erpregout(data, tlim, reglim);

    if nargin < 1
        help erpregout;
        return;
    end
    if nargin < 2
        tlim = [0 1];
    end
    if nargin < 3
        reglim = tlim;
    end
    if ndims(data) == 2
        data = reshape(data, 1, size(data,1), size(data,2));
        redim = 1;
    else
        redim = 0;
    end
    
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
            if ~mod(trial, 10) , fprintf('%d ', trial); end
            if ~mod(trial, 200), fprintf('\n', trial); end
            [factors(chan, trial) tmpf exitflag] = fminbnd('erpregoutfunc', 0, 10, [], ...
                                   data(chan, begpoint:endpoint, trial), erp(chan, begpoint:endpoint));
            data(chan,:,trial) = data(chan,:,trial) - factors(chan, trial)*erp(chan, :);
        end
        fprintf('\n');
    end
    
    if redim
        data = squeeze(data);
    end
