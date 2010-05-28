% erpregoutfunc() - sub function of erpregout() used to regress 
%                   out the ERP from the data
%
% Usage:
%   totdiff = erpregout(fact, data, erp);
%
% Inputs:
%   fact    - factor 
%   data    - [float] 1-D data (time points).
%   erp     - [float] 1-D data (time points).
%
% Outputs:
%   totdif  - residual difference
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

function totdiff = erpregoutfunc(fact, data, erp);

    totdiff = mean(abs(data - fact*erp));
    