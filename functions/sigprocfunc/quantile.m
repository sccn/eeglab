% quantile() - computes the quantiles of the data sample from a distribution X
%
% Description:
%      If F is the cumulative distribution function (CDF) of X,
%      the p-th-quantile Qp of distribution X is the value for which holds
%                        F(x) <  p,  for x < Qp, and 
%                        F(x) >= p,  for x >= Qp.
%      for example, for p = 0.5, Qp is the median of X. p must be in [0..1].
%
% Usage:
%   >>  q = quantile( data, pc );
%
% Inputs:
%   data   - vector of observations
%   pc     - the quantiles will be estimated at the values in pc [0..1]
%    
% Outputs:
%   q      - pc-th-quantiles of the distribution generating the observation
%
% Authors: Scott Makeig & Luca Finelli, CNL/Salk Institute-SCCN, August 21, 2002
%
% Note: this function overload the function from the statistics toolbox. In
%       case the statistic toolbox is present, the function from the
%       statistics toolbox is being called instead.
%
% See also: 
%   pop_sample(), eeglab() 

% Copyright (C) Scott Makeig & Luca Finelli, CNL/Salk Institute-SCCN, August 21, 2002
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

function q = quantile(data,varargin); 

% detect overloaded method in stat toolbox
curPath = fileparts(which('quantile'));
rmpath(curPath);
path2   = fileparts(which('quantile'));
addpath(curPath);
if ~isempty(path2)
    addpath(path2);
    q = quantile(data,varargin{:});
    return;
else
    pc = varargin{1};
end;

if nargin < 2
	help quantile;
	return;
end;	

[prows pcols] = size(pc);
if prows ~= 1 & pcols ~= 1
    error('pc must be a scalar or a vector.');
end
if any(pc > 1) | any(pc < 0)
    error('pc must be between 0 and 1');
end
[i,j] = size(data);
sortdata = sort(data);
if i==1 | j==1 % if data is a vector
    i = max(i,j); j = 1;
    if i == 1,
        fprintf('  quantile() note: input data is a single scalar!\n')
        q = data*ones(length(pc),1); % if data is scalar, return it
        return;
    end
    sortdata = sortdata(:);
end
pt = [0 ((1:i)-0.5)./i 1];
sortdata = [min(data); sortdata; max(data)];
q = interp1(pt,sortdata,pc);
