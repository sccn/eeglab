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
end

if nargin < 2
	help quantile;
	return;
end;	

[prows pcols] = size(pc);
if prows ~= 1 && pcols ~= 1
    error('pc must be a scalar or a vector.');
end
if any(pc > 1) || any(pc < 0)
    error('pc must be between 0 and 1');
end
[i,j] = size(data);
sortdata = sort(data);
if i==1 || j==1 % if data is a vector
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
