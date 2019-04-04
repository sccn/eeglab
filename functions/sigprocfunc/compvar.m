% compvar()   - project selected components and compute the variance of
%               the original data they account for.
%
% Usage:
%   >> [proj, variance] = compvar( data, wts_or_act, winv, components);
%
% Required Inputs:
% data        - 2-D (channels, points) or 3-D (channels, frames, trials)
%               data array.
% wts_or_act  - {sphere weights} cell array containing the ICA sphere 
%               and weights matrices. May also be a 2-D (channels, points) 
%               or 3-D (channels, frames, trials) array of component 
%               activations.
% winv        - inverse or pseudo-inverse of the product of the weights
%               and sphere matrices returned by the ICA decompnumsition,
%               i.e., inv(weights*sphere) or pinv(weights*sphere).
% components  - array of component indices to back-project
%
% Outputs:
%  proj       - summed back-projections of the specified components
%  pvaf       - percent variance of the data that the selected 
%               components account for (range: 100% to -Inf%).
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [ compproj, varegg ] = compvar( data, act, winv, compnums)

if nargin < 4
   help compvar;
   return;
end;   

if iscell(act) && length(act) == 1
    act = act{1};
end

data = reshape(data, size(data,1), size(data,2)*size(data,3));
act  = reshape(act , size(act ,1), size(act ,2)*size( act,3));
squaredata  = sum(sum(data.^2));             % compute the grand sum-squared data

if iscell(act)
    sphere = act{1};
    weight = act{2};
    act = (weight*sphere)*data;
end

compproj   = winv(:,compnums)*act(compnums,:)-data; % difference between data and back-projection
squarecomp = sum(sum(compproj.^2));                 % the summed-square difference
varegg     = 100*(1- squarecomp/squaredata);        % compute pvaf of components in data
compproj   = compproj+data;                         % restore back-projected data

return;

