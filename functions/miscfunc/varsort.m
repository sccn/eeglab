% varsort() - reorder ICA components, largest to smallest, by 
%             the size of their MEAN projected variance 
%             across all time points
% Usage:
%   >> [windex,meanvar] = varsort(activations,weights,sphere);
%
% Inputs:
%   activations = (chans,framestot) the runica() activations
%   weights     = ica weight matrix from runica() 
%   sphere      = sphering matrix from runica() 
%
% Outputs:
%   windex   = order of projected component mean variances (large to small)
%   meanvar  = projected component mean variance (in windex order)
%
% Author: Scott Makeig & Martin McKeown, SCCN/INC/UCSD, La Jolla, 09-01-1997 
%
% See also: runica()

% Copyright (C) 9-01-1997 Scott Makeig & Martin McKeown, SCCN/INC/UCSD, 
% scott@sccn.ucsd.edu
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

% 03-19-97 simplified, replaced grandmean with datamean info in calculation, 
%          made function return mean projected variance across the data, 
%          changed var() to diag(cov()) -sm
% 05-20-97 use sum-of-squares instead of diag() to allow long data sets -sm
% 06-07-97 changed order of args to conform to runica, fixed meanvar computation -sm
% 07-25-97 removed datamean -sm
% 01-25-02 reformated help & license, added link -ad 

function [windex,meanvar] = varsort(activations,weights,sphere)
%
if nargin ~= 3     % needs all 3 args
     help varsort
     return
end
[chans,framestot] = size(activations);
if framestot==0,
    fprintf('Gvarsort(): cannot process an empty activations array.\n\n');
    return
end

[srows,scols] = size(sphere);
[wrows,wcols] = size(weights);

if nargin<3,
    fprintf('Gvarsort(): needs at least 3 arguments.\n\n');
    return
end

% activations = (wrows,wcols)X(srows,scols)X(chans,framestot)
if chans ~= scols || srows ~= wcols,
   fprintf('varsort(): input data dimensions do not match.\n');
   fprintf('              i.e., Either %d ~= %d or %d ~= %d\n',...
                                     chans,scols,srows,wcols);
   return
end

%%%%%%%%%%%%%%%%%%%% Find mean variances %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Computing mean projected variance for all %d components:\n',wrows);
meanvar  = zeros(wrows,1);  % size of the projections
winv = inv(weights*sphere);
for s=1:wrows
     fprintf('%d ',s);      % construct single-component data matrix
                            % project to scalp, then add row means 
    compproj = winv(:,s)*activations(s,:);
    meanvar(s) = mean(sum(compproj.*compproj)/(size(compproj,1)-1));
                            % compute mean variance 
end                         % at all scalp channels

%%%%%%%%%%%%%%%%%%% sort by mean variance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sortvar, windex] = sort(meanvar);
windex = windex(wrows:-1:1);% order large to small 
meanvar = meanvar(windex);
fprintf('\n');
