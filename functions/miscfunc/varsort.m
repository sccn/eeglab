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
end;

[srows,scols] = size(sphere);
[wrows,wcols] = size(weights);

if nargin<3,
    fprintf('Gvarsort(): needs at least 3 arguments.\n\n');
    return
end;

% activations = (wrows,wcols)X(srows,scols)X(chans,framestot)
if chans ~= scols | srows ~= wcols,
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
