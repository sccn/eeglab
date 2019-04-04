% icaact() - compute ICA activation waveforms = weights*sphere*(data-meandata)
%
% Usage: >> [activations] = icaact(data,weights,datamean);
%
% Inputs:  
%     data     = input data (chans,frames)
%     weights  = unmixing matrix (runica() weights*sphere)
%     datamean = 0 or mean(data')  (default 0);
%
% Note:  If datamean==0, data means are distributed over activations.
%        Use this form for plotting component projections.
%
% Output:  
%        activations = ICA component activation waveforms 
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 4-3-97 
%
% See also: runica(), icaproj(), icavar()

% Copyright (C) 4-3-97 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 6-17-97 extended to non-square weight matrices -sm
% 1-12-01 removed sphere argument -sm
% 01-25-02 reformated help & license, added links -ad 

function [activations] = icaact(data,weights,datamean)

if nargin < 4
    datamean = 0;
elseif nargin < 3
    help icaact
    return
end

[chans, framestot] = size(data);

if datamean == 0,
    datamean = zeros(chans,1); % single-epoch 0s
end

if size(datamean,1) == 1    % if row vector
    datamean = datamean';   % make a column vector
end
[meanchans,epochs] = size(datamean);
if epochs < 1,
	fprintf('icaact(): datamean empty.\n');
	return
end
frames = fix(framestot/epochs);

if frames < 1,
	fprintf('icaact(): data empty.\n');
	return
end

if frames*epochs ~= framestot
	fprintf(...
   'icaact(): datamean epochs %d does not divide data length %d.\n',...
                          epochs,                           framestot);
	return
end

if size(datamean,1) ~= chans
	fprintf('icaact(): datamean channels ~= data channels.\n');
	return
end

w = weights;
activations = zeros(size(w,1),size(data,2));
for e=1:epochs
	activations(:,(e-1)*frames+1:e*frames) =  ...
        w*(data(:,(e-1)*frames+1:e*frames) - datamean(:,e)*ones(1,frames));
end
