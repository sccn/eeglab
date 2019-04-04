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
