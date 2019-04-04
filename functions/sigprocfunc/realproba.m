% realproba() - compute the effective probability of the value 
%               in the sample.
%
% Usage: 
%   >> [probaMap, probaDist ] = realproba( data, discret);
%
% Inputs:
%   data       - the data onto which compute the probability
%   discret    - discretisation factor (default: (size of data)/5)
%                if 0 base the computation on a Gaussian 
%                approximation of the data 
%
% Outputs:
%   probaMap   - the probabilities associated with the values
%   probaDist  - the probabilities distribution 
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

function [ probaMap, sortbox ] = realproba( data, bins );

if nargin < 1
	help realproba;
	return;
end
if nargin < 2
	bins = round(size(data,1)*size(data,2)/5);
end;	

if bins > 0
	% COMPUTE THE DENSITY FUNCTION
	% ----------------------------
	SIZE = size(data,1)*size(data,2);
	sortbox = zeros(1,bins);
	minimum =  min(data(:));
	maximum =  max(data(:));
	data = floor((data - minimum )/(maximum - minimum)*(bins-1))+1;
    if any(any(isnan(data))), warning('Binning failed - could be due to zeroed out channel'); end
	for index=1:SIZE
		sortbox(data(index)) = sortbox(data(index))+1;
	end
	probaMap = sortbox(data) / SIZE;
	sortbox  = sortbox / SIZE;
else
	% BASE OVER ERROR FUNCTION
	% ------------------------
	data     = (data-mean(data(:)))./std(data(:));
	probaMap = exp(-0.5*( data.*data ))/(2*pi);
	probaMap = probaMap/sum(probaMap); % because the total surface under a normalized Gaussian is 2
	sortbox  = probaMap/sum(probaMap);
end
return;
