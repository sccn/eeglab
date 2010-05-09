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

function [ probaMap, sortbox ] = realproba( data, bins );

if nargin < 1
	help realproba;
	return;
end;
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
    if any(any(isnan(data))), warning('Binning failed - could be due to zeroed out channel'); end;
	for index=1:SIZE
		sortbox(data(index)) = sortbox(data(index))+1;
	end;
	probaMap = sortbox(data) / SIZE;
	sortbox  = sortbox / SIZE;
else
	% BASE OVER ERROR FUNCTION
	% ------------------------
	data     = (data-mean(data(:)))./std(data(:));
	probaMap = exp(-0.5*( data.*data ))/(2*pi);
	probaMap = probaMap/sum(probaMap); % because the total surface under a normalized Gaussian is 2
	sortbox  = probaMap/sum(probaMap);
end;
return;
