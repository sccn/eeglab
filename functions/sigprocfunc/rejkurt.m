% rejkurt()  - calculation of kutosis of a 1D, 2D or 3D array and
%              rejection of outliers values of the input data array   
%              using the discrete kutosis of the values in that dimension.
%
% Usage:
%   >>  [kurtosis rej] = rejkurt( signal, threshold, kurtosis, normalize);
%
% Inputs:
%   signal     - one dimensional column vector of data values, two 
%                dimensional column vector of values of size 
%                sweeps x frames or three dimensional array of size 
%                component x sweeps x frames. If three dimensional, 
%                all components are treated independently. 
%   threshold  - Absolute threshold. If normalization is used then the 
%                threshold is expressed in standard deviation of the
%                mean. 0 means no threshold.
%   kurtosis   - pre-computed kurtosis (only perform thresholding). Default
%                is the empty array [].
%   normalize  - 0 = do not not normalize kurtosis. 1 = normalize kurtosis.
%                2 is 20% trimming (10% low and 10% high) kurtosis before 
%                normalizing. Default is 0.
% 
% Outputs:
%   kurtosis    - normalized joint probability  of the single trials 
%                (same size as signal without the last dimension)
%   rej         - rejected matrix (0 and 1, size: 1 x sweeps)
%
% Remarks:
%   The exact values of kurtosis depend on the size of a time 
%   step and thus cannot be considered as absolute.
%   This function uses the kurtosis function from the statistival
%   matlab toolbox. If the statistical toolbox is not installed, 
%   it uses the 'kurt' function of the ICA/EEG toolbox.
%
% See also: kurt(), kurtosis()

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

function [kurto, rej] = rejkurt( signal, threshold, oldkurtosis, normalize);

if nargin < 1
	help rejkurt;
	return;
end;	
if nargin < 2
	threshold = 0;
end;	
if nargin < 4
	normalize = 0;
end;	
if nargin < 3
	oldkurtosis = [];
end;	

if size(signal,2) == 1 % transpose if necessary
	signal = signal';
end;

nbchan = size(signal,1);
pnts = size(signal,2);
sweeps = size(signal,3);
kurto = zeros(nbchan,sweeps);

if ~isempty( oldkurtosis ) % speed up the computation
	kurto = oldkurtosis;
else
	for rc = 1:nbchan
		% compute all kurtosis
		% --------------------
		for index=1:sweeps
			try 
			    kurto(rc, index) = kurtosis(signal(rc,:,index));
			catch
				kurto(rc, index) = kurt(signal(rc,:,index));
			end;	
		end;
	end;

	% normalize the last dimension
	% ----------------------------	
	if normalize
        tmpkurt = kurto;
        if normalize == 2,
            tmpkurt = sort(tmpkurt);
            minind  = max(round(length(tmpkurt)*0.1),1);
            maxind  = round(length(tmpkurt)-round(length(tmpkurt)*0.1));
            if size(tmpkurt,2) == 1
                 tmpkurt = tmpkurt(minind:maxind);
            else tmpkurt = tmpkurt(:,minind:maxind);
            end;
        end;
	    switch ndims( signal )
	    	case 2,	kurto = (kurto-mean(tmpkurt)) / std(tmpkurt);
	    	case 3,	kurto = (kurto-mean(tmpkurt,2)*ones(1,size(kurto,2)))./ ...
				        (std(tmpkurt,0,2)*ones(1,size(kurto,2)));
		end;
	end;
end;

% reject
% ------	
if threshold(1) ~= 0 
    if length(threshold) > 1
    	rej = (threshold(1) > kurto) | (kurto > threshold(2));
    else
    	rej = abs(kurto) > threshold;
    end;
else
	rej = zeros(size(kurto));
end;	

return;
