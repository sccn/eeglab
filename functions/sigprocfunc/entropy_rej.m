% entropy_rej() - calculation of entropy of a 1D, 2D or 3D array and
%             rejection of odd last dimension values of the input data array   
%             using the discrete entropy of the values in that dimension
%             (and using the probability distribution of all columns).
%
% Usage:
%   >>  [entropy rej] = entropy_rej( signal, threshold, entropy, normalize, discret);
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
%   entropy    - pre-computed entropy_rej (only perform thresholding). Default
%                is the empty array [].
%   normalize  - 0 = do not not normalize entropy. 1 = normalize entropy.
%                Default is 0.
%   discret    - discretization variable for calculation of the 
%                discrete probability density. Default is 1000 points. 
% 
% Outputs:
%   entropy    - entropy (normalized or not) of the single data trials 
%                (same size as signal without the last dimension)
%   rej        - rejection matrix (0 and 1, size of number of rows)
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: realproba()

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

function [ent, rej] = entropy_rej( signal, threshold, oldentropy_rej, normalize, discret );

if nargin < 1
	help entropy_rej;
	return;
end;	
if nargin < 2
	threshold = 0;
end;	
if nargin < 3
	oldentropy_rej = [];
end;	
if nargin < 4
	normalize = 0;
end;	
if nargin < 5
	discret = 1000;
end;	
%	threshold = erfinv(threshold);

if size(signal,2) == 1 % transpose if necessary
	signal = signal';
end;

[nbchan pnts sweeps] = size(signal);
ent  = zeros(nbchan,sweeps);

if ~isempty( oldentropy_rej ) % speed up the computation
	ent = oldentropy_rej;
else
	for rc = 1:nbchan

		% COMPUTE THE DENSITY FUNCTION
		% ----------------------------
		[ dataProba sortbox ] = realproba( signal(rc, :), discret );

		% compute all entropy
		% -------------------
		for index=1:sweeps
			datatmp = dataProba((index-1)*pnts+1:index*pnts);
			ent(rc, index) = - sum( datatmp .* log( datatmp ) ); 
		end;
	end;

	% normalize the last dimension
	% ----------------------------	
	if normalize
	    switch ndims( signal )
	    	case 2,	ent = (ent-mean(ent)) / std(ent);
	    	case 3,	ent = (ent-mean(ent,2)*ones(1,size(ent,2)))./ ...
				        (std(ent,0,2)*ones(1,size(ent,2)));
		end;
	end;
end	

% reject
% ------	
if threshold ~= 0 
	rej = abs(ent) > threshold;
else
	rej = zeros(size(ent));
end;	

return;
