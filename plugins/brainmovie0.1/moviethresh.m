% moviethresh() - threshold cell arrays that will be sent to brainmovie
%
% Usage:
%   >> data = moviethresh( data, threshold, continous, dimcont);
%
% Inputs:
%   data       - input cell array containing data arrays
%   threshold  - absolute threshold (0=none)
%   continuous - impair number of data points that must be sequentially
%                organized (3,5,7,9 ...)
%   dimcont    - dimention (1=rows, 2=columns) of the data arrays to 
%                consider for sequential continuity (default: 2) 
%
% Outputs:
%   data       - modified cell array
%
% Example:
%   d = moviethresh( data, 0.1, 3, 2);
%   % remove (i.e) any (absolute) value in the array contained in data
%   % that are below 0.1 and remove any value that is flanked by  
%
% See also: BRAINMOVIE, TIMECROSSF

% arno@salk.edu, Arnaud Delorme, CNL / Salk Institute, 2001

% This program is free software; you can redistribute it and/or
% modify it.  
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

function data = moviethresh( data, threshold, cnt, dimcontinous);

if nargin < 4
	help moviethresh;
end;	

% scan datas
% ----------
nbremoved = 0;
total     = 0;
for index = 1:length(data(:))
	tmp = data{ index };

	% remove values below threshold
	% -----------------------------
	tmp( find(abs(tmp) < threshold) ) = 0;

	% %%%%%%%%%%%%%%%%% TRANSPOSE if dimention 2 
	if dimcontinous == 2, tmp = tmp'; end;

	% remove values that are sparse along dimention 1
	% ----------------------------------------------- 
	s = size(tmp, 1);
	for dimcont = 1:s
		nbcontrigth = zeros(1, size(tmp,2));
		nbcontleft  = zeros(1, size(tmp,2));
		for dim2 = 1:size(tmp,2)
			for ind=1:cnt 									% cont elements on the left
				if dimcont-ind > 0  
					if tmp(dimcont-ind , dim2) ~= 0, nbcontleft(dim2)  = nbcontleft(dim2)  + 1;
					else, break; end;
				else, break; end;
			end;
			for ind=1:cnt  									% count elements on the right
				if dimcont+ind <= s
					if tmp(dimcont+ind , dim2) ~= 0, nbcontrigth(dim2)  = nbcontrigth(dim2)  + 1;
					else, break; end;
				else, break; end;
			end;
		end;	

		total     = total + sum(tmp(dimcont,:) ~= 0); % non nul values
		
		nbremoved = nbremoved + sum( (tmp(dimcont,:) .*  ((nbcontleft+nbcontrigth) < cnt)) > 0);
		tmp(dimcont,:) = tmp(dimcont,:) .* ( (nbcontleft+nbcontrigth) >= cnt); 
		%fprintf( 'Pos %d: %s ** %s\n', dimcont, int2str( nbcontleft), int2str( nbcontrigth));
	end;	

	% %%%%%%%%%%%%%%%%% TRANSPOSE if dimention 2 
	if dimcontinous == 2, tmp = tmp'; end;

	data{ index } = tmp;

end;	
fprintf('%d removed values out of %d\n', nbremoved, total)

%debug fprintf('\tval:%f\tdestval:%f\tdimcont:%d dim2:%d ind:%d left:%d\n', tmp(dimcont, dim2), tmp(dimcont-ind , dim2), dimcont, dim2, ind, nbcontleft(dim2));
