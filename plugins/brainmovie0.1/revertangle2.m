%                 angles 0->180 and 180->0) to maximize synchronization
%                 between components. 
%
% Usage:
%   >> [crossfangle, inversion] = revertangle( crossfangle, crossfampl );
%
% Inputs:
%   crossfangle        - crossf angle for all pair of components and all
%                        conditions (cell array nbcomp x nbconp x 
%                        conditions)
%   crossfampl         - crossf amplitidue (same structure as above). 
%                        Used to determine non-nul coherences
% Outputs:
%   crossfangle        - idem input, some element in the cell array
%                        having their angle reversed 
%   inversion          - array of 0 and 1 indicating components that
%                        had to be reversed
%
% See also: TIMECROSSF, BRAINMOVIE

% arno@salk.edu, Arnaud Delorme, CNL / Salk Institute, 2001

% This program is free software; you can redistribute it and/or
% modify it.  
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

% revertangle( ALLCROSSFANGLE, ALLCROSSF );
% ALLCROSSFANGLE
% ALLCROSSF

function [ALLCROSSFANGLE, inversion] = revertangle( ALLCROSSFANGLE, ALLCROSSF );

nbcompo = size(ALLCROSSFANGLE,2);
conditions = size(ALLCROSSFANGLE,3);

% invert coherences if necessary
% ------------------------------
combinations = dec2bin( 0:(2^nbcompo-1) );
sumangle = zeros( size(combinations,1), 1);
for comb = 1:size(combinations,1)
	for cnd = 1:conditions
		for index1 = 1:nbcompo
			for index2 = 1:nbcompo
				if index2 > index1

					% significant angles	
					% ------------------
					tmp = ALLCROSSF{ index1, index2, cnd };
					I = find(tmp(:) > 0);
					tmpangle = ALLCROSSFANGLE{ index1, index2, cnd };
					tmpangle = tmpangle(I);
					
					if combinations( comb, index1) == '1'
						if combinations( comb, index2) == '1'
							sumangle(comb) = sumangle(comb) + sum(sum(abs(tmpangle))); % 2 inverted
						else
							sumangle(comb) = sumangle(comb) + sum(sum(abs(tmpangle+180 - (tmpangle > 0).*360))); 
						end;
					else
						if combinations( comb, index2) == '1'
							sumangle(comb) = sumangle(comb) + sum(sum(abs(tmpangle+180 - (tmpangle > 0).*360))); 
							%sumangle(comb) = sumangle(comb) + sum(sum(abs(mod(180+tmpangle,180)))); 
						else
							sumangle(comb) = sumangle(comb) + sum(sum(abs(tmpangle))); % 2 inverted
						end;
					end;						
				end;
			end;
		end;
	end;	
end;

% find the max and revert the specified angles
% --------------------------------------------
for index = 1:length(sumangle)
	fprintf('%s %f\n', combinations(index, :), sumangle(index, :));
end;	
[tmp minI] = min(sumangle);
inversion = combinations( minI, :);
fprintf('Inversion: %s\n', inversion);
for cnd = 1:conditions
	for index1 = 1:nbcompo
		for index2 = 1:nbcompo
			if index2 > index1

				% significant angles	
				% ------------------
				tmp = ALLCROSSF{ index1, index2, cnd };
				I = find(tmp(:) > 0);
				tmpangle = ALLCROSSFANGLE{ index1, index2, cnd };

				if combinations( minI, index1) == '1'
					tmpangle(I) = tmpangle(I)+180 - (tmpangle(I) > 0).*360;
				end;	
				if combinations( minI, index2) == '1'
					tmpangle(I) = tmpangle(I)+180 - (tmpangle(I) > 0).*360;
					% orignally 2steps 
					% tmpangle(I) = 180+tmpangle(I);
					% tmpangle(I) = tmpangle(I) - (tmpangle(I) > 180).*360;
				end;	
				ALLCROSSFANGLE{ index1, index2, cnd } = tmpangle;
			end;
		end;
	end;	
end;

return;		
