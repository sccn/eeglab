% trial2eegplot() - convert eeglab format to eeplot format of rejection window
%
% Usage:
%   >> eegplotarray = trial2eegplot(rej, rejE, points, color);
%
% Inputs:
%   rej   - rejection vector (0 and 1) with one value per trial
%   rejE  - electrode rejection array (size nb_elec x trials) also
%           made of 0 and 1.
%   points - number of points per trials
%   color  - color of the windows for eegplot()
%
% Outputs:
%   eegplotarray - array defining windows which is compatible with 
%                  the function eegplot()
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eegtresh(), eeglab(), eegplot(), pop_rejepoch()

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

% ---------------------------------------------------------- 
function rejeegplot = trial2eegplot( rej, rejE, pnts, color)
	rej  = find(rej>0);
	rejE = rejE(:, rej)';
   	rejeegplot = zeros(length(rej), size(rejE,2)+5);
   	rejeegplot(:, 6:end) = rejE;
   	rejeegplot(:, 1) = (rej(:)-1)*pnts;
   	rejeegplot(:, 2) = rej(:)*pnts-1;
   	rejeegplot(:, 3:5) = ones(size(rejeegplot,1),1)*color;
return
