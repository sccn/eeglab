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
