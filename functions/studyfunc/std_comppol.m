% std_comppol() - inverse component polarity in a component cluster
%
% Usage: [compout pol] = std_comppol(compin);
%
% Inputs:
%    compin  - component scalp maps, one per column.
%
% Outputs:
%    compout - component scalp maps some of them with inverted
%              polarities, one per column.
%    pol     - logical vector of component with inverted 
%              polarities (same length as the number of rows in 
%              compin)
%
% Author: Arnaud Delorme & Hilit Serby, SCCN, INC, UCSD, 2004

% Copyright (C) 2004 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [compin, pol] = std_comppol(compin);

if nargin < 1
    help std_comppol;
    return;
end;

% remove the NaN
% --------------
for index = 1:size(compin,2)
    compin(isnan(compin(:,index)),:) =[];
end;

% run several iterations
% ----------------------
pol     = ones(1,size(compin,2));
for repeat=1:3
    compave = mean(compin,2);
    for index = 1:size(compin,2)
        
        % remove diagonal and put 0 and 1
        % -------------------------------
        if ~all(compin(:,index) == 0)
             r = corrcoef(compave, compin(:,index) );
        else r = zeros(2,2);
        end;
        
        % invert component polarities
        % ---------------------------
        if r(2) < 0
            compin(:,index) = -compin(:,index);
            pol(index)      = -pol(index);
        end;
    end;
end;
