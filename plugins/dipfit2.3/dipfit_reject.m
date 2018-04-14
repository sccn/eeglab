% dipfit_reject() - remove dipole models with a poor fit
%
% Usage: 
%  >> dipout = dipfit_reject( model, reject )
%
% Inputs:
%   model	struct array with a dipole model for each component
%
% Outputs:
%   dipout	struct array with a dipole model for each component
%
% Author: Robert Oostenveld, SMI/FCDC, Nijmegen 2003

% SMI, University Aalborg, Denmark http://www.smi.auc.dk/
% FC Donders Centre, University Nijmegen, the Netherlands http://www.fcdonders.kun.nl/

% Copyright (C) 2003 Robert Oostenveld, SMI/FCDC roberto@miba.auc.dk
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dipout] = dipfit_reject(model, reject)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
   help dipfit_reject;
   return;
end;

fields = fieldnames(model);
indxtmp = strfind(fields, 'rv');
rvindx = find(not(cellfun('isempty', indxtmp))); % Do not use contain(). It is not backward-compatible.

for i=1:length(model)
  if model(i).rv>reject
    % reject this dipole model by replacing it by an empty model  
    for ifield =1:length(fields)
        if ifield ~= rvindx
            dipout(i).(fields{ifield})  = [];
        else
            dipout(i).rv  = 1;
        end
    end
  else
    for ifield =1:length(fields)
        dipout(i).(fields{ifield})  = model(i).(fields{ifield});
    end
  end
end