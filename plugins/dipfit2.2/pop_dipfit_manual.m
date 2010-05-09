% pop_dipfit_manual() - interactively do dipole fit of selected ICA components
%                       Function deprecated. Use pop_dipfit_nonlinear()
%                       instead
% Usage: 
%  >> OUTEEG = pop_dipfit_manual( INEEG )
%
% Inputs:
%   INEEG       input dataset
%
% Outputs:
%   OUTEEG      output dataset
%
% Author: Robert Oostenveld, SMI/FCDC, Nijmegen 2003
%         Arnaud Delorme, SCCN, La Jolla 2003

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
function [OUTEEG, com] = pop_dipfit_manual( varargin )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1
  help pop_dipfit_manual;
  return
else
    disp('Warning: pop_dipfit_manual is outdated. Use pop_dipfit_nonlinear instead');
    [OUTEEG, com] = pop_dipfit_nonlinear( varargin{:} );
end;
