% pop_dipfit_batch() - interactively do batch scan of all ICA components
%                      with a single dipole
%                      Function deprecated. Use pop_dipfit_gridsearch()
%                      instead
%
% Usage: 
%  >> OUTEEG = pop_dipfit_batch( INEEG ); % pop up interactive window
%  >> OUTEEG = pop_dipfit_batch( INEEG, comps );
%  >> OUTEEG = pop_dipfit_batch( INEEG, comps, xgrid, ygrid, zgrid, thresh )
%
% Inputs:
%   INEEG     - input dataset
%   comps     - [integer array] component indices
%   xgrid     - [float array] x-grid. Default is 10 elements between
%               -1 and 1.
%   ygrid     - [float array] y-grid. Default is 10 elements between
%               -1 and 1.
%   zgrid     - [float array] z-grid. Default is 10 elements between
%               -1 and 1.
%   threshold - [float] threshold in percent. Default 40.
%
% Outputs:
%   OUTEEG      output dataset
%
% Authors: Robert Oostenveld, SMI/FCDC, Nijmegen 2003
%          Arnaud Delorme, SCCN, La Jolla 2003

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

% $Log: not supported by cvs2svn $
% Revision 1.15  2003/11/25 01:06:30  arno
% header edit
%
% Revision 1.14  2003/11/25 00:58:58  arno
% interface text
%
% Revision 1.13  2003/10/29 22:44:44  arno
% wording of GUI
%
% Revision 1.12  2003/10/29 16:39:27  arno
% fixing return command
%
% Revision 1.11  2003/10/29 03:13:32  arno
% constrain electrodes to sphere
%
% Revision 1.10  2003/10/08 17:35:19  arno
% fixing thresh command line bug
%
% Revision 1.9  2003/10/07 00:45:03  arno
% typo
%
% Revision 1.8  2003/09/29 22:17:07  arno
% debuging new select option
%
% Revision 1.7  2003/08/05 17:50:14  arno
% fixing dipole index for batch
%
% Revision 1.6  2003/07/01 22:23:01  arno
% nothing
%
% Revision 1.5  2003/06/30 01:46:30  arno
% fixing waitbar
%
% Revision 1.4  2003/06/30 01:43:54  arno
% *** empty log message ***
%
% Revision 1.3  2003/03/06 15:58:02  roberto
% fixed bug with channel selection of EEG data
% changed rejection threshold into percent
%
% Revision 1.1  2003/02/24 10:05:58  roberto
% Initial revision
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OUTEEG, com] = pop_dipfit_batch( varargin )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1
  help pop_dipfit_batch;
  return
else
    disp('Warning: pop_dipfit_manual is outdated. Use pop_dipfit_nonlinear instead');
    [OUTEEG, com] = pop_dipfit_gridsearch( varargin{:} );
end;
