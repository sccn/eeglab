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

% $Log: not supported by cvs2svn $
% Revision 1.28  2004/08/12 22:03:40  arno
% same
%
% Revision 1.27  2004/08/12 22:02:20  arno
% same
%
% Revision 1.26  2004/08/12 22:00:45  arno
% topoplot with dipole
%
% Revision 1.25  2004/08/12 21:52:59  arno
% do not crash
%
% Revision 1.24  2004/03/26 01:33:29  arno
% plot only active dipoles
%
% Revision 1.23  2003/12/04 17:09:07  arno
% num2str -> str2num
%
% Revision 1.22  2003/12/04 16:05:58  arno
% setting values when modifying gui
%
% Revision 1.21  2003/10/31 16:45:52  arno
% worind from Scott's feedback
%
% Revision 1.20  2003/10/30 02:39:14  arno
% gui typo
%
% Revision 1.19  2003/10/29 23:21:01  arno
% more checking for component index
%
% Revision 1.18  2003/10/29 23:03:37  arno
% [Adecrease window width
%
% Revision 1.17  2003/10/29 16:08:00  arno
% removing debug msg
%
% Revision 1.16  2003/10/29 03:12:11  arno
% contrain electrode to sphere
%
% Revision 1.15  2003/09/30 15:42:34  roberto
% minor bug fix, related to the change in dipole constraint handling
%
% Revision 1.14  2003/09/12 08:43:18  roberto
% changed symmetry constraint into optional input argument
%
% Revision 1.12  2003/07/01 23:49:08  arno
% implementing dipole flipping
%
% Revision 1.11  2003/06/30 02:11:54  arno
% *** empty log message ***
%
% Revision 1.10  2003/06/16 10:10:11  roberto
% added interruptible dialog (gui) for non-linear fit
%
% Revision 1.9  2003/06/13 16:48:57  arno
% remove chanlocs conversion
%
% Revision 1.8  2003/06/13 01:26:01  arno
% automatic conversion of channel location files
%
% Revision 1.7  2003/06/13 01:17:49  arno
% adding plotting button
%
% Revision 1.6  2003/03/12 10:32:32  roberto
% fixed dialog title
%
% Revision 1.5  2003/03/06 15:58:21  roberto
% fixed bug with channel selection of EEG data
%
% Revision 1.3  2003/03/03 16:52:27  roberto
% modified for posxyz/momxyz instead of dip.pos/dip.mom
% changed large listbox into edit field
%
% Revision 1.2  2003/02/28 23:03:14  arno
% making the listbox to scroll component
%
% Revision 1.1  2003/02/24 10:06:15  roberto
% Initial revision
%

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
