% pop_rejepoch() - Reject pre-labeled trials in a EEG dataset. 
%                   Ask for confirmation and accept the rejection
%
% Usage:
%         >> OUTEEG = pop_rejepoch( INEEG, trialrej, confirm)
%
% Inputs:
%   INEEG      - Input dataset
%   trialrej   - Array of 0s and 1s (depicting rejected trials) (size is 
%                number of trials)
%   confirm    - Display rejections and ask for confirmation. (0=no. 1=yes;
%                default is 1).
% Outputs:
%   OUTEEG     - output dataset
%
% Example:
%          >> data2 = pop_rejepoch( EEG, [1 1 1 0 0 0] );
%             % reject the 3 first trials of a six-trial dataset
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab(), eegplot()

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.7  2003/12/04 23:25:48  arno
% simplifyaing for history
%
% Revision 1.6  2002/10/11 21:35:54  arno
% debugging function call
%
% Revision 1.5  2002/10/11 01:13:25  arno
% confirm->0 by default
%
% Revision 1.4  2002/10/10 22:31:04  arno
% debugging command call
%
% Revision 1.3  2002/08/12 02:34:11  arno
% questdlg2
%
% Revision 1.2  2002/06/25 01:57:00  arno
% debuging trial select
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 

function [EEG, com] = pop_rejepoch( EEG, tmprej, confirm);

com = '';
if nargin < 2
   help pop_rejepoch;
   return;
end;
if nargin < 3
   confirm = 1;
end;
   
if all(ismember(sort(unique(tmprej)), [0 1]))
    format0_1 = 1;
    fprintf('%d/%d trials rejected\n', sum(tmprej), EEG.trials);
else 
    format0_1 = 0;
    fprintf('%d/%d trials rejected\n', length(tmprej), EEG.trials);
end;

if confirm ~= 0
    ButtonName=questdlg2('Are you sure, you want to reject the labeled trials ?', ...
                         'Reject pre-labelled epochs -- pop_rejepoch()', 'NO', 'YES', 'YES');
    switch ButtonName,
        case 'NO', 
        	disp('Operation cancelled');
 			return;
       case 'YES',
       		disp('Compute new dataset');
    end % switch

end;

% create a new set if set_out is non nul 
% --------------------------------------
if format0_1
    EEG = pop_select( EEG, 'notrial', find(tmprej > 0));
else
    EEG = pop_select( EEG, 'notrial', tmprej);
end;

com = sprintf( '%s = pop_rejepoch( %s, %s);', inputname(1), ...
			inputname(1), vararg2str({ find(tmprej>0) 0}));		
return;
