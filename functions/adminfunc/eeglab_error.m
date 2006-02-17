% eeglab_error() - generate an eeglab error.
%
% Usage: >>  eeglab_error;
%
% Inputs: none, simply capture the last error generated.
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2006-
%
% see also: eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2006 Arnaud Delorme, Salk Institute, arno@salk.edu
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
% Revision 1.2  2006/01/31 19:59:42  arno
% fixing matlab 5 etc ...
%
% Revision 1.1  2006/01/31 00:11:31  arno
% Initial revision
%
% Revision 1.1  2006/01/26 21:25:24  arno
% Initial revision
%

function eeglab_error;

    % handling errords
    % ----------------
    tmplasterr = lasterr;
    try
        tmp = lasterror;
        tmperr = [ 'EEGLAB error in function ' tmp.stack(1).name '() at line ' int2str(tmp.stack(1).line) ':' 10 10 lasterr ];
    catch % Matlab 5 and when the stack is empty
        tmperr = [ 'EEGLAB error:' 10 10 tmplasterr ];
    end;
    tmperr = [ tmperr 10 10 'If you think this is a bug, please send a detailed' 10 ...
                            'description of how to reproduce the problem, with' 10 ...
                            'a (small) test dataset, to eeglab@sccn.ucsd.edu' ];
    errordlg2(tmperr, 'EEGLAB error');
    
    
