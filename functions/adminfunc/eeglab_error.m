% eeglab_error() - generate an eeglab error.
%
% Usage: >>  eeglab_error;
%
% Inputs: none, simply capture the last error generated.
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2006-
%
% see also: eeglab()

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

function eeglab_error

    % handling errors
    % ----------------
    tmplasterr = lasterr;
    [iseeglaberror tmplasterr header] = testeeglaberror(tmplasterr);
    ft_error = false;
    if iseeglaberror
        tmp = lasterror; % Add more information to eeglab errors JRI, RMC
        tmplasterr = sprintf('%s,\n\n  (Error occurred in function %s() at line %d)',...
          tmplasterr, tmp.stack(1).name, tmp.stack(1).line);
        errordlg2(tmplasterr, header);
    else
        try
            tmp = lasterror;
            if length(tmp.stack(1).name) && all(tmp.stack(1).name(1:3) == 'ft_')
                ft_error = true;
            end;
            tmperr = [ 'EEGLAB error in function ' tmp.stack(1).name '() at line ' int2str(tmp.stack(1).line) ':' 10 10 lasterr ];
        catch % Matlab 5 and when the stack is empty
            tmperr = [ 'EEGLAB error:' 10 10 tmplasterr ];
        end;
        if ~ft_error
            tmperr = [ tmperr 10 10 'If you think this is a bug, please submit a detailed ' 10 ...
                                    'description of how to reproduce the problem (and upload' 10 ...
                                    'a small dataset) at http://sccn.ucsd.edu/eeglab/bugzilla' ];
        else
            tmperr = [ tmperr 10 10 'This is a problem with FIELDTRIP. The Fieldtrip version you downloaded' 10 ...
                                    'is corrupted. Please manually replace Fieldtrip with an earlier version' 10 ...
                                    'and/or email the Fieldtrip developpers so they can fix the issue.' ];
        end;
        errordlg2(tmperr, header);
    end;
    
    function [ val str header ] = testeeglaberror(str)
        header = 'EEGLAB error';
        val = 0;
        try,
            if strcmp(str(1:5), 'Error'), 
                val = 1; 
                str = str(12:end); %Corrected extraction of function name, occurs after 'Error using ' JRI, RMC
                indendofline = find(str == 10);
                funcname = str(1:indendofline(1)-1);
                str = str(indendofline(1)+1:end);
                header = [ funcname ' error' ];
            end;
        catch
        end;