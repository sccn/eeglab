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

function eeglab_error

    % handling errors
    % ----------------
    tmplasterr = lasterr;
    [iseeglaberror tmplasterr header] = testeeglaberror(tmplasterr);
    ft_error = false;
    if iseeglaberror
        tmp = lasterror; % Add more information to eeglab errors JRI, RMC
        if ~isempty(tmp.stack)
              tmplasterr = sprintf('%s,\n\nThis is not a bug (Error occurred in function %s() at line %d)',...
              tmplasterr, tmp.stack(1).name, tmp.stack(1).line);
        else
            tmplasterr = sprintf('%s,\n', tmplasterr);
        end
        errordlg2(tmplasterr, header);
    else
        try
            tmp = lasterror;
            if length(tmp.stack(1).name) && all(tmp.stack(1).name(1:3) == 'ft_')
                ft_error = true;
            end
            tmperr = [ 'EEGLAB error in function ' tmp.stack(1).name '() at line ' int2str(tmp.stack(1).line) ':' 10 10 lasterr ];
        catch % Matlab 5 and when the stack is empty
            tmperr = [ 'EEGLAB error:' 10 10 tmplasterr ];
        end
        if ~ft_error
            tmperr = [ tmperr 10 10 'If you think this is a bug, please submit a detailed ' 10 ...
                                    'description of how to reproduce the problem (and upload' 10 ...
                                    'a small dataset) at https://github.com/sccn/eeglab/issues' ];
        else
            tmperr = [ tmperr 10 10 'This is a problem with FIELDTRIP. The Fieldtrip version you downloaded' 10 ...
                                    'is corrupted. Please manually replace Fieldtrip with an earlier version' 10 ...
                                    'and/or email the Fieldtrip developers so they can fix the issue.' ];
        end
        errordlg2(tmperr, header);
    end
    
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
            end
        catch
        end
