% eeg_readoptions() - Read EEGLAB memory options file (eeg_options) into a
%                     structure variable (opt).
%
% Usage:
%   [ header, opt ] = eeg_readoptions( filename, opt );
%
% Input:
%   filename    - [string] name of the option file
%   opt         - [struct] option structure containing backup values
%
% Outputs:
%   header      - [string] file header.
%   opt         - [struct] option structure containing an array of 3 fields
%                 varname     -> all variable names.
%                 value       -> value for each variable name
%                 description -> all description associated with each variable
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2006-
%
% See also: eeg_options(), eeg_editoptions()

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, 2006-
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

function [ header, opt ] = eeg_readoptions( filename, opt_backup );

    if nargin < 1
        help eeg_readoptions;
        return;
    end
    
    if nargin < 2
        opt_backup = [];
    end
    
    if ischar(filename)
         fid = fopen(filename, 'r');
    else fid = filename;
    end
    
    % skip header
    % -----------
    header = '';
    str = fgets( fid );
    while (str(1) == '%')
        header = [ header str];
        str = fgets( fid );
    end
    
    % read variables values and description
    % --------------------------------------
    str = fgetl( fid ); % jump a line
    index = 1;
    opt = [];
    while (str(1) ~= -1)
        if str(1) == '%'
            opt(index).description = str(3:end-1);
            opt(index).value       = [];
            opt(index).varname     = '';
        else
            [ opt(index).varname str ] = strtok(str); % variable name
            [ equal              str ] = strtok(str); % =
            [ opt(index).value   str ] = strtok(str); % value
            [ tmp                str ] = strtok(str); % ;
            [ tmp                dsc ] = strtok(str); % comment
            dsc = deblank( dsc(end:-1:1) );
            opt(index).description = deblank( dsc(end:-1:1) );
            opt(index).value       = str2num(  opt(index).value );
        end
        
        str = fgets( fid ); % jump a line
        index = index+1;
    end
    fclose(fid);

    % replace in backup structure if any
    % ----------------------------------
    if ~isempty(opt_backup)
        if ~isempty(opt)
            for index = 1:length(opt_backup)
                ind = strmatch(opt_backup(index).varname, { opt.varname }, 'exact');
                if ~isempty(ind) && ~isempty(opt_backup(index).varname)
                    opt_backup(index).value = opt(ind).value;
                end
            end
        end
        opt = opt_backup;
    end
