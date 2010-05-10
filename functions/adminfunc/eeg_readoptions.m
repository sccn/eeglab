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

function [ header, opt ] = eeg_readoptions( filename, opt_backup );

    if nargin < 1
        help eeg_readoptions;
        return;
    end;
    
    if nargin < 2
        opt_backup = [];
    end;
    
    if isstr(filename)
         fid = fopen(filename, 'r');
    else fid = filename;
    end;
    
    % skip header
    % -----------
    header = '';
    str = fgets( fid );
    while (str(1) == '%')
        header = [ header str];
        str = fgets( fid );
    end;
    
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
        end;
        
        str = fgets( fid ); % jump a line
        index = index+1;
    end;
    fclose(fid);

    % replace in backup structure if any
    % ----------------------------------
    if ~isempty(opt_backup)
        if ~isempty(opt)
            for index = 1:length(opt_backup)
                ind = strmatch(opt_backup(index).varname, { opt.varname }, 'exact');
                if ~isempty(ind) & ~isempty(opt_backup(index).varname)
                    opt_backup(index).value = opt(ind).value;
                end;
            end;
        end;
        opt = opt_backup;
    end;
