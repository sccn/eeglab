% eeg_readoptions() - handle EEGLAB options. This script (not function)
%                    set the various options in the eeg_options() file.
%
% Usage:
%   [ header, varname, value, description ] = eeg_readoptions( filename );
%
% Input:
%   filename    - [string] name of the option file
%
% Outputs:
%   header      - [string] file header.
%   varname     - [cell array] all variable names.
%   value       - [cell array] value for each variable name
%   description - [cell array] all description associated with each variable
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2006-

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% read option file
% ----------------
function [ header, varname, value, description ] = eeg_readoptions( filename );
    
    if nargin < 1
        help eeg_readoptions;
        return;
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
    while (str(1) ~= -1)
        if str(1) == '%'
            description{index} = str(3:end-1);
            value{index}       = [];
            varname{index}     = '';
        else
            [ varname{index} str ] = strtok(str); % variable name
            [ equal          str ] = strtok(str); % =
            [ value{index}   str ] = strtok(str); % value
            [ tmp            str ] = strtok(str); % ;
            [ tmp            dsc ] = strtok(str); % comment
            dsc = deblank( dsc(end:-1:1) );
            description{index} = deblank( dsc(end:-1:1) );
            value{index}       = str2num( value{index} );
        end;
        
        str = fgets( fid ); % jump a line
        index = index+1;
    end;
    fclose(fid);
