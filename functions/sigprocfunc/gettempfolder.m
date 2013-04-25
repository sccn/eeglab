% gettempfolder() - return the temporary folder
%
% Usage: >> folder = gettempfolder;
%
% Output: a string containing the folder if a temporary folder can be found. 
%         Empty if the folder cannot be found.
%
% Author: Arnaud Delorme, SCCN, UCSD, 2012
%
% Copyright (C) Arnaud Delorme

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

function tmp_fld = gettempfolder(errorflag);

if nargin < 1
    errorflag = 0;
end;

tmp_fld = getenv('TEMP');
if isempty(tmp_fld) && isunix
    if is_sccn && exist('/var/tmp')
        tmp_fld = '/var/tmp';
    elseif exist('/tmp') == 7
        tmp_fld = '/tmp';
    else
        try
            mkdir('/tmp');
            tmp_fld = '/tmp';
        catch, end;
    end;
end;
if isempty(tmp_fld) && errorflag
    error('Cannot find a temporary folder to store data files');
end;