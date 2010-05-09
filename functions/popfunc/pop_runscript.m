% pop_runscript() - Run Matlab script
%
% Usage: >> pop_runscript;
%        >> pop_runscript( filename );
%
% Input:
%   filename - [string] name of the file.
%
% Author: Arnaud Delorme, SCCN / INC / UCSD, August 2009

% Copyright (C) Arnaud Delorme, August 2009
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

function com = pop_runscript(filename);

com = [];
if nargin <1
    [filename filepath] = uigetfile('*.*', 'Please select input script -- pop_runscript()');
    
    if filename(1) == 0, return; end;

    filename = fullfile(filepath, filename);
end;

str = readtxtfile(filename);
try
    evalin('base', str);
catch
    lasterr
end;
com = sprintf('pop_runscript(''%s'');', filename);
