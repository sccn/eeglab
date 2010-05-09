% memmapdata() - create a memory-mapped data class
%
% Usage:
%   >> data_class = memmapdata(data);
%
% Inputs:
%   data         - input data or data file
%
% Outputs:
%   data_class    - output dataset class
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, Nov. 2008

% Copyright (C) 2008 Arnaud Delorme, SCCN, INC, UCSD
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

function dataout = memmapdata(data, datadims);

    dataout.fileformat = 'normal';
    if isstr(data)
        if length(data) > 3
            if strcmpi('.dat', data(end-3:end))
                dataout.fileformat = 'transposed';
            end;
        end;
        
        % check that the file is not empty
        % --------------------------------
        a = dir(data);
        if isempty(a)
            error([ 'Data file ''' data '''not found' ]);
        elseif a(1).bytes == 0
            error([ 'Empty data file ''' data '''' ]);
        end;
        
        if ~strcmpi(dataout.fileformat, 'transposed')
            dataout.data = memmapfile(data, 'writable', false, 'format', { 'single' datadims 'x' });
        else
            dataout.data = memmapfile(data, 'writable', false, 'format', { 'single' [ datadims(2:end) datadims(1) ] 'x' });            
        end;
        dataout = class(dataout,'memmapdata');
    else 
        dataout = data;
    end;

