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

function dataout = memmapdata(data, datadims);

    dataout.fileformat = 'normal';
    if ischar(data)
        if length(data) > 3
            if strcmpi('.dat', data(end-3:end))
                dataout.fileformat = 'transposed';
            end
        end
        
        % check that the file is not empty
        % --------------------------------
        a = dir(data);
        if isempty(a)
            error([ 'Data file ''' data '''not found' ]);
        elseif a(1).bytes == 0
            error([ 'Empty data file ''' data '''' ]);
        end
        
        if ~strcmpi(dataout.fileformat, 'transposed')
            dataout.data = memmapfile(data, 'writable', false, 'format', { 'single' datadims 'x' });
        else
            dataout.data = memmapfile(data, 'writable', false, 'format', { 'single' [ datadims(2:end) datadims(1) ] 'x' });            
        end
        dataout = class(dataout,'memmapdata');
    else 
        dataout = data;
    end

