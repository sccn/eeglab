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

function tmp_fld = gettempfolder(errorflag);

if nargin < 1
    errorflag = 0;
end

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
        catch, end
    end
end
if isempty(tmp_fld) && errorflag
    error('Cannot find a temporary folder to store data files');
end
