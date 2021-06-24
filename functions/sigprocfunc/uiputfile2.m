% uiputfile2() - same as uigputfile but remember folder location.
%
% Usage: >> uiputfile2(...)
%
% Inputs: Same as uiputfile
%
% Author: Arnaud Delorme & Hilit Serby, Scott Makeig, SCCN, UCSD, 2004
%         Thanks to input from Bas Kortmann
%
% Copyright (C) Arnaud Delorme & Hilit Serby, Scott Makeig, SCCN, UCSD, 2004

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


function varargout = uiputfile2(varargin);
mlock
persistent tmp_fld

if nargin < 1
    help uiputfile2;
    return;
end;

% remember old folder
%% search for the (mat) file which contains the latest used directory
% -------------------
olddir = pwd;
try
    eeglab_options;
    if option_rememberfolder
        if isempty(tmp_fld)
            tmp_fld = getpref('eeglab','lastdir',olddir);
        end
        cd(tmp_fld);
    end
end

%% Show the open dialog and save the latest directory to the file
% ---------------------------------------------------------------
[varargout{1} varargout{2}] = uiputfile(varargin{:});
try
    if option_rememberfolder
        if ischar(varargout{1})
            tmp_fld = varargout{2};
            tmp_fld = varargout{2};
            setpref('eeglab','lastdir',tmp_fld);
        end
    end
end
cd(olddir)


    
