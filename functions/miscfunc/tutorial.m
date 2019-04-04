% tutorial() - Bring up the ICA / electrophysiology toolbox tutorial
%              in a browser window (see docopt.m in the toolbox dir).
%              Tutorial URL: http://www.sccn.ucsd.edu/tutorial/
%              Download: See http://www.sccn.ucsd.edu/ica.html
%
% Authors: Scott Makeig & Tzyy-Ping Jung, SCCN/INC/UCSD, 12/29/00

% Copyright (C) 12/29/00 Scott Makeig & Tzyy-Ping Jung, SCCN/INC/UCSD,
% scott@sccn.ucsd.edu
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

% 01-25-02 reformated help & license -ad 

icadefs % load icadefs.m globals including TUTDIR

%TUTDIR = which('eeglab');
%TUTDIR = TUTDIR(1:findstr(TUTDIR, 'eeglab')-1);

%if exist([TUTDIR 'index.html'])
%   eval(['web file://' TUTDIR 'index.html' ]);
%else
%   fprintf('ICA Matlab Toolbox Tutorial not found in the toolbox directory.\n');
   fprintf('Opening the toolbox www site ...\n\n');
   web(TUTORIAL_URL, '-browser');
%end

