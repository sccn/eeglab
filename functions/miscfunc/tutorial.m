% tutorial() - Bring up the ICA / electrophysiology toolbox tutorial
%              in a browser window (see docopt.m in the toolbox dir).
%              Tutorial URL: http://www.sccn.ucsd.edu/tutorial/
%              Download: See http://www.sccn.ucsd.edu/ica.html
%
% Authors: Scott Makeig & Tzyy-Ping Jung, SCCN/INC/UCSD, 12/29/00

% Copyright (C) 12/29/00 Scott Makeig & Tzyy-Ping Jung, SCCN/INC/UCSD,
% scott@sccn.ucsd.edu
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

