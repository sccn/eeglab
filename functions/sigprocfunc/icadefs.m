% icadefs() - file of default filenames and constants to source in the ICA /ERP
%             package functions.  Insert local dir reference below. 
%
% Note: Edit this file to change local directories under Unix and Windows 
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 05-20-97 

% Copyright (C) 05-20-97 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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
% Revision 1.13  2002/08/14 16:58:59  arno
% update for new release
%
% Revision 1.12  2002/08/12 14:50:47  arno
% [6~color
%
% Revision 1.11  2002/08/12 14:40:25  arno
% color
%
% Revision 1.10  2002/08/12 14:36:02  arno
% color
%
% Revision 1.9  2002/08/12 14:30:17  arno
% color
%
% Revision 1.8  2002/08/12 01:19:24  arno
% change colors
%
% Revision 1.7  2002/08/11 22:46:44  arno
% backcolor
%
% Revision 1.6  2002/08/11 19:26:11  arno
% editing
%
% Revision 1.5  2002/07/31 23:30:10  arno
% updating
% ,
%
% Revision 1.4  2002/06/25 02:36:27  scott
% removed outdated defs, clarified comments and layout -sm
%
% Revision 1.3  2002/05/01 18:22:55  arno
% making binica available from everywhere
%
% Revision 1.2  2002/04/09 02:17:27  arno
% adding comments for unused variables
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% Version history:
% 08-04-00  added ICA and SC -sm
% 01-24-02  added directory check -ad
% 01-25-02  reformated help & license -ad 

% ------------------------------------------------------
% -------------- EEGLAB DEFINITION (V 4.0) -------------
% ------------------------------------------------------

TUTORIAL_URL = 'http://sccn.ucsd.edu/eeglab/eeglab.html'; % online version

ICABINARY = 'ica_linux2.4'; % <=INSERT name of ica executable for binica.m

% COLORS
% ------

BACKCOLOR           = [.93 .96 1];    % Background figure color 
BACKEEGLABCOLOR     = [.66 .76 1];    % EEGLAB main window background
GUIBUTTONCOLOR      = [.66 .76 1];    % Buttons colors in figures
GUIPOPBUTTONCOLOR   = [.93 .96 1];    % Buttons colors in GUI windows
GUIBACKCOLOR        = [.66 .76 1];    % GUI background color
GUITEXTCOLOR        = [0 0 0.4];      % GUI foreground color for text

% THE FOLLOWING PARAMETERS WILL BE DEPRECATED IN LATER VERSIONS
% -------------------------------------------------------------

MAXENVPLOTCHANS   = 256;  % maximum number of channels to plot in envproj.m
MAXPLOTDATACHANS  = 256;  % maximum number of channels to plot in dataplot.m
MAXPLOTDATAEPOCHS = 256;  % maximum number of epochs to plot in dataplot.m
MAXEEGPLOTCHANS   = 256;  % maximum number of channels to plot in eegplot.m
MAXTOPOPLOTCHANS  = 256;  % maximum number of channels to plot in topoplot.m

DEFAULT_SRATE = 256.0175; % default sampling rate <-- RESET TO LOCAL RATE
DEFAULT_ELOC  = 'chan.locs'; % default electrode location file for topoplot.m
DEFAULT_EPOCH = 10;       % default epoch width to plot in eegplot(s) (in sec)

SC  =  ['binica.sc'];           % Master .sc script file for binica.m
                                % MATLAB will use first such file found
                                % in its path of script directories.
                                % Copy to pwd to alter ica defaults
