% icadefs() - function to read in a set of EEGLAB system-wide (i.e. lab-wide)
%             or working directory-wide constants and preferences. Change the 
%             way these are defined in the master icadefs.m file (usually
%             in dir eeglab/functions/sigprocfunc) or make a custom copy of 
%             the icadefs.m file in a project directory. Then, calling functions 
%             that call icadefs from an EEGLAB session in that working directory 
%             will read the local copy, which may set preferences different from 
%             the system-wide copy.
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
% Revision 1.53  2007/11/30 00:26:44  arno
% now using Jason binica
%
% Revision 1.52  2007/11/01 23:06:31  arno
% change version number
%
% Revision 1.51  2007/08/23 19:24:20  nima
% _
%
% Revision 1.50  2007/08/23 19:19:22  arno
% put semicolom at the end
%
% Revision 1.49  2007/08/23 19:18:35  nima
%
% Revision 1.48  2007/08/13 16:14:07  arno
% new ICA bin
%
% Revision 1.47  2007/04/20 15:19:48  scott
% added DEFAULT_TIMLIM
%
% Revision 1.46  2007/04/20 14:51:36  scott
% undeprecated DEFAULT_SRATE
%
% Revision 1.45  2007/03/20 20:14:12  arno
% version
%
% Revision 1.44  2007/01/02 01:16:49  scott
% added debug print
%
% Revision 1.43  2006/11/12 19:15:38  arno
% version number
%
% Revision 1.42  2006/09/26 20:20:37  scott
% changed help message (1st change since '97 !)
%
% Revision 1.41  2006/09/26 20:10:37  scott
% added HZDIR to set frequency axis direction (up/down) in timef/newtimef
%
% Revision 1.40  2006/09/07 17:31:30  arno
% change revision number
%
% Revision 1.39  2006/05/27 03:23:26  toby
% temp bug fix in search of a better solution
%
% Revision 1.38  2006/05/26 01:30:52  toby
% JR says ica_linux is more recent than ica_linux2.4
%
% Revision 1.37  2006/05/12 17:53:05  arno
% changing version number
%
% Revision 1.36  2006/03/13 22:18:42  arno
% revision number
%
% Revision 1.35  2006/02/07 19:35:45  arno
% changing revision
%
% Revision 1.34  2005/09/13 21:46:00  arno
% version 4.6b
%
% Revision 1.33  2005/03/10 01:30:17  arno
% nothing
%
% Revision 1.32  2004/11/22 17:55:10  scott
% nothing
%
% Revision 1.31  2004/11/22 17:54:24  scott
% nothing
%
% Revision 1.30  2004/11/20 02:00:19  arno
% revision number
%
% Revision 1.29  2004/11/12 18:07:04  arno
% EEG_VERSION -> EEGLAB_VERSION
%
% Revision 1.28  2004/11/05 18:00:49  scott
% added variable EEG_VERSION giving current version of EEGLAB (for SCCN, currently 4.7a = alpha). -sm
%
% Revision 1.27  2004/10/15 15:17:49  scott
% added YDIR to get postive up/down standard -sm
%
% Revision 1.26  2004/07/07 22:22:53  arno
% add shrink warning
%
% Revision 1.25  2004/03/21 15:43:50  scott
% changing back - EEG not imported from workspace to functions.
%
% Revision 1.24  2004/03/21 15:34:57  scott
% debug same
%
% Revision 1.23  2004/03/21 15:33:29  scott
% made DEFAULT_ELOC = 'EEG.chanlocs'    for topoplot
%
% Revision 1.22  2003/11/25 23:01:29  arno
% changing pluginmenucolor
%
% Revision 1.21  2003/11/25 22:56:32  arno
% plugin menu color
%
% Revision 1.20  2003/07/18 14:22:44  scott
% commented, increased max channel defaults
%
% Revision 1.19  2003/02/12 23:33:41  arno
% changing default color if less or equal to 256 colors
%
% Revision 1.18  2002/11/15 17:57:33  arno
% update example tutorials
% URLs
%
% Revision 1.17  2002/11/15 15:52:57  arno
% more help for local tutorial copies
%
% Revision 1.16  2002/11/15 15:44:37  arno
% warning for tutorial path
%
% Revision 1.15  2002/11/13 17:57:03  arno
% tutorial link update
%
% Revision 1.14  2002/08/18 23:54:06  arno
% changing ica executable
%
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

EEGLAB_VERSION = '6.01b'; % EEGLAB version s=stable, b=beta, a=alpha (SCCN only)

TUTORIAL_URL = 'http://sccn.ucsd.edu/eeglab/eeglabdocs.html'; % online version
% NB: If there is a local copy of the web site, 
%     replace the line above with a line like the following: 
% TUTORIAL_URL = 'file://C:\folder\eeglabtutorial\eeglabdocs.html'; % Windows
% TUTORIAL_URL = 'file:///home/user/eeglabtutorial/eeglabdocs.html'; % Unix
% TUTORIAL_URL = 'file://::disk:folder:eeglabtutorial:eeglabdocs.html'; % Mac

%ICABINARY = '/home/duann/matlab/fmrlab4.0/ica_linux2.4'; 
ICABINARY = '/data/common/matlab/eeglab/functions/resources/ica_linux'; 
%                           % <=INSERT location of ica executable for binica.m above
%                           % If none, use []

YDIR  = 1;                  % positive potential up = 1; negative up = -1
HZDIR = 'up';               % ascending freqs = 'up'; descending = 'down' 
                            % (e.g., timef/newtimef frequency direction)
DEFAULT_SRATE = 256.0175;   % default local sampling rate 
DEFAULT_TIMLIM = [-1000 2000]; % default local epoch limits (ms)




% Set EEGLAB figure and GUI colors
% --------------------------------
if get(0, 'screendepth') <=8 % if mono or 8-bit color
    fprintf('icadefs(): Setting display parameters for mono or 8-bit color\n');
    BACKCOLOR           = [1 1 1];    % Background figure color 
    BACKEEGLABCOLOR     = [1 1 1];    % EEGLAB main window background
    GUIBUTTONCOLOR      = [1 1 1];    % Buttons colors in figures
    GUIPOPBUTTONCOLOR   = [1 1 1];    % Buttons colors in GUI windows
    GUIBACKCOLOR        = [1 1 1];    % GUI background color
    GUITEXTCOLOR        = [0 0 0];      % GUI foreground color for text    
    PLUGINMENUCOLOR     = [.5 0 .5];  % plugin menu color

else % if full color screen
    BACKCOLOR           = [.93 .96 1];    % EEGLAB Background figure color 
    BACKEEGLABCOLOR     = [.66 .76 1];    % EEGLAB main window background
    GUIBUTTONCOLOR      = BACKEEGLABCOLOR;% Buttons colors in figures
    GUIPOPBUTTONCOLOR   = BACKCOLOR;      % Buttons colors in GUI windows
    GUIBACKCOLOR        = BACKEEGLABCOLOR;% EEGLAB GUI background color <---------
    GUITEXTCOLOR        = [0 0 0.4];      % GUI foreground color for text
    PLUGINMENUCOLOR     = [.5 0 .5];      % plugin menu color
end;


% THE FOLLOWING PARAMETERS WILL BE DEPRECATED IN LATER VERSIONS
% -------------------------------------------------------------

SHRINKWARNING = 1;          % Warn user about the shrink factor in topoplot() (1/0)

MAXENVPLOTCHANS   = 264;  % maximum number of channels to plot in envproj.m
MAXPLOTDATACHANS  = 264;  % maximum number of channels to plot in dataplot.m
MAXPLOTDATAEPOCHS = 264;  % maximum number of epochs to plot in dataplot.m
MAXEEGPLOTCHANS   = 264;  % maximum number of channels to plot in eegplot.m
MAXTOPOPLOTCHANS  = 264;  % maximum number of channels to plot in topoplot.m

DEFAULT_ELOC  = 'chan.locs'; % default electrode location file for topoplot.m
DEFAULT_EPOCH = 10;       % default epoch width to plot in eegplot(s) (in sec)

SC  =  ['binica.sc'];           % Master .sc script file for binica.m
                                % MATLAB will use first such file found
                                % in its path of script directories.
                                % Copy to pwd to alter ICA defaults
