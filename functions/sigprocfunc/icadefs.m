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

% -----------------------------------------------------------
% ------------- START OF PATH DEFINITION --------------------
% -----------------------------------------------------------
if isunix
 ICADIR = [ '/home/scott/matlab/' ];    % <=INSERT Unix Matlab ICA dirname here
                                        %    Include trailing /
 TUTDIR = ['/home/scott/matlab/tutorial'];% <=INSERT ica tutorial dirname here

else % assume Windows                  
 ICADIR = [ 'f:\scott\matlab\' ];       % <=INSERT PC matlab ICA dirname here
                                        %    Include trailing /
 TUTDIR = ['f:\scott\matlab\tutorial']; % <=INSERT ica tutorial dirname here
end
                                       
TUTORIAL_URL = 'http://sccn.ucsd.edu/eeglab/icatutorial/'; % online version

ICABINARY = 'ica_linux.bin'; % <=INSERT name of ica executable for binica.m

% -----------------------------------------------------------
% ------------- END OF PATH DEFINITION ----------------------
% -----------------------------------------------------------

SC  =  ['binica.sc'];           % Master .sc script file for binica.m
                                % MATLAB will use first such file found
                                % in its path of script directories.
                                % Copy to pwd to alter ica defaults

BACKCOLOR  = [0.8 0.8 0.8];     % grey background color for plotting

MAXENVPLOTCHANS   = 256;  % maximum number of channels to plot in envproj.m
MAXPLOTDATACHANS  = 256;  % maximum number of channels to plot in dataplot.m
MAXPLOTDATAEPOCHS = 256;  % maximum number of epochs to plot in dataplot.m
MAXEEGPLOTCHANS   = 256;  % maximum number of channels to plot in eegplot.m
MAXTOPOPLOTCHANS  = 256;  % maximum number of channels to plot in topoplot.m

DEFAULT_SRATE = 256.0175; % default sampling rate <-- RESET TO LOCAL RATE
DEFAULT_ELOC  = 'chan.locs'; % default electrode location file for topoplot.m
DEFAULT_EPOCH = 10;       % default epoch width to plot in eegplot(s) (in sec)

if strcmp(ICADIR,'XXX/') | ~exist(ICADIR)
    fprintf('===============================================\n');
    fprintf('You have not set the ICADIR variable in icadefs.m\n');
    fprintf('===============================================\n');
    return
end
