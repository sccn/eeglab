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
if isunix
 ICADIR = [ '/home/scott/matlab/' ];   % <=INSERT Unix Matlab ICA dirname here
else                                   %    Include trailing /
 ICADIR = [ 'f:\scott\matlab\' ];      % <=INSERT PC matlab ICA dirname here
end                                    %    Include trailing /
                                       % Change this if you install the
TUTDIR = '/home/scott/matlab/tutorial';% toolbox tutorial elsewere. If
                                       % you choose not to install the tutorial
                                       % leave this as is.
TUTORIAL_URL = 'http://sccn.ucsd.edu/eeglab/icatutorial/';
ICABINARY = 'ica_linux.bin'; % <=INSERT path of ica executable for binica.m

% ------------- END OF PATH DEFINITION ----------------------
% -----------------------------------------------------------

SC  =  ['binica.sc'];           % master .sc file for binica.m

%ENVCOLORS  = [ ICADIR 'envproj.col' ]; % default color-order
ENVCOLORS  = [ 'envproj.col' ]; % default color-order
% THIS ENVCOLORS IS NOT USED BY ANY PROGRAM ANY MORE
%                                        filename for envproj.m here.

%PROJCOLORS = [ ICADIR 'white1st.col' ];% default color-order
PROJCOLORS = [ 'white1st.col' ];% default color-order
% THIS PROJCOLORS IS NOT USED BY ANY PROGRAM ANY MORE
%                                         filename for plotproj.m here.
BACKCOLOR  = [0.8 0.8 0.8];            % background color for plotting

MAXENVPLOTCHANS   = 256;  % maximum number of channels to plot in envproj.m
MAXPLOTDATACHANS  = 256;  % maximum number of channels to plot in dataplot.m
                          %         and functions that use it.
MAXPLOTDATAEPOCHS = 256;  % maximum number of epochs to plot in dataplot.m
MAXEEGPLOTCHANS   = 256;  % maximum number of channels to plot in eegplot.m
MAXTOPOPLOTCHANS  = 256;  % maximum number of channels to plot in topoplot.m
DEFAULT_ELOC  = 'chan.locs'; % default electrode location file for topoplot.m
DEFAULT_SRATE = 256.0175; % default sampling rate
DEFAULT_EPOCH = 10;       % default epoch width to plot in eegplot(s) (in sec)

if strcmp(ICADIR,'XXX/') | ~exist(ICADIR)
    fprintf('===============================================\n');
    fprintf('You have not set the ICADIR variable in icadefs.m\n');
    fprintf('===============================================\n');
    return
end
