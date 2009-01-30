function [ctf] = ctf_read_badchannels(folder,ctf);

% ctf_read_badchannels - read BadChannels file from a CTF .ds folder
%
% [ctf] = ctf_read_badchannels(folder,ctf)
%
% This function reads the BadChannels text file in a .ds folder.
%
% The BadChannels file is a text file that contains the channel names of
% any bad channels.  The channel names are assumed to be listed in a single
% column, one channel per line.
%
% INPUTS
%
% If you do not wish to specify an input option, use [], but keep the order
% of the input options as above.  Only specify as many input options as
% required.  With no input options, the function will prompt for a folder
% and read the BadChannels file (if it exists).
%
% folder - the directory of the .ds data set to read.  By
% default, a gui prompts for the folder.
%
% ctf - a struct with folder, setup, sensor and data fields.
%
% OUTPUTS
%
% ctf.sensor.bad - a string cell array of channel labels
%
%      <>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %
%      <                                                       > %  
%      <                      DISCLAIMER:                      > %
%      <                                                       > %
%      < THIS PROGRAM IS INTENDED FOR RESEARCH PURPOSES ONLY.  > %
%      < THIS PROGRAM IS IN NO WAY INTENDED FOR CLINICAL OR    > %
%      <                     OFFICIAL USE.                     > %
%      <                                                       > %
%      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<> %
%


% $Revision: 1.1 $ $Date: 2009-01-30 03:49:27 $

% Copyright (C) 2005  Darren L. Weber
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% Modified: 01/2005, Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-----------------------------------------------
% ensure we have the folder

if exist('folder','var'),
  if exist('ctf','var'),
    ctf = ctf_folder(folder,ctf);
  else
    ctf = ctf_folder(folder);
  end
else
  if exist('ctf','var'),
    ctf = ctf_folder([],ctf);
  else
    ctf = ctf_folder;
  end
end


%--------------------------------------------------------------
ver = '$Revision: 1.1 $ $Date: 2009-01-30 03:49:27 $';
fprintf('\nCTF_READ_BADCHANNELS [v %s]\n',ver(11:15)); tic;


%----------------------------------------------------------------
% read the file

[folderPath,folderName,folderExt] = fileparts(ctf.folder);
BadChannelsFile = findBadChannelsFile( ctf.folder );

ctf.sensor.bad = textread(BadChannelsFile,'%s','delimiter','\n');

t = toc; fprintf('...done (%6.2f sec)\n\n',t);

return



% -------------------------------------------------------
function BadChannelsFile = findBadChannelsFile( folder )

filename = dir([ folder filesep 'BadChannels' ]);
if isempty(filename)
    filename = dir([ folder filesep 'badchannels' ]);
end;

if isempty(filename)
    error('No BadChannels file in selected folder');
else
    BadChannelsFile = [ folder filesep filename.name ];
end;

return
