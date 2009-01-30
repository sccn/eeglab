function [ctf] = ctf_folder(folder,ctf);

% ctf_folder - get and check CTF .ds folder name
%
% [ctf] = ctf_folder( [folder], [ctf] );
% 
% folder:  The .ds directory of the dataset.  It should be a complete path
% or given relative to the current working directory (given by pwd).  The
% returned value will ensure the complete path is identified.  If this
% argument is empty or not given, a graphical prompt for the folder
% appears.
%
% eg,
%     ctf = ctf_folder;
% 
% ctf.folder is returned (as a complete path).
%
%      <>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %
%      <                                                      > %  
%      <                    DISCLAIMER:                       > %
%      <                                                      > %
%      < THIS PROGRAM IS INTENDED FOR RESEARCH PURPOSES ONLY. > %
%      < THIS PROGRAM IS IN NO WAY INTENDED FOR CLINICAL OR   > %
%      <                    OFFICIAL USE.                     > %
%      <                                                      > %
%      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<> %
%


% $Revision: 1.1 $ $Date: 2009-01-30 03:49:27 $

% Copyright (C) 2004  Darren L. Weber
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

% Modified: 01/2004, Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~exist('folder','var'), folder = []; end
if ~exist('ctf','var'), ctf = []; end

if isempty(folder),
    if isfield(ctf,'folder'),
        folder = ctf.folder;
    else
        folder = [];
    end
end

if exist(folder) ~= 7,
  fprintf('...folder inputs are invalid\n');
  folder = getfolder;
end

ctf.folder = folder;

% ensure we get the folder path
current_dir = pwd;
cd(ctf.folder);
cd ..
folderPath = pwd;
cd(current_dir);

% check whether the folder path is in the folder already
[path,file] = fileparts(ctf.folder);

% if findstr(folderPath,ctf.folder),
  % OK the path is already in the folder
% else
if isempty(path)
  % Add the folderPath
  ctf.folder = fullfile(folderPath,ctf.folder);
end


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function folder = getfolder,

% This does not work on linux systems (only windows)
folder = uigetdir(pwd,'locate CTF .ds folder');

if ~folder,
    % try to get it this way
    [filename, folder, filterindex] = uigetfile('*.res4', 'Pick a .res4 file');
    if ~folder,
      error('invalid folder');
    end
end

return
