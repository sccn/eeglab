function HeadModel_Info = ctf_read_fiducials(file)

% ctf_read_fiducials - read a CTF .fid file
%
% HeadModel_Info = ctf_read_fiducials(fileName)
%
% The *.fid file contains a simple matrix of MRI fiducial
% coordinates in the voxel MRI coordinate system (see ctf_read_mri
% for more about that coordinate system).  The file is an ascii
% text file that contains scalar values, in this order:
% Nasion:    Sagittal <spaces> Coronal <spaces> Axial
% Left Ear:  Sagittal <spaces> Coronal <spaces> Axial
% Right Ear: Sagittal <spaces> Coronal <spaces> Axial
% For example,
% 128   43  158
%  41  137  196
% 210  131  199
%
% This function simply reads these values into the appropriate
% fields of the HeadModel_Info struct:
%
%>> HeadModel_Info
% 
%ans =
% 
%             Nasion_Sag: 128
%             Nasion_Cor: 43
%             Nasion_Axi: 158
%            LeftEar_Sag: 41
%            LeftEar_Cor: 137
%            LeftEar_Axi: 196
%           RightEar_Sag: 210
%           RightEar_Cor: 131
%           RightEar_Axi: 199
%         defaultSphereX: 0
%         defaultSphereY: 0
%         defaultSphereZ: 50
%    defaultSphereRadius: 75
%
% This function initialises the defaultSphere values.  The
% struct returned can be allocated into the mri.hdr.HeadModel_Info
% struct returned by ctf_read_mri, but the defaultSphere values
% may be incorrect.  The ctf_read_mri function actually reads all
% these fiducial locations directly from the .mri file, so it
% should not be necessary to use this function very often (unless
% the fiducials are missing from the .mri file, in which case it is
% also unlikely there will be an .fid file!).
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

% History:  02/2004, Darren.Weber_at_radiology.ucsf.edu
%                    - adapted from an appendex to CTF document
%                    MRIViewer.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ver = '[$Revision: 1.1 $]';
fprintf('\nCTF_READ_FIDUCIALS [v%s]\n',ver(12:16));  tic;

fprintf('...checking file input parameter\n');

if ~exist('file','var'),
  [fileName, filePath, filterIndex] = uigetfile('*.fid', 'Locate CTF .fid file');
  file = fullfile(filePath, fileName);
elseif isempty(file),
  fprintf('...file is empty\n');
  [fileName, filePath, filterIndex] = uigetfile('*.fid', 'Locate CTF .fid file');
  file = fullfile(filePath, fileName);
end
if ~exist(file,'file'),
  fprintf('...file does not exist\n');
  [fileName, filePath, filterIndex] = uigetfile('*.fid', 'Locate CTF .fid file');
  file = fullfile(filePath, fileName);
end

%-----------------------------------------------------
fprintf('...reading fiducials values.\n');

ctf_fiducials = load(file,'ascii');

HeadModel_Info.Nasion_Sag   = ctf_fiducials(1,1);
HeadModel_Info.Nasion_Cor   = ctf_fiducials(1,2);
HeadModel_Info.Nasion_Axi   = ctf_fiducials(1,3);
HeadModel_Info.LeftEar_Sag  = ctf_fiducials(2,1);
HeadModel_Info.LeftEar_Cor  = ctf_fiducials(2,2);
HeadModel_Info.LeftEar_Axi  = ctf_fiducials(2,3);
HeadModel_Info.RightEar_Sag = ctf_fiducials(3,1);
HeadModel_Info.RightEar_Cor = ctf_fiducials(3,2);
HeadModel_Info.RightEar_Axi = ctf_fiducials(3,3);

%-----------------------------------------------------
% These values are not in the .fid file, so we initialise them to
% default values here
fprintf('...initialising defaultSphere values.\n');
HeadModel_Info.defaultSphereX = 0;
HeadModel_Info.defaultSphereY = 0;
HeadModel_Info.defaultSphereZ = 0;
HeadModel_Info.defaultSphereRadius = 75;

t=toc; fprintf('...done (%5.2f sec).\n\n',t);

return
