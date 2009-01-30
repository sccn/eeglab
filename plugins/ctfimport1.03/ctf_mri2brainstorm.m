function [Comment, Cube, PCS, Segment, Voxsize] = ctf_mri2brainstorm(mri)

% ctf_mri2brainstorm - create subjectimage variables from ctf .mri file
%
% [Comment, Cube, PCS, Segment, Voxsize] = ctf_mri2brainstorm(mri)
%
% The purpose of this function is to convert a CTF .mri volume into 
% brainstorm subjectimage variables.  The returned variables can be saved
% into *_subjectimage.mat using the matlab save command.
%
% The mri struct input for this function is returned by ctf_read_mri.
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

% History:  05/2004, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '[$Revision: 1.1 $]';
fprintf('\nCTF_MRI2BRAINSTORM [v%s]\n',ver(12:16));  tic;


Comment = 'CTF';

% CTF .mri files are always 256x256x256 voxels, 1mm^3
Voxsize = [1 1 1];

% Initialize the Segment variable; As of 05/2004, the codes for this field
% are not clear to me.  They are defined in David Shattuck's
% BrainSuite MRI processing software.
Segment = [];

% the ctf mri.img orientation is
% +X right
% +Y posterior
% +Z inferior

% the brainstorm Cube orientation is
% +X right
% +Y anterior
% +Z superior

temp = mri.img;
temp = flipdim(temp,2);
Cube = flipdim(temp,3);

% -----
% Patient Coordinate System

PCS.Comment = 'CTF';

PCS.R = zeros(3,3);  % not clear how to get this
PCS.t = zeros(3,1);  % not clear how to get this

PCS.PCSFiducial = zeros(3,4); % not clear how to get this

PCS.CubeFiducial = zeros(3,4);
PCS.CubeFiducial(1,1) =       mri.hdr.HeadModel_Info.Nasion_Sag;
PCS.CubeFiducial(2,1) = 256 - mri.hdr.HeadModel_Info.Nasion_Cor;
PCS.CubeFiducial(3,1) = 256 - mri.hdr.HeadModel_Info.Nasion_Axi;

PCS.CubeFiducial(1,2) =       mri.hdr.HeadModel_Info.LeftEar_Sag;
PCS.CubeFiducial(2,2) = 256 - mri.hdr.HeadModel_Info.LeftEar_Cor;
PCS.CubeFiducial(3,2) = 256 - mri.hdr.HeadModel_Info.LeftEar_Axi;

PCS.CubeFiducial(1,3) =       mri.hdr.HeadModel_Info.RightEar_Sag;
PCS.CubeFiducial(2,3) = 256 - mri.hdr.HeadModel_Info.RightEar_Cor;
PCS.CubeFiducial(3,3) = 256 - mri.hdr.HeadModel_Info.RightEar_Axi;

% set the origin at half way between the left and right ears
PCS.CubeFiducial(1,4) = PCS.CubeFiducial(1,3) - PCS.CubeFiducial(1,2);
PCS.CubeFiducial(2,4) = PCS.CubeFiducial(2,3) - PCS.CubeFiducial(2,2);
PCS.CubeFiducial(3,4) = PCS.CubeFiducial(3,2) + (PCS.CubeFiducial(3,3) - PCS.CubeFiducial(3,2));

PCS.FiducialName = {'NAS'  'LPA'  'RPA'  'Origin'};

PCS.PERM = [1 2 3]; % not sure what this is
PCS.FLP  = [1 1 0]; % not sure what this is

t=toc; fprintf('...done (%5.2f sec).\n\n',t);

return
