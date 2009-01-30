function mri = ctf_mri_make

% ctf_mri_make - create a CTF mri struct with zero image
%
% mri = ctf_mri_make
%
% mri is a data struct similar to those returned from ctf_read_mri,
% although this one returns a zero matrix.
%
% see also: ctf_read_mri, ctf_write_mri
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

% History:  08/2003, Darren.Weber_at_radiology.ucsf.edu
%                    - adapted from an appendex to CTF document
%                    MRIConverter.pdf, which is copied at the end of this
%                    function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ver = '[$Revision: 1.1 $]';
fprintf('\nCTF_MRI_MAKE [v%s]\n',ver(12:16));  tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the header
fprintf('...creating header\n');
mri.hdr = Version_2_Header_make;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ensure the data is 16 bits (it can be 8 or 16)
mri.hdr.dataSize = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a zero image matrix, CTF mri files are always 256 x 256 x 256
fprintf('...creating zero image\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% imageOrientation, 0 = left on left, 1 = left on right
mri.hdr.imageOrientation = 0;
fprintf('...creating sagittal slices in neurological orientation (left is on the left)\n');
fprintf('...+X left to right, +Y anterior to posterior, +Z superior to inferior\n');


mri.img = zeros(256,256,256);

t=toc; fprintf('...done (%5.2f sec).\n\n',t);

return






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Version_2_Header = Version_2_Header_make,

Version_2_Header.identifierString = 'CTF_MRI_FORMAT VER 2.2          ';


Version_2_Header.imageSize           = 256; % always = 256
Version_2_Header.dataSize            = 2; % 1 or 2 (bytes), 8 or 16 bits
Version_2_Header.clippingRange       = 65536; % max. integer value of data
Version_2_Header.imageOrientation    = 1; % eg., 0 = left on left, 1 = left on right

% voxel dimensions in mm
Version_2_Header.mmPerPixel_sagittal = 1.0;
Version_2_Header.mmPerPixel_coronal  = 1.0;
Version_2_Header.mmPerPixel_axial    = 1.0;

Version_2_Header.HeadModel_Info = headModel_make;  % defined below...
Version_2_Header.Image_Info = imageInfo_make;      % defined below...

% voxel location of head origin
Version_2_Header.headOrigin_sagittal = 128;
Version_2_Header.headOrigin_coronal  = 128;
Version_2_Header.headOrigin_axial    = 128;

% euler angles to align MR to head coordinate system (angles in degrees!)
% 1. rotate in coronal plane by this angle
% 2. rotate in sagittal plane by this angle
% 3. rotate in axial plane by this angle
Version_2_Header.rotate_coronal  = 0;
Version_2_Header.rotate_sagittal = 0;
Version_2_Header.rotate_axial    = 0;

Version_2_Header.orthogonalFlag   = 0; % if set then image is orthogonal
Version_2_Header.interpolatedFlag = 0; % if set than image was interpolated

% original spacing between slices before interpolation to CTF format
Version_2_Header.originalSliceThickness = 0;

% transformation matrix head->MRI [column][row]
Version_2_Header.transformMatrix = eye(4);

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function HeadModel_Info = headModel_make,

% this function is called from Version_2_Header_make

% fid. point coordinate (in voxels)
HeadModel_Info.Nasion_Sag = 0; % nasion - sagittal
HeadModel_Info.Nasion_Cor = 0; % nasion - coronal
HeadModel_Info.Nasion_Axi = 0; % nasion - axial
HeadModel_Info.LeftEar_Sag = 0; % left ear - sagittal
HeadModel_Info.LeftEar_Cor = 0; % left ear - coronal
HeadModel_Info.LeftEar_Axi = 0; % left ear - axial
HeadModel_Info.RightEar_Sag = 0; % right ear - sagittal
HeadModel_Info.RightEar_Cor = 0; % right ear - coronal
HeadModel_Info.RightEar_Axi = 0; % right ear - axial

% fid. point coordinate (in voxels)
% HeadModel_Info.Nasion_Sag = 128; % nasion - sagittal
% HeadModel_Info.Nasion_Cor = 64; % nasion - coronal
% HeadModel_Info.Nasion_Axi = 128; % nasion - axial
% HeadModel_Info.LeftEar_Sag = 64; % left ear - sagittal
% HeadModel_Info.LeftEar_Cor = 128; % left ear - coronal
% HeadModel_Info.LeftEar_Axi = 128; % left ear - axial
% HeadModel_Info.RightEar_Sag = 192; % right ear - sagittal
% HeadModel_Info.RightEar_Cor = 128; % right ear - coronal
% HeadModel_Info.RightEar_Axi = 128; % right ear - axial

% default sphere origin
HeadModel_Info.defaultSphereX = 0.0; % sphere origin x coordinate ( in mm )
HeadModel_Info.defaultSphereY = 0.0; % sphere origin y coordinate ( in mm )
HeadModel_Info.defaultSphereZ = 0.0; % sphere origin z coordinate ( in mm )
HeadModel_Info.defaultSphereRadius = 75.0; % default sphere radius ( in mm )

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Image_Info = imageInfo_make,

% this function is called from Version_2_Header_make

Image_Info.modality         = 0; % 0 = MRI, 1 = CT, 2 = PET, 3 = SPECT, 4 = OTHER
Image_Info.manufacturerName = char(repmat(double(' '),1,64));
Image_Info.instituteName    = char(repmat(double(' '),1,64));
Image_Info.patientID        = char(repmat(double(' '),1,32));
Image_Info.dateAndTime      = char(repmat(double(' '),1,32));
Image_Info.scanType         = char(repmat(double(' '),1,32));
Image_Info.contrastAgent    = char(repmat(double(' '),1,32));
Image_Info.imagedNucleus    = char(repmat(double(' '),1,32));

Image_Info.Frequency        = 0;
Image_Info.FieldStrength    = 0;
Image_Info.EchoTime         = 0;
Image_Info.RepetitionTime   = 0;
Image_Info.InversionTime    = 0;
Image_Info.FlipAngle        = 0;
Image_Info.NoExcitations    = 0;
Image_Info.NoAcquisitions   = 0;

Image_Info.commentString    = char(repmat(double(' '),1,256));
Image_Info.forFutureUse     = char(repmat(double(' '),1,64));

return
