function mri = ctf_read_mri(file)

% ctf_read_mri - read a CTF .mri file
%
% mri = ctf_read_mri(fileName)
%
% The CTF MRI File format used by MRIViewer consists of a binary file with
% a 1,028 byte header. The MRI data can be in 8-bit (unsigned character) or
% 16-bit (unsigned short integer) format and consists of 256 x 256 pixel
% slices, stored as 256 contiguous sagittal slices from left to right (or
% right to left if head orientation is left-on-right). Each slice is stored
% as individual pixels starting at the left, anterior, superior
% corner and scanning downwards row by row. Therefore the coronal
% position is fastest changing, axial position second fastest
% changing and sagittal position slowest changing value in the
% file, always in the positive direction for each axis (see section
% on Head Coordinate System for axis definitions). By default CTF
% MRI files have the file extension .mri 
%

% $Revision: 1.1 $ $Date: 2009-01-30 03:49:27 $

% Copyright (C) 2003  Darren L. Weber
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
fprintf('\nCTF_READ_MRI [v%s]\n',ver(12:16));  tic;

if ~exist('file','var'),
  [fileName, filePath, filterIndex] = uigetfile('*.mri', 'Locate CTF .mri file');
  file = fullfile(filePath, fileName);
elseif isempty(file),
  fprintf('...file is empty\n');
  [fileName, filePath, filterIndex] = uigetfile('*.mri', 'Locate CTF .mri file');
  file = fullfile(filePath, fileName);
end
if ~exist(file,'file'),
  fprintf('...file does not exist\n');
  [fileName, filePath, filterIndex] = uigetfile('*.mri', 'Locate CTF .mri file');
  file = fullfile(filePath, fileName);
end

mri.file = file;
fprintf('...reading %s\n', file);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open the file for reading
%    'ieee-be.l64' or 's' - IEEE floating point with big-endian byte
%                            ordering and 64 bit long data type.
[fid,message] = fopen(mri.file,'rb','s');
if fid < 0, error('cannot open file'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the file header
fprintf('...reading header ');
mri.hdr = Version_2_Header_read(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check header size, header should be 1028 bytes
header_bytes = ftell(fid);
fprintf('(%d bytes)\n',header_bytes);
if header_bytes ~= 1028,
  msg = sprintf('failed to read 1028 bytes from the header, read %d bytes',header_bytes);
  error(msg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seek beyond the header, to the beginning of the data matrix
fseek(fid,1028,'bof');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if the data is 8 or 16 bits
switch mri.hdr.dataSize,
  case 1, % we have 8 bit data
    fprintf('...reading 8 bit image data\n');
    precision = 'uchar';
  case 2, % we have 16 bit data
    fprintf('...reading 16 bit image data\n');
    precision = 'int16';
  otherwise,
    msg = sprintf('unknown mri.hdr.dataSize: %g',mri.hdr.dataSize);
    error(msg);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the image data array

% adjust for matlab version
ver = version;
ver = str2num(ver(1));
if ver < 6,
  data = fread(fid,inf,sprintf('%s',precision));
else,
  data = fread(fid,inf,sprintf('%s=>double',precision));
end
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now we have the data array, allocate it to a 3D matrix

% The CTF MRI File format used by MRIViewer consists of a binary file with a
% 1,028 byte header. The MRI data can be in 8-bit (unsigned character) or 16-bit
% (unsigned short integer) format and consists of 256 x 256 pixel slices, stored as
% 256 contiguous sagittal slices from left to right (or right to left if head orientation
% is "left-on-right"). Each slice is stored as individual pixels starting at the
% top left corner and scanning downwards row by row. Therefore the coronal
% position is fastest changing, axial position second fastest changing and sagittal
% position slowest changing value in the file, always in the positive direction for
% each axis (see section on Head Coordinate System for axis definitions). By
% default CTF MRI files have the file extension ".mri"

% MRIViewer uses these cardinal directions as axes in an internal coordinate system
% where sagittal = X, coronal = Y and axial = Z forming an additional
% right-handed coordinate system which is translated and rotated with respect to
% the Head Coordinate System and has its origin at the upper left anterior corner
% of the volume.

PixelDim = 256;
RowDim   = 256;
SliceDim = 256;

% imageOrientation, 0 = left on left, 1 = left on right
switch mri.hdr.imageOrientation,
  case 0,
    fprintf('...sagittal slices are neurological orientation (left is on the left)\n');
    fprintf('...+X left to right, +Y anterior to posterior, +Z superior to inferior\n');
  case 1,
    fprintf('...sagittal slices are radiological orientation (left is on the right)\n');
    fprintf('...+X right to left, +Y anterior to posterior, +Z superior to inferior\n');
  otherwise,
    msg = sprintf('...unknown mri.hdr.imageOrientation: %d\n',mri.hdr.imageOrientation);
    error(msg);
end

mri.img = zeros(SliceDim,PixelDim,RowDim);

n = 1;
y = 1:PixelDim; % +Y is from anterior to posterior

for x = 1:SliceDim, % +X is from left to right
  for z = 1:RowDim, % +Z is from superior to inferior
    mri.img(x,y,z) = data(n:n+(PixelDim-1));
    n = n + PixelDim;
  end
end



t=toc; fprintf('...done (%5.2f sec).\n\n',t);

return






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Version_2_Header = Version_2_Header_read(fid),

tmp = fread(fid,[1,32],'char');
tmp(find(tmp<0)) = 0;
Version_2_Header.identifierString = char( tmp );

% check the header format version is 2.2
if strmatch('CTF_MRI_FORMAT VER 2.2',Version_2_Header.identifierString),
  % OK this function should read this
else
  msg = sprintf('this function is not designed to read this format:\n%s\n',Version_2_Header.identifierString);
  warning(msg);
end

Version_2_Header.imageSize           = fread(fid,1,'short'); % always = 256
Version_2_Header.dataSize            = fread(fid,1,'short'); % 1 or 2 (bytes), 8 or 16 bits
Version_2_Header.clippingRange       = fread(fid,1,'short'); % max. integer value of data
Version_2_Header.imageOrientation    = fread(fid,1,'short'); % eg., 0 = left on left, 1 = left on right

% voxel dimensions in mm
Version_2_Header.mmPerPixel_sagittal = fread(fid,1,'float');
Version_2_Header.mmPerPixel_coronal  = fread(fid,1,'float');
Version_2_Header.mmPerPixel_axial    = fread(fid,1,'float');

Version_2_Header.HeadModel_Info = headModel(fid); % defined below...
Version_2_Header.Image_Info = imageInfo(fid);     % defined below...

% voxel location of head origin
Version_2_Header.headOrigin_sagittal = fread(fid,1,'float');
Version_2_Header.headOrigin_coronal  = fread(fid,1,'float');
Version_2_Header.headOrigin_axial    = fread(fid,1,'float');

% euler angles to align MR to head coordinate system (angles in degrees!)
% 1. rotate in coronal plane by this angle
% 2. rotate in sagittal plane by this angle
% 3. rotate in axial plane by this angle
Version_2_Header.rotate_coronal  = fread(fid,1,'float');
Version_2_Header.rotate_sagittal = fread(fid,1,'float');
Version_2_Header.rotate_axial    = fread(fid,1,'float');

Version_2_Header.orthogonalFlag   = fread(fid,1,'short'); % if set then image is orthogonal
Version_2_Header.interpolatedFlag = fread(fid,1,'short'); % if set than image was interpolated

% original spacing between slices before interpolation to CTF format
Version_2_Header.originalSliceThickness = fread(fid,1,'float');

% transformation matrix head->MRI [column][row]
Version_2_Header.transformMatrix = fread(fid,[4 4],'float')';


Version_2_Header.transformMatrixHead2MRI = Version_2_Header.transformMatrix;
Version_2_Header.transformMatrixMRI2Head = inv(Version_2_Header.transformMatrix);

% padding to 1028 bytes
%tmp = fread(fid,[1,202],'uchar'); % according to CTF manual, this should
%be 202, but it reads out to the 1028 bytes with 204.  So maybe there are a
%few bytes in the file read operations above that are missed?
tmp = fread(fid,[1,204],'uchar');
tmp = zeros(size(tmp));
Version_2_Header.unused = char( tmp ); % uchar

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function HeadModel_Info = headModel(fid),

% this function is called from Version_2_Header_read

% fid. point coordinate (in voxels)
HeadModel_Info.Nasion_Sag = fread(fid,1,'short'); % nasion - sagittal
HeadModel_Info.Nasion_Cor = fread(fid,1,'short'); % nasion - coronal
HeadModel_Info.Nasion_Axi = fread(fid,1,'short'); % nasion - axial
HeadModel_Info.LeftEar_Sag = fread(fid,1,'short'); % left ear - sagittal
HeadModel_Info.LeftEar_Cor = fread(fid,1,'short'); % left ear - coronal
HeadModel_Info.LeftEar_Axi = fread(fid,1,'short'); % left ear - axial
HeadModel_Info.RightEar_Sag = fread(fid,1,'short'); % right ear - sagittal
HeadModel_Info.RightEar_Cor = fread(fid,1,'short'); % right ear - coronal
HeadModel_Info.RightEar_Axi = fread(fid,1,'short'); % right ear - axial

fread(fid,2,'char'); % padding to 4 byte boundary - from Robert Oostenveld

% default sphere origin
HeadModel_Info.defaultSphereX = fread(fid,1,'float'); % sphere origin x coordinate ( in mm )
HeadModel_Info.defaultSphereY = fread(fid,1,'float'); % sphere origin y coordinate ( in mm )
HeadModel_Info.defaultSphereZ = fread(fid,1,'float'); % sphere origin z coordinate ( in mm )
HeadModel_Info.defaultSphereRadius = fread(fid,1,'float'); % default sphere radius ( in mm )

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Image_Info = imageInfo(fid),

% this function is called from Version_2_Header_read

Image_Info.modality         = fread(fid,1,'short'); % 0 = MRI, 1 = CT, 2 = PET, 3 = SPECT, 4 = OTHER
Image_Info.manufacturerName = char( fread(fid,[1,64],'char') );
Image_Info.instituteName    = char( fread(fid,[1,64],'char') );
Image_Info.patientID        = char( fread(fid,[1,32],'char') );
Image_Info.dateAndTime      = char( fread(fid,[1,32],'char') );
Image_Info.scanType         = char( fread(fid,[1,32],'char') );
Image_Info.contrastAgent    = char( fread(fid,[1,32],'char') );
Image_Info.imagedNucleus    = char( fread(fid,[1,32],'char') );

fread(fid,2,'char'); % padding to 4 byte boundary - from Robert Oostenveld

Image_Info.Frequency        = fread(fid,1,'float');
Image_Info.FieldStrength    = fread(fid,1,'float');
Image_Info.EchoTime         = fread(fid,1,'float');
Image_Info.RepetitionTime   = fread(fid,1,'float');
Image_Info.InversionTime    = fread(fid,1,'float');
Image_Info.FlipAngle        = fread(fid,1,'float');
Image_Info.NoExcitations    = fread(fid,1,'short');
Image_Info.NoAcquisitions   = fread(fid,1,'short');

tmp = fread(fid,[1,256],'char');
tmp = zeros(size(tmp));
Image_Info.commentString = char( tmp );

tmp = fread(fid,[1,64],'char');
tmp = zeros(size(tmp));
Image_Info.forFutureUse  = char( tmp );

return








% The CTF MRI File format used by MRIViewer consists of a binary file with a
% 1,028 byte header. The MRI data can be in 8-bit (unsigned character) or 16-bit
% (unsigned short integer) format and consists of 256 x 256 pixel slices, stored as
% 256 contiguous sagittal slices from left to right (or right to left if head orientation
% is “left-on-right”). Each slice is stored as individual pixels starting at the
% top left corner and scanning downwards row by row. Therefore the coronal
% position is fastest changing, axial position second fastest changing and sagittal
% position slowest changing value in the file, always in the positive direction for
% each axis (see section on Head Coordinate System for axis definitions). By
% default CTF MRI files have the file extension ".mri"
% The following is the C language header definitions for the CTF Format MRI
% File (current version is 2.2). Note that the Version_2_Header structure comprises
% the file header and contains other structures also defined below and is
% padded out to 1028 bytes. This is followed immediately by the voxel data
% itself.
% typedef struct Version_2_Header {
% char identifierString[32]; // "CTF_MRI_FORMAT VER 2.2"
% short imageSize; // always = 256
% short dataSize; // 1 or 2 (bytes)
% short clippingRange; // max. integer value of data
% short imageOrientation; // eg., 0 = left on left, 1 = left on right
% float mmPerPixel_sagittal; // voxel dimensions in mm
% float mmPerPixel_coronal; // voxel dimensions in mm
% float mmPerPixel_axial; // voxel dimensions in mm
% HeadModel_Info headModel; // defined below...
% Image_Info imageInfo; // defined below..
% float headOrigin_sagittal; // voxel location of head origin
% float headOrigin_coronal; // voxel location of head origin
% float headOrigin_axial; // voxel location of head origin
% // euler angles to align MR to head coordinate system (angles in degrees!)
% float rotate_coronal; // 1. rotate in coronal plane by this angle
% float rotate_sagittal; // 2. rotate in sagittal plane by this angle
% float rotate_axial; // 3. rotate in axial plane by this angle
% short orthogonalFlag; // if set then image is orthogonal
% short interpolatedFlag // if set than image was interpolated
% float originalSliceThickness // original spacing between slices before interpolation
% float transformMatrix[4][4] // transformation matrix head->MRI [column][row]
% unsigned char unused[202]; // padding to 1028 bytes
% } Version_2_Header;
% typedef struct HeadModel_Info {
% short Nasion_Sag; // fid. point coordinate (in voxels) for nasion - sagittal
% short Nasion_Cor; // nasion - coronal
% short Nasion_Axi; // nasion - axial
% short LeftEar_Sag; // left ear - sagittal
% short LeftEar_Cor; // left ear - coronal
% short LeftEar_Axi; // left ear - axial
% short RightEar_Sag; // right ear - sagittal
% short RightEar_Cor; // right ear - coronal
% short RightEar_Axi; // right ear - axial
% float defaultSphereX; // default sphere origin x coordinate ( in mm )
% float defaultSphereY; // sphere origin y coordinate ( in mm )
% float defaultSphereZ; // sphere origin z coordinate ( in mm )
% float defaultSphereRadius; // default sphere radius ( in mm )
% } HeadModel_Info;
% typedef struct Image_Info {
% short modality; // 0 = MRI, 1 = CT, 2 = PET, 3 = SPECT, 4 = OTHER
% char manufacturerName[64];
% char instituteName[64];
% char patientID[32];
% char dateAndTime[32];
% char scanType[32];
% char contrastAgent[32];
% char imagedNucleus[32];
% float Frequency;
% float FieldStrength;
% float EchoTime;
% float RepetitionTime;
% float InversionTime;
% float FlipAngle;
% short NoExcitations;
% short NoAcquisitions;
% char commentString[256];
% char forFutureUse[64];
% } Image_Info;




% From MRIViewer.pdf

% CTF MRI Coordinate System
% In contrast, the MRI coordinate system, shown in Figure 3, is defined by the
% orthogonal viewing directions typically used in radiology -- sagittal, coronal
% (also termed frontal) and axial (also termed horizontal).
% 
% MRIViewer uses these cardinal directions as axes in an internal coordinate system
% where sagittal = X, coronal = Y and axial = Z forming an additional
% right-handed coordinate system which is translated and rotated with respect to
% the Head Coordinate System and has its origin at the upper left anterior corner
% of the volume.
% Upon defining the landmark or “fiduciary” locations set by the placement of
% the head localization coils an internal representation of the MEG Head Coordinate
% System is defined using the same convention as that of the MEG acquisition
% software described above. This information combined with the voxel
% resolution in each direction is used to construct an internal transformation
% matrix (T) that is used to automatically convert any 3-dimensional location in
% the Head Coordinate System to any voxel location in the MRI. Thus, a point
% P(x, y, z) specified in the Head Coordinate System (see Figure 3) can be converted
% to the MRI voxel P'(x', y', z') by multiplication by the transformation
% matrix T where P' = TP. This immediate conversion can be seen by moving
% the cursor through any given view (hold left mouse button down and drag)
% which will display at the top of the image the voxel location and x, y, z location
% (in centimeters) in the Head Coordinate System for the current cursor position.
% In a similar manner an inverse transformation, P = T^-1 P' can be used to convert
% any MRI voxel to the corresponding location in the Head Coordinate System.  For a 
% transformation matrix of this type, T^-1 = T'. The fiduciary location information 
% is stored with the MRI data file and retrieved upon rereading
% the image file into MRIViewer.
% 
% The command-line program 'mrihead' can be used to print the transformation matrix and other MRI
% parameters to the screen for any given CTF format MRI file, based on the default fiduciary points
% last saved with the file.
