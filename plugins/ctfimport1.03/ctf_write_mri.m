function ctf_write_mri(mri, fileName, force)

% ctf_write_mri - write a CTF .mri file
%
% ctf_write_mri(mri,fileName,force)
%
% mri is a data struct returned from ctf_read_mri.  It may contain
% mri.file, in which case, you do not need the fileName input here.  If the
% file exists, you are prompted for a new file name, unless force = 1.
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
fprintf('\nCTF_WRITE_MRI [v%s]\n',ver(12:16));  tic;

if ~exist('mri','var'), error('no input mri data struct'); end
if isempty(mri), error('empty input mri data struct'); end

if ~exist('force','var'), force = 0; end % don't overwrite
if isempty(force), force = 0; end

if ~exist('fileName','var'),
    if isfield(mri,'file'),
        file = mri.file;
    else
        [fileName, filePath, filterIndex] = uigetfile('*.mri', 'Locate CTF .mri file');
        file = fullfile(filePath, fileName);
    end
else
    [filePath, fileName, fileExt] = fileparts(fileName);
    file = fullfile(filePath, [fileName,'.mri']);
end

if isempty(file),
  error('...file is empty\n');
end

if exist(file,'file'),
    if force,
        fprintf('...file already exists, overwriting it.\n');
    else
        fprintf('...file already exists\n');
        [fileName, pathName] = uiputfile('*.mri', 'Specify CTF .mri file to write');
        file = fullfile(filePath, [fileName,'.mri']);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% open the file for writing
%    'ieee-be.l64' or 's' - IEEE floating point with big-endian byte
%                            ordering and 64 bit long data type.
[fid,message] = fopen(file,'wb','s');
if fid < 0, error('cannot open file for writing'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write the file header
fprintf('...writing header ');
Version_2_Header_write(fid, mri.hdr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check header size, header should be 1028 bytes
header_bytes = ftell(fid);
fprintf('(wrote %d bytes)\n',header_bytes);
if header_bytes ~= 1028,
  msg = sprintf('failed to write 1028 bytes into the header, wrote %d bytes',header_bytes);
  error(msg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seek beyond the header, to the beginning of the data matrix
%fseek(fid,1028,'bof');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if the data is 8 or 16 bits
switch mri.hdr.dataSize,
  case 1, % we have 8 bit data
    fprintf('...writing unsigned char (8 bit) image data\n');
    precision = 'uchar';
  case 2, % we have 16 bit data
    fprintf('...writing unsigned short (16 bit) image data\n');
    precision = 'ushort';
  otherwise,
    msg = sprintf('unknown mri.hdr.dataSize: %g',mri.hdr.dataSize);
    error(msg);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now have the data array in a 3D matrix, we have to write it out in the
% correct byte order

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

% output into sagittal slices, with the fastest moving index being Y,
% from anterior to posterior, then Z, from superior to inferior, then X,
% from left to right (or vice versa; depending on input mri struct).

n = 1;
y = 1:PixelDim; % +Y is from anterior to posterior

for x = 1:SliceDim, % +X is from left to right (or vice versa)
  for z = 1:RowDim, % +Z is from superior to inferior
      
      count = fwrite(fid,mri.img(x,y,z),precision);
      
      if count ~= PixelDim,
          error('failed to output 256 data points');
      end
      
  end
end

fclose(fid);

t=toc; fprintf('...done (%5.2f sec).\n\n',t);

return






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Version_2_Header_write(fid,Version_2_Header),

identifierString = sprintf('%-32s', Version_2_Header.identifierString);
if length(identifierString) < 32,
    paddingN = 32 - length(identifierString);
    padding = char(repmat(double(' '),1,paddingN));
    identifierString = [identifierString,padding];
end
    
fwrite(fid,identifierString(1:32),'char');

fwrite(fid,Version_2_Header.imageSize           ,'short'); % always = 256
fwrite(fid,Version_2_Header.dataSize            ,'short'); % 1 or 2 (bytes), 8 or 16 bits
fwrite(fid,Version_2_Header.clippingRange       ,'short'); % max. integer value of data
fwrite(fid,Version_2_Header.imageOrientation    ,'short'); % eg., 0 = left on left, 1 = left on right

% voxel dimensions in mm
fwrite(fid,Version_2_Header.mmPerPixel_sagittal ,'float');
fwrite(fid,Version_2_Header.mmPerPixel_coronal  ,'float');
fwrite(fid,Version_2_Header.mmPerPixel_axial    ,'float');

headModel_write(fid,Version_2_Header.HeadModel_Info); % defined below...
imageInfo_write(fid,Version_2_Header.Image_Info);     % defined below...

% voxel location of head origin
fwrite(fid,Version_2_Header.headOrigin_sagittal ,'float');
fwrite(fid,Version_2_Header.headOrigin_coronal  ,'float');
fwrite(fid,Version_2_Header.headOrigin_axial    ,'float');

% euler angles to align MR to head coordinate system (angles in degrees!)
% 1. rotate in coronal plane by this angle
% 2. rotate in sagittal plane by this angle
% 3. rotate in axial plane by this angle
fwrite(fid,Version_2_Header.rotate_coronal  ,'float');
fwrite(fid,Version_2_Header.rotate_sagittal ,'float');
fwrite(fid,Version_2_Header.rotate_axial    ,'float');

fwrite(fid,Version_2_Header.orthogonalFlag   ,'short'); % if set then image is orthogonal
fwrite(fid,Version_2_Header.interpolatedFlag ,'short'); % if set than image was interpolated

% original spacing between slices before interpolation to CTF format
fwrite(fid,Version_2_Header.originalSliceThickness ,'float');

% transformation matrix head->MRI [column][row]
fwrite(fid,Version_2_Header.transformMatrix' ,'float')';

% padding to 1028 bytes
% according to CTF manual, this should
%be 202, but it works out to the 1028 bytes with 204.
spaces = char(repmat(double(' '),1,204));
fwrite(fid,spaces,'uchar');

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function headModel_write(fid,HeadModel_Info),

% this function is called from Version_2_Header_read

% fid. point coordinate (in voxels)
fwrite(fid,HeadModel_Info.Nasion_Sag ,'short'); % nasion - sagittal
fwrite(fid,HeadModel_Info.Nasion_Cor ,'short'); % nasion - coronal
fwrite(fid,HeadModel_Info.Nasion_Axi ,'short'); % nasion - axial
fwrite(fid,HeadModel_Info.LeftEar_Sag ,'short'); % left ear - sagittal
fwrite(fid,HeadModel_Info.LeftEar_Cor ,'short'); % left ear - coronal
fwrite(fid,HeadModel_Info.LeftEar_Axi ,'short'); % left ear - axial
fwrite(fid,HeadModel_Info.RightEar_Sag ,'short'); % right ear - sagittal
fwrite(fid,HeadModel_Info.RightEar_Cor ,'short'); % right ear - coronal
fwrite(fid,HeadModel_Info.RightEar_Axi ,'short'); % right ear - axial

fwrite(fid,'  ','char'); % padding to 4 byte boundary - from Robert Oostenveld

% default sphere origin
fwrite(fid,HeadModel_Info.defaultSphereX ,'float'); % sphere origin x coordinate ( in mm )
fwrite(fid,HeadModel_Info.defaultSphereY ,'float'); % sphere origin y coordinate ( in mm )
fwrite(fid,HeadModel_Info.defaultSphereZ ,'float'); % sphere origin z coordinate ( in mm )
fwrite(fid,HeadModel_Info.defaultSphereRadius ,'float'); % default sphere radius ( in mm )

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imageInfo_write(fid,Image_Info),

% this function is called from Version_2_Header_read

fwrite(fid,Image_Info.modality,'short'); % 0 = MRI, 1 = CT, 2 = PET, 3 = SPECT, 4 = OTHER


% check the length of char variables

% Image_Info.manufacturerName = char( fread(fid,[1,64],'char') );
% Image_Info.instituteName    = char( fread(fid,[1,64],'char') );
% Image_Info.patientID        = char( fread(fid,[1,32],'char') );
% Image_Info.dateAndTime      = char( fread(fid,[1,32],'char') );
% Image_Info.scanType         = char( fread(fid,[1,32],'char') );
% Image_Info.contrastAgent    = char( fread(fid,[1,32],'char') );
% Image_Info.imagedNucleus    = char( fread(fid,[1,32],'char') );

manufacturerName = sprintf('%-64s', Image_Info.manufacturerName);
if length(manufacturerName) < 64,
    paddingN = 64 - length(manufacturerName);
    padding = char(repmat(double(' '),1,paddingN));
    manufacturerName = [manufacturerName,padding];
end

instituteName = sprintf('%-64s', Image_Info.instituteName);
if length(instituteName) < 64,
    paddingN = 64 - length(instituteName);
    padding = char(repmat(double(' '),1,paddingN));
    instituteName = [instituteName,padding];
end

patientID = sprintf('%-32s', Image_Info.patientID);
if length(patientID) < 32,
    paddingN = 32 - length(patientID);
    padding = char(repmat(double(' '),1,paddingN));
    patientID = [patientID,padding];
end

dateAndTime = sprintf('%-32s', Image_Info.dateAndTime);
if length(dateAndTime) < 32,
    paddingN = 32 - length(dateAndTime);
    padding = char(repmat(double(' '),1,paddingN));
    dateAndTime = [dateAndTime,padding];
end

scanType = sprintf('%-32s', Image_Info.scanType);
if length(scanType) < 32,
    paddingN = 32 - length(scanType);
    padding = char(repmat(double(' '),1,paddingN));
    scanType = [scanType,padding];
end

contrastAgent = sprintf('%-32s', Image_Info.contrastAgent);
if length(contrastAgent) < 32,
    paddingN = 32 - length(contrastAgent);
    padding = char(repmat(double(' '),1,paddingN));
    contrastAgent = [contrastAgent,padding];
end

imagedNucleus = sprintf('%-32s', Image_Info.imagedNucleus);
if length(imagedNucleus) < 32,
    paddingN = 32 - length(imagedNucleus);
    padding = char(repmat(double(' '),1,paddingN));
    imagedNucleus = [imagedNucleus,padding];
end

% output these char variables

fwrite(fid,manufacturerName(1:64),'char');
fwrite(fid,instituteName(1:64),'char');
fwrite(fid,patientID(1:32),'char');
fwrite(fid,dateAndTime(1:32),'char');
fwrite(fid,scanType(1:32),'char');
fwrite(fid,contrastAgent(1:32),'char');
fwrite(fid,imagedNucleus(1:32),'char');

fwrite(fid,'  ','char'); % padding to 4 byte boundary - from Robert Oostenveld

fwrite(fid,Image_Info.Frequency        ,'float');
fwrite(fid,Image_Info.FieldStrength    ,'float');
fwrite(fid,Image_Info.EchoTime         ,'float');
fwrite(fid,Image_Info.RepetitionTime   ,'float');
fwrite(fid,Image_Info.InversionTime    ,'float');
fwrite(fid,Image_Info.FlipAngle        ,'float');
fwrite(fid,Image_Info.NoExcitations    ,'short');
fwrite(fid,Image_Info.NoAcquisitions   ,'short');


commentString = sprintf('%-256s', Image_Info.commentString);
if length(commentString) < 256,
    paddingN = 256 - length(commentString);
    padding = char(repmat(double(' '),1,paddingN));
    commentString = [commentString,padding];
end

forFutureUse = sprintf('%-64s', Image_Info.forFutureUse);
if length(forFutureUse) < 64,
    paddingN = 64 - length(forFutureUse);
    padding = char(repmat(double(' '),1,paddingN));
    forFutureUse = [forFutureUse,padding];
end

fwrite(fid,commentString(1:256),'char');
fwrite(fid,forFutureUse(1:64),'char');

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
