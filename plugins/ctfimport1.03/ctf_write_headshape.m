function ctf_write_headshape(HeadShape,mri,file)

% ctf_write_headshape - write a CTF .shape file
%
% ctf_write_headshape(HeadShape,mri,file)
%
% - HeadShape is Nx3 list of coordinates in the CTF head shape coordinates.
% - mri is a struct returned by ctf_read_mri; if it is empty, this function
% will prompt with a file browser and read in the file identified.  The mri
% struct contains the mri.file field, which is used to define the name of
% the .shape and the .shape_info output files (unless the optional 'file'
% input is given).  Also, the fiducials and other information in this
% struct are required to output the .shape_info file.  The .shape and
% .shape_info files are required to import the data into CTF's MRIViewer.
%
% The *.shape file is an ascii text file in the following format:
%
% Number of Points
% x1 y1 z1
% x2 y2 z2
% .
% .
% .
% xn yn zn
%
% These vertex coordinates are contained in the input HeadShape
% matrix (Nx3).  The coordinate values must be in centimeters in
% either the voxel MRI coordinate system or the MEG Head Coordinate
% System (see ctf_read_mri for more about the coordinate system).
%
% This function assumes the MEG head coordinate system is required.
% The values in the input HeadShape matrix (Nx3) should be in
% centimeteres, with the column values representing X, Y, Z, such
% that +X is toward the nasion, +Y is left preauricular and +Z is
% superior.  The values are given with respect to the origin, which
% lies half way between the left and right preauricular fiducials.  It is
% essential that the fiducials are defined and saved in the .mri file for
% this function to operate correctly.
%
% see also ctf_write_mrishape ctf_read_mri ctf_mri2head ctf_head2mri
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
fprintf('CTF_WRITE_HEADSHAPE [v%s]\n',ver(12:16));  tic;



%--------------------------------------------------
% check the mri input variable

if exist('mri','var'),
  if isfield(mri,'file'),
    % great!
  else
    fprintf('...cannot find mri.file\n');
    [fileName, filePath] = uigetfile('*.mri', 'Pick a CTF .mri file');
    mri.file = fullfile(filePath,fileName);
    mri = ctf_read_mri(mri.file);
  end
else
  [fileName, filePath] = uigetfile('*.mri', 'Pick a CTF .mri file');
  mri.file = fullfile(filePath,fileName);
  mri = ctf_read_mri(mri.file);
end

if exist('file','var'),
    shapeFile = file;
else
    shapeFile = [mri.file,'.head.shape'];
end
shapeInfoFile = [shapeFile,'_info'];    


%--------------------------------------------------
% output the .shape file

fid = fopen(shapeFile,'w');

if isequal(fid,-1),
    S=sprintf('Could not open file: "%s"',shapeFile);
    error(S);
else
    
    fprintf('...writing to file:\n   %s\n',shapeFile);
    fprintf('...writing CTF .shape file, with head coordinates\n');
    tic;
    
    % Write vertices
    Nvertices = size(HeadShape,1);
    fprintf(fid,'%d\n',Nvertices);
    
    for v = 1:Nvertices,
      fprintf(fid,'%6.3f %6.3f %6.3f\n',HeadShape(v,1),HeadShape(v,2),HeadShape(v,3));
    end
    
    fclose(fid);
    
    t = toc;
    fprintf('...done (%6.2f sec).\n',t);
    
end

%--------------------------------------------------
% output the .shape_info file

fid = fopen(shapeInfoFile,'w');

if isequal(fid,-1),
    S=sprintf('Could not open file: "%s"',shapeInfoFile);
    error(S);
else
  
  fprintf('...writing to file:\n   %s\n',shapeInfoFile);
  fprintf('...writing CTF .shape_info file, with head coordinates\n');
  tic;
  
  [mriFilePath,mriFileName,mriFileExt] = fileparts(mri.file);
  
  nasion = [ mri.hdr.HeadModel_Info.Nasion_Sag,   mri.hdr.HeadModel_Info.Nasion_Cor,   mri.hdr.HeadModel_Info.Nasion_Axi ];
  left   = [ mri.hdr.HeadModel_Info.LeftEar_Sag,  mri.hdr.HeadModel_Info.LeftEar_Cor,  mri.hdr.HeadModel_Info.LeftEar_Axi ];
  right  = [ mri.hdr.HeadModel_Info.RightEar_Sag, mri.hdr.HeadModel_Info.RightEar_Cor, mri.hdr.HeadModel_Info.RightEar_Axi ];
  
  fprintf(fid,'// *************************************\n');
  fprintf(fid,'// 	 Headshape File Information \n');
  fprintf(fid,'// *************************************\n');
  fprintf(fid,'\n');
  fprintf(fid,'MRI_Info\n');
  fprintf(fid,'{\n');
  fprintf(fid,'	VERSION:	1.00\n');
  fprintf(fid,'\n');
  fprintf(fid,'	FILENAME:	%s\n',[mriFileName,mriFileExt]);
  fprintf(fid,'\n');
  fprintf(fid,'	// Fid. Points	Sag	Cor	Axi\n');
  fprintf(fid,'	NASION:\t\t\t%5.1f %5.1f %5.1f\n',nasion(1),nasion(2),nasion(3));
  fprintf(fid,'	LEFT_EAR:\t\t%5.1f %5.1f %5.1f\n',left(1),left(2),left(3));
  fprintf(fid,'	RIGHT_EAR:\t%5.1f %5.1f %5.1f\n',right(1),right(2),right(3));
  fprintf(fid,'\n');
  fprintf(fid,'	MM_PER_VOXEL_SAGITTAL:\t\t%10.8f\n', mri.hdr.mmPerPixel_sagittal);
  fprintf(fid,'	MM_PER_VOXEL_CORONAL:\t\t%10.8f\n',  mri.hdr.mmPerPixel_coronal);
  fprintf(fid,'	MM_PER_VOXEL_AXIAL:\t\t\t%10.8f\n',  mri.hdr.mmPerPixel_axial);
  fprintf(fid,'\n');
  fprintf(fid,'	COORDINATES:	HEAD\n');
  fprintf(fid,'\n');
  fprintf(fid,'}\n');
  
  fclose(fid);
  
  t = toc;
  fprintf('...done (%6.2f sec).\n\n',t);
  
end


return






%// *************************************
%//       Headshape File Information
%// *************************************
% 
%MRI_Info
%{
%        VERSION:        1.00
% 
%        FILENAME:       ucsf_mh_orig_axial_las2ctf.mri
% 
%        // Fid. Points          Sag     Cor     Axi
%        NASION:                 128.0   35.0    130.0
%        LEFT_EAR:               53.0    123.0   152.0
%        RIGHT_EAR:              207.0   123.0   152.0
% 
%        MM_PER_VOXEL_SAGITTAL:          1.00000000
%        MM_PER_VOXEL_CORONAL:           1.00000000
%        MM_PER_VOXEL_AXIAL:             1.00000000
% 
%        COORDINATES:    MRI
% 
%}
