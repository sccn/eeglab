function HeadVertices = ctf_mri2head(MRIVertices,mri)

% ctf_mri2head - convert a CTF voxel into a head coordinate
%
% HeadVertices = ctf_mri2head(MRIVertices,mri)
%
% This function converts CTF MRI voxel coordinates into the CTF head
% coordinates.  The input 'MRIVertices' are MRI voxel locations; they are
% Nx3 voxel index (or slice) coordinates in the CTF MRI volume. The 'mri'
% input is a struct returned by ctf_read_mri, which contains the fiducials
% and a transformation matrix from MRI voxel coordinates into the head
% space coordinates.  CTF MRI volumes are 256^3 voxels, 1mm^3 each. The
% volume index has an origin at the left, anterior, superior voxel, such
% that:
%
% Sag increases from left to right (+X Right)
% Cor increases from anterior to posterior (+Y Posterior)
% Axi increases from superior to inferior (+Z Inferior)
%
% The head space coordinates are defined in relation to the MRI fiducial
% locations (nasion, left preauricular and right preauricular).  The origin
% lies half way between the left and right preauricular points and the
% coordinate axes are given as:
%
% +X is through the nasion, 
% +Y is left 
% +Z is superior
%
% CTF head coordinate axes are othogonalized by first taking the cross
% product of the +X and +Y vectors (from the origin through the nasion and
% left preauricular, respectively) to find the +Z vector.  Then, the Y axis
% is orthogonalized to the X/Z plane using the cross product and right hand
% rule.  Hence, +Y can be slightly offset from the left preauricular.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ver = '$Revision: 1.1 $';
fprintf('CTF_MRI2HEAD [v %s]\n',ver(11:15)); tic;

%--------------------------------------------------------------------
% Extract the coordinate transform matrices from the mri struct.  These are
% designed to be left multiplied into the vertices (I don't like this
% because it means the column arrays need to be transposed below, but
% that's how it is!)
%T.ctfHEAD2MRI = mri.hdr.transformMatrixHead2MRI;
%T.ctfMRI2HEAD = mri.hdr.transformMatrixMRI2Head;

% OK, so I got fed up with trying to figure out how to use these
% bloody matrices and decided to calculate them myself, so
% that the entries would work in a way that I can understand!
[trans,rot] = calc_tranfer_matrix(mri);


%--------------------------------------------------------------------
% Pad out the MRIVertices to Nx4 matrix

Nvertices = size(MRIVertices,1);

right_column = ones( Nvertices, 1 );

MRIVertices = [ MRIVertices right_column ];

%--------------------------------------------------------------------
% Convert CTF voxel indices into CTF head space coordinates

HeadVertices = (-1 * MRIVertices * trans) * inv(rot);

fprintf('...converting from mm to cm\n');
HeadVertices = HeadVertices(:,1:3) / 10;


t = toc; fprintf('...done (%6.2f sec).\n\n',t);

return







%-------------------------------------------------------
function [trans,rot] = calc_tranfer_matrix(mri)

% these are the translations!
%ctf.mri.hdr.headOrigin_sagittal: 130
%ctf.mri.hdr.headOrigin_coronal: 123
%ctf.mri.hdr.headOrigin_axial: 152

% This is how they are calculated from the fiducials:

nas(1) = mri.hdr.HeadModel_Info.Nasion_Sag;
nas(2) = mri.hdr.HeadModel_Info.Nasion_Cor;
nas(3) = mri.hdr.HeadModel_Info.Nasion_Axi;
lpa(1) = mri.hdr.HeadModel_Info.LeftEar_Sag;
lpa(2) = mri.hdr.HeadModel_Info.LeftEar_Cor;
lpa(3) = mri.hdr.HeadModel_Info.LeftEar_Axi;
rpa(1) = mri.hdr.HeadModel_Info.RightEar_Sag;
rpa(2) = mri.hdr.HeadModel_Info.RightEar_Cor;
rpa(3) = mri.hdr.HeadModel_Info.RightEar_Axi;

headOffset_Sag = ( rpa(1) - lpa(1) ) / 2;
headOrigin_Sag = lpa(1) + headOffset_Sag;

headOffset_Cor = ( rpa(2) - lpa(2) ) / 2;
headOrigin_Cor = lpa(2) + headOffset_Cor;

headOffset_Axi = ( rpa(3) - lpa(3) ) / 2;
headOrigin_Axi = lpa(3) + headOffset_Axi;

headOrigin = [headOrigin_Sag headOrigin_Cor headOrigin_Axi];

%-------------from here to...

% calculate voxel space vectors, in slices
voxNASvector = nas - headOrigin;
voxLPAvector = lpa - headOrigin;
voxRPAvector = rpa - headOrigin;

% calculate head space vectors, in mm
headNASvector = -1 * [ voxNASvector(2), voxNASvector(1), voxNASvector(3) ];
headLPAvector = -1 * [ voxLPAvector(2), voxLPAvector(1), voxLPAvector(3) ];
headRPAvector = -1 * [ voxRPAvector(2), voxRPAvector(1), voxRPAvector(3) ];

%-------------here; is encapsulated in trans:

trans = eye(3);
trans(1,1) = 0;
trans(2,1) = 1;
trans(1,2) = 1;
trans(2,2) = 0;
trans(4,1:3) = -1 * [ headOrigin(2) headOrigin(1) headOrigin(3) ];

%headVector = -1 * ctfVox * trans;


% At this point, the voxel indices are now in the head space, in
% mm; however, the head space at this point is not orthogonal

%--------------now the head space is orthogonalized

headNASunit = headNASvector / norm(headNASvector);
headLPAunit = headLPAvector / norm(headLPAvector);
headRPAunit = headRPAvector / norm(headRPAvector);

headVERTEXunit = cross( headNASunit, headLPAunit );

% revise the LPA unit vector, so it is orthogonal to Nasion
headLPAunit = cross( headVERTEXunit, headNASunit );

% these dot products = 0
%dot( headNASunit, headLPAunit )
%dot( headNASunit, headVERTEXunit )
%dot( headLPAunit, headVERTEXunit )

% Note that the LPA/RPA moves!  This has the effect of
% rotating the coordinates in the XY plane.

rot = [ headNASunit; headLPAunit; headVERTEXunit ];

return
