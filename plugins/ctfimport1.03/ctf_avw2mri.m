function mri = ctf_avw2mri(avw)

% mri = ctf_avw2mri(avw)
%
% The purpose of this function is to convert an Analyze volume into a CTF
% .mri volume.  It currently requires that the Analyze volume is
% 256x256x256, 1x1x1 mm voxels.  The avw_read function can
% handle various Analyze orientations, so this function assumes that the
% avw.img volume has axial unflipped orientation (it always is when
% returned by avw_read, regardless of the format in the .img file). The
% returned mri struct in the matlab workspace can be saved using
% ctf_write_mri.
%
% This function depends on the avw* functions availabe in the mri_toolbox
% at http://eeg.sf.net/.  In particular, the avw_read and associated
% functions for reading an Analyze volume into a matlab avw struct, used as
% input to this function.
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

% $Revision: 1.1 $ $Date: 2009-01-30 03:49:26 $

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


mri = ctf_mri_make;

mri.file = [avw.fileprefix,'.mri'];


ver = '[$Revision: 1.1 $]';
fprintf('\nCTF_AVW2MRI [v%s]\n',ver(12:16));  tic;


% these checks for the volume dims/pixel size could be replaced with an
% interpolation function that could take any Analyze volume and convert it
% to 256x256x256, 1x1x1 mm
if avw.hdr.dime.dim(2) ~= 256,
    error('avw.hdr.dime.dim(2) ~= 256');
end
if avw.hdr.dime.dim(3) ~= 256,
    error('avw.hdr.dime.dim(3) ~= 256');
end
if avw.hdr.dime.dim(4) ~= 256,
    error('avw.hdr.dime.dim(4) ~= 256');
end

if avw.hdr.dime.pixdim(2) ~= 1,
    error('avw.hdr.dime.dim(2) ~= 256');
end
if avw.hdr.dime.pixdim(3) ~= 1,
    error('avw.hdr.dime.dim(3) ~= 256');
end
if avw.hdr.dime.pixdim(4) ~= 1,
    error('avw.hdr.dime.dim(4) ~= 256');
end


% mri.hdr.dataSize = 1 or 2 (bytes), 8 or 16 bits
if avw.hdr.dime.bitpix == 8,
    mri.hdr.dataSize = 1;
    mri.hdr.clippingRange = 255;
else
    mri.hdr.dataSize = 2;
    mri.hdr.clippingRange = 65536;
end

% This next step should always work correctly, given that avw_read always
% converts any Analyze orientation into axial unflipped in the matlab
% workspace variable avw.img; the axial unflipped orientation is 
% (+X left,  +Y anterior,  +Z superior); the ctf orientation is
% (+X right, +Y posterior, +Z inferior), the complete opposite.
temp = flipdim(avw.img,1);
temp = flipdim(temp,2);
temp = flipdim(temp,3);
mri.img = temp;

t=toc; fprintf('...done (%5.2f sec).\n\n',t);

return
