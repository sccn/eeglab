function avw = ctf_mri2avw(mri)

% avw = ctf_mri2avw(mri)
%
% The purpose of this function is to convert a CTF .mri volume into 
% an Analyze volume.  The returned avw struct can be saved using avw_write.
%
% This function depends on the avw* functions availabe in the mri_toolbox
% at http://eeg.sf.net/.
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

% History:  04/2004, Darren.Weber_at_radiology.ucsf.edu
%                    - adapted from an appendex to CTF document
%                    MRIConverter.pdf, which is copied at the end of this
%                    function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

avw = avw_hdr_make;

avw.fileprefix = strrep(mri.file,'.mri','');

ver = '[$Revision: 1.1 $]';
fprintf('\nCTF_MRI2AVW [v%s]\n',ver(12:16));  tic;

% CTF .mri files are always 256x256x256 voxels, 1mm^3
avw.hdr.dime.dim(2) = 256;
avw.hdr.dime.dim(3) = 256;
avw.hdr.dime.dim(4) = 256;
avw.hdr.dime.pixdim(2) = 1;
avw.hdr.dime.pixdim(3) = 1;
avw.hdr.dime.pixdim(4) = 1;

% mri.hdr.dataSize = 1 or 2 (bytes), 8 or 16 bits
if mri.hdr.dataSize == 1,
    avw.hdr.dime.bitpix = 8;
    avw.hdr.dime.datatype = 2;
else
    avw.hdr.dime.bitpix = 16;
    avw.hdr.dime.datatype = 4;
end

% This next step should always work correctly, given that avw structs
% are always axial unflipped in the matlab workspace; the axial 
% unflipped orientation is 
% (+X left,  +Y anterior,  +Z superior); the ctf orientation is
% (+X right, +Y posterior, +Z inferior), the complete opposite.
temp = flipdim(mri.img,1);
temp = flipdim(temp,2);
temp = flipdim(temp,3);
avw.img = temp;

t=toc; fprintf('...done (%5.2f sec).\n\n',t);

return
