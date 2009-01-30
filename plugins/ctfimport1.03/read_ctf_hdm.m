function [vol] = read_ctf_hdm(filename);

% READ_CTF_HDM reads the head volume conductor model from a *.hdm file
%
% vol = read_ctf_hdm(filename)

% Copyright (C) 2003, Robert Oostenveld
% 
% $Log: not supported by cvs2svn $
% Revision 1.1  2005/12/06 06:24:23  psdlw
% Alternative functions from the FieldTrip package, which is now released under GPL (so I assume these functions can be committed to the sourceforge cvs)
%
% Revision 1.1  2003/03/24 12:30:42  roberto
% new implementation
%

ascii = read_ctf_ascii(filename);

if isfield(ascii, 'MultiSphere_Data')
  chans = fieldnames(ascii.MultiSphere_Data);
  % remove the first two fields: SEARCH_RADIUS and HEADSHAPE_FILE
  chans = chans(3:end);
  for i=1:length(chans)
    tmp = getfield(ascii.MultiSphere_Data, chans{i});
    vol.label{i} = chans{i};
    vol.r(i) = tmp(4);
    vol.o(i, :) = tmp(1:3);
  end
  vol.r = vol.r(:);	% ensure column vector
elseif isfield(ascii, 'MEG_Sphere')
  vol.r    = ascii.MEG_Sphere.RADIUS;
  vol.o(1) = ascii.MEG_Sphere.ORIGIN_X;
  vol.o(2) = ascii.MEG_Sphere.ORIGIN_Y;
  vol.o(3) = ascii.MEG_Sphere.ORIGIN_Z;
else
  error('no headmodel information found');
end

% add the fiducials, these are in raw MRI coordinates
if isfield(ascii, 'Fid_Points')
  vol.mri.nas = ascii.Fid_Points.NASION;
  vol.mri.lpa = ascii.Fid_Points.LEFT_EAR;
  vol.mri.rpa = ascii.Fid_Points.RIGHT_EAR;
end

% add the units in which the volume conductor is defined
vol.unit = 'cm';

