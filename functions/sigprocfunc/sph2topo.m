% sph2topo() - Convert from a 3-column headplot file in spherical coordinates
%              to 3-column topoplot() locs file in polar (not cylindrical) coords.
%              Used for topoplot() and other 2-D topographic plotting programs.
%              Assumes a spherical coordinate system in which horizontal angles 
%              have a range [-180,180] deg,  with zero pointing to the right ear. 
%              In the output polar coordinate system, zero points to the nose.
%              See  >> help readlocs
% Usage:
%          >> [chan_num,angle,radius] = sph2topo(input,shrink_factor,method);
%
% Inputs:
%   input         = [channo,az,horiz] = chan_number, azumith (deg), horiz. angle (deg)
%                   When az>0, horiz=0 -> right ear, 90 -> nose 
%                   When az<0, horiz=0 -> left ear, -90 -> nose
%   shrink_factor = arc_length shrinking factor>=1 (deprecated).
%                   1 -> plot edge is 90 deg azimuth {default};
%                   1.5 -> plot edge is +/-135 deg azimuth See 
%                   >> help topoplot(). 
%   method        = [1|2], optional. 1 is for Besa compatibility, 2 is for
%                   compatibility with Matlab function cart2sph(). Default is 2
%
% Outputs:
%   channo  = channel number (as in input)
%   angle   = horizontal angle (0 -> nose; 90 -> right ear; -90 -> left ear)
%   radius  = arc_lengrh from vertex (Note: 90 deg az -> 0.5/shrink_factor);
%             By topoplot() convention, radius=0.5 is the nasion-ear_canal plane.
%             Use topoplot() 'plotrad' to plot chans with abs(az) > 90 deg.
%
% Author: Scott Makeig & Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 6/12/98 
%
% See also: cart2topo(), topo2sph()

% Copyright (C) 6/12/98 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

% corrected left/right orientation mismatch, Blair Hicks 6/20/98
% changed name sph2pol() -> sph2topo() for compatibility -sm
% 01-25-02 reformated help & license -ad 
% 01-25-02 changed computation so that it works with sph2topo -ad 

function [channo,angle,radius] = sph2topo(input,factor, method)

chans = size(input,1);
angle = zeros(chans,1);
radius = zeros(chans,1);

if nargin < 1
   help sph2topo
   return
end
   
if nargin< 2
  factor = 0;
end
if factor==0
  factor = 1;
end
if factor < 1
  help sph2topo
  return
end

if size(input,2) ~= 3
   help sph2topo
   return
end

channo = input(:,1);
az = input(:,2);
horiz = input(:,3);

if exist('method')== 1 && method == 1
  radius = abs(az/180)/factor;
  i = find(az>=0);
  angle(i) = 90-horiz(i);
  i = find(az<0);
  angle(i) = -90-horiz(i);
else
  angle  = -horiz;
  radius = 0.5 - az/180;
end
