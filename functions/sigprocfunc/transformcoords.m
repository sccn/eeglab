% transformcoords() - Select nazion and inion in anatomical MRI images.
%
% Usage:
%   mewcoords = transformcoords(coords, rotate, scale, center, reverse);
%
% Inputs:
%   coords    - array of 3-D coordinates (3 by N or N by 3)
%   rotate    - [pitch roll yaw] rotate in 3-D using pitch (x plane), 
%               roll (y plane) and yaw (z plane). An empty array does
%               not perform any rotation.
%   scale     - [scalex scaley scalez] scale axis. A single numeric
%               input scale all the dimentions the same. Default 1
%               does not scale.
%   shifts    - [x y z] shift coordinates (after rotation and scaling). 
%               Default [0 0 0] does not move the center.
%   reverse   - [0|1] when set to 1 perform the reverse transformation,
%               first moving to the old center, unscaling, and unrotating.
%               Default is 0.
%
% Output:
%   newcoords - coordinates after rotating, scaling and recentering
%
% Author: Arnaud Delorme, Salk, SCCN, UCSD, CA, March 23, 2004

% Copyright (C) 2004 Arnaud Delorme
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

function coords = transformcoords(coords, rotate, scale, center, reverse);
    
    if nargin < 2
        help transformcoords;
        return;
    end
    if nargin < 3
        scale =1;
    end
    if nargin < 4
        center = [0 0 0];
    end
    if nargin < 5
        reverse = 0;
    end
    if size(coords, 1) ~= 3
        coords = coords';
        trp    = 1;
    else
        trp = 0;
    end
    if size(coords, 1) ~= 3
        error('Number of columns must be 3 for the coordinate input');
    end
    if length(rotate) > 0 && length(rotate) ~= 3
        error('rotate parameter must have 3 values');
    end
    
    % decode parameters
    % -----------------
    centx = -center(1);
    centy = -center(2);
    centz = -center(3);
    if length(scale) == 1
        scale = [scale scale scale];
    end
    scalex = scale(1);
    scaley = scale(2);
    scalez = scale(3);
    if length(rotate) < 3
        rotate = [0 0 0]
    end
    pitch = rotate(1);
    roll  = rotate(2);
    yaw   = rotate(3);
    
    
    if ~reverse
        % pitch roll yaw rotation
        % -----------------------
        % pitch (x-axis); roll = y axis rotation; yaw = z axis
        % see http://bishopw.loni.ucla.edu/AIR5/homogenous.html
        cp = cos(pitch); sp = sin(pitch);
        cr = cos(roll);  sr = sin(roll);
        cy = cos(yaw);   sy = sin(yaw);
        rot3d = [ cy*cr+sy*sp*sr    sy*cr-cy*sp*sr     cp*sr  ;
                  -sy*cp            cy*cp               sp     ;
                  sy*sp*cr-cy*sr    -cy*sp*cr-sy*sr     cp*cr  ];
        coords = rot3d*coords;
        
        % scaling and centering
        % ---------------------
        coords(1,:) = coords(1,:)*scalex-centx;
        coords(2,:) = coords(2,:)*scaley-centy;
        coords(3,:) = coords(3,:)*scalez-centz;
    else
        % unscaling and uncentering
        % -------------------------
        coords(1,:) = (coords(1,:)+centx)/scalex;
        coords(2,:) = (coords(2,:)+centy)/scaley;
        coords(3,:) = (coords(3,:)+centz)/scalez;
        
        % pitch roll yaw rotation
        % -----------------------
        cp = cos(-pitch); sp = sin(-pitch);
        cr = cos(-roll);  sr = sin(-roll);
        cy = cos(-yaw);   sy = sin(-yaw);
        rot3d = [ cy*cr+sy*sp*sr    sy*cr-cy*sp*sr     cp*sr  ;
                  -sy*cp            cy*cp               sp     ;
                  sy*sp*cr-cy*sr    -cy*sp*cr-sy*sr     cp*cr  ];
        coords = rot3d*coords;
    end
    
    if trp
        coords = coords';
    end
