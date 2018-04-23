% sph2spm() - compute homogenous transformation matrix from 
%             BESA spherical coordinates  to SPM 3-D coordinate
%
% Usage:
%  >> trans = sph2spm;
%
% Outputs:
%   trans - homogenous transformation matrix
%
% Note: head radius for spherical model is assumed to be 85 mm.
%
% Author: Robert Oostenveld, SMI/FCDC, Nijmegen 2005
%         Arnaud Delorme, SCCN, La Jolla 2005

% SMI, University Aalborg, Denmark http://www.smi.auc.dk/
% FC Donders Centre, University Nijmegen, the Netherlands http://www.fcdonders.kun.nl/

% Copyright (C) 2003 Robert Oostenveld, SMI/FCDC roberto@smi.auc.dk
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function besa2SPM_result = besa2SPM;

if 0
    % original transformation: problem occipital part of the haed did not
    % fit
    
    % NAS, Left EAR, Right EAR coordinates in BESA
    besa_NAS	= [0.0000	0.0913	-0.0407];
    besa_LPA	= [-0.0865	0.0000	-0.0500];
    besa_RPA	= [0.0865	0.0000	-0.0500];
    
    % NAS, Left EAR, Right EAR coordinates in SPM average
    SPM_NAS = [0 84 -48];
    SPM_LPA = [-82 -32 -54];
    SPM_RPA = [82  -32 -54];
    
    % transformation to CTF coordinate system
    % ---------------------------------------
    SPM2common  = headcoordinates(SPM_NAS , SPM_LPA , SPM_RPA,  0);
    besa2common = headcoordinates(besa_NAS, besa_LPA, besa_RPA, 0);
    
    nazcommon1 = besa2common * [ besa_NAS 1]';
    nazcommon2 = SPM2common  * [ SPM_NAS 1]';
    ratiox     = nazcommon1(1)/nazcommon2(1);
    
    lpacommon1 = besa2common * [ besa_LPA 1]';
    lpacommon2 = SPM2common  * [ SPM_LPA 1]';
    ratioy     = lpacommon1(2)/lpacommon2(2);
    
    scaling = eye(4);
    scaling(1,1) = 1/ratiox;
    scaling(2,2) = 1/ratioy;
    scaling(3,3) = mean([ 1/ratioy 1/ratiox]);
    
    besa2SPM_result = inv(SPM2common) * scaling * besa2common;
end;

if 0
    % using electrodenormalize to fit standard BESA electrode (haed radius
    % has to be 85) to BEM electrodes
    % problem: fit not optimal for temporal electrodes

    % traditional takes as input the .m field returned in the output from
    % electrodenormalize
    besa2SPM_result = traditionaldipfit([0.5588  -14.5541    1.8045    0.0004    0.0000   -1.5623    1.1889    1.0736  132.6198])
end;
   
% adapted manualy from above for temporal electrodes (see factor 0.94
% instead of 1.1889 and x shift of -18.0041 instead of -14.5541)
%traditionaldipfit([0.5588  -18.0041    1.8045    0.0004    0.0000   -1.5623    1.1889    0.94  132.6198])

besa2SPM_result = [
   0.0101   -0.9400         0    0.5588
    1.1889    0.0080    0.0530  -18.0041
   -0.0005   -0.0000    1.1268    1.8045
         0         0         0    1.0000
          ];
