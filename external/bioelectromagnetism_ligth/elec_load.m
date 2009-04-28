function [elec,type,X,Y,Z,theta,phi,r] = elec_load(filename,coordinates,xo,yo,zo,N)

% elec_load - Read ascii electrode file
% 
% Load an ascii electrode file and return 
% the electrode labels and coordinates.  This
% function can read and return Cartesian (x,y,z) 
% and/or spherical (theta,phi,r) coordinates.
%
% The ascii file format is that of NeuroScan 
% 3Dspace export files.  Each row of the file 
% comprises an electrode label, an electrode 
% type code (see below), and the x,y,z 
% coordinates (cm).  Each field is separated 
% by spaces.  For example,
%
% Fz    69  Xcm Ycm Zcm
% 
% Usage:
% 
% [elec,type,X,Y,Z,theta,phi,r] = elec_load(file,[coordinates],[xo],[yo],[zo],[N])
%
% where:
% 
% file = 'path\filename' with row format:
%
% '%s %d %f %f %f' = [elec_label, type, (X,Y,Z) or (theta,phi,r)]
%
% N is how many electrodes to load (129 default).
% (xo,yo,zo) are the origin {(0,0,0) default}.
%
% coordinates is a string option for input data type:
% 'Cartesian'  = cartesian in Neuroscan 3Dspace format (default);
%                format is ['label' type X Y Z].
%
% 							 As of 08/2003, it is better to use 
% 							 elec_load_scan_3ddasc for Neuroscan 3Dspace
% 							 ascii export files (*.dat).  Otherwise, try
% 							 this function.
%
% 'Spherical1' = spherical, theta and phi in degrees,
% 'Spherical2' = spherical, theta and phi in radians,
%                format is ['label' type theta phi r]
%           
%  theta is counterclockwise rotation from +x in x-y plane (azimuth),
%  phi is elevation with respect to z-axis, r is radius in cm,
%  +ve x-axis from origin through T4 (right ear)
%  +ve y-axis from origin through Nasion (theta = 90degrees = pi/2 rads)
% 
% Result:
%
% elec is an electrode label, as cellstr array (Nx1)
% type is an electrode type, as uint8 (unsigned integer, Nx1),
% (X,Y,Z) & (theta,phi,r) as double floating point (Nx1), with
% theta and phi in radians.  The origin is (0,0,0), unless forced
% otherwise with the input arguments.
%           
% Notes:
% 
% i) Type is defined as follows (from Neuroscan 3Dspace ascii export):
%
%           Electrode               Type
%           ---------               ----
%           Nasion                   110 (or 78)
%           Left                     108 (or 76)
%           Right                    114 (or 82)
%           electrodes                69
%           Centroid                  99 (or 67; eg <0,0,0>)
%           Ref                      120 (or 88)
%           Scalp Points              32
%
% Open the 3Dspace ascii export (*.dat file) and check these type
% values.  If they are not consistent with the above, either hack
% the code of this function or edit the type values in the ascii
% file.
% 
% From the return values, the electrodes can be selected by:
% index = find(type == 69); % or whatever the type is in your file
% X = X(index); Y = Y(index); Z = Z(index); elec = elec(index);
%
% ii) Conversion from spherical to Cartesian coordinates is:
% 
% x = r .* sin(phi) .* cos(theta);
% y = r .* sin(phi) .* sin(theta);
% z = r .* cos(phi);
% 
% Phi here is elevation with respect to z-axis.  The matlab
% function sph2cart uses elevation with respect to xy plane.
%
% iii) Conversion from Cartesian to spherical coordinates is:
%
% theta = atan2( y, x );
% phi = atan2( sqrt( x.^2 + y.^2 ), z );
% r = sqrt( x.^2 + y.^2 + z.^2);
% 
% Phi here is the elevation with respect to z-axis.  The matlab
% function cart2sph uses elevation with respect to xy plane.
%
% Bugs:
% For NeuroScan 3Dspace export files, reading beyond the electrode
% points, into the scalp points (type 32), causes all fields to
% shift one to the left.  To avoid this error, include a short 
% string before each '32' in the electrode file or use the N
% option to only read the electrode points.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:55 $

% Licence:  GNU GPL, no express or implied warranties
% History:  08/1999, Darren.Weber_at_radiology.ucsf.edu
%           08/2003, Darren.Weber_at_radiology.ucsf.edu
%                    new version of 3Dspace exports different type numbers ;-)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('N', 'var'),  N = 129;
elseif isempty(N),      N = 129;
end
if ~exist('coordinates', 'var'),
  coordinates = 'Cartesian';
elseif isempty(coordinates),
  coordinates = 'Cartesian';
end
if ~exist('xo', 'var'), xo = 0;
elseif isempty(xo),     xo = 0;
end
if ~exist('yo', 'var'), yo = 0;
elseif isempty(xo),     yo = 0;
end
if ~exist('zo', 'var'), zo = 0;
elseif isempty(xo),     zo = 0;
end

tic;

[path,name,ext] = fileparts(filename);
file = fullfile(path,[name ext]);

eegversion = '$Revision: 1.1 $';
fprintf('ELEC_LOAD [v %s]\n',eegversion(11:15));
fprintf('...loading ''%s'' electrodes from:\n\t%s\n', coordinates, file);

switch coordinates
  
  case 'Cartesian'
    [elec,type,X,Y,Z] = textread(file,'%s %d %f %f %f', N);
    type = uint8(type);
    
    fprintf('...converting from cm to meters.\n');
    % Convert from Neuroscan 3Dspace coordinate
    % metric (cm) to meters.
    X = X ./ 100;        Y = Y ./ 100;        Z = Z ./ 100;
    
    % If possible, adjust all coordinates so that origin is (0,0,0)
    index = find(type == 99);
    if ~isempty(index),
      xo = X(index);    yo = Y(index);    zo = Z(index);
    end
    fprintf('...centering origin at (0,0,0).\n');
    X = X - xo;        Y = Y - yo;       Z = Z - zo;
    xo = X(index);    yo = Y(index);    zo = Z(index);
    
    % Convert to spherical
    fprintf('...calculating spherical coordinates.\n');
    theta = atan2( (Y-yo), (X-xo) );
    phi = atan2( sqrt( (X-xo).^2 + (Y-yo).^2 ), (Z-zo) );
    r = sqrt( (X-xo).^2 + (Y-yo).^2 + (Z-zo).^2);
    
  case 'Spherical1' %degrees
    [elec,type,theta,phi,r] = textread(file,'%s %d %f %f %f', N);
    % convert theta and phi to radians
    theta = theta * (pi/180); phi = phi * (pi/180);
    [X, Y, Z] = elec_sph2cart(theta,phi,r,1);
    fprintf('...converting from cm to meters.\n');
    X = X ./ 100;        Y = Y ./ 100;        Z = Z ./ 100;
    
  case 'Spherical2' %radians
    [elec,type,theta,phi,r] = textread(file,'%s %d %f %f %f', N);
    [X, Y, Z] = elec_sph2cart(theta,phi,r,0);
    fprintf('...converting from cm to meters.\n');
    X = X ./ 100;        Y = Y ./ 100;        Z = Z ./ 100;
    
  otherwise
    doc elec_load;
    msg = sprintf('...invalid coordinate type: ''%s'', see ''help elec_load''\n', coordinates);
    error(msg);
    
end

fprintf('...loaded %d electrodes\n', size(elec,1));

t = toc; fprintf('...done (%6.2f sec).\n\n',t);

return
