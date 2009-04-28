function [theta,phi,r] = elec_cart2sph(x,y,z,xo,yo,zo)

% elec_cart2sph - Convert Cartesian (X,Y,Z) to spherical (theta,phi,r)
%
% Useage: [theta,phi,r] = elec_cart2sph(X,Y,Z,xo,yo,zo)
%
% Result:   (theta,phi,r) are double floating point (Nx1);
%           r is radius,
%           theta is counterclockwise rotation from +x in x-y plane (azimuth),
%           phi is elevation with respect to Cz-axis,
%           theta & phi in radians,
%
%           +ve x-axis from origin through T4 (right ear)
%           +ve y-axis from origin through Nasion (at theta = +90degrees)
%           
% Note:     Conversion from Cartesian to spherical coordinates, with
%           origin or centroid (xo,yo,zo) is:
%           theta = atan( (y-yo) ./ (x-xo) );
%           phi = atan( sqrt( (x-xo).^2 + (y-yo).^2 ) ./ (z-zo) );
%           r = sqrt( (x-xo).^2 + (y-yo).^2 + (z-zo).^2);
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:54 $

% Licence:  GNU GPL, no implied or express warranties
% History:  07/2001, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('xo', 'var') xo = 0; end
if ~exist('yo', 'var') yo = 0; end
if ~exist('zo', 'var') zo = 0; end

theta = atan( (y-yo) ./ (x-xo) );
phi = atan( sqrt( (x-xo).^2 + (y-yo).^2 ) ./ (z-zo) );
r = sqrt( (x-xo).^2 + (y-yo).^2 + (z-zo).^2 );

return
