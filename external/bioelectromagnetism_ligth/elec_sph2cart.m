function [X,Y,Z] = elec_sph2cart(theta,phi,r,degrees)

% elec_sph2cart - convert spherical to Cartesian coordinates
%
% Convert from spherical (theta,phi,r) to Cartesian rectangular (x,y,z)
% coordinates.
%
% Useage: [X,Y,Z] = elec_sph2cart(theta,phi,r,[degrees])
%
%           theta is counterclockwise rotation from +x in x-y plane (azimuth),
%           phi is elevation with respect to z-axis (eg, Cz in 10-20 system),
%           r is radius,
%           degrees: 0=radians (default) or 1=degrees for theta & phi
% 
% Result:   (X,Y,Z) are double floating point
%           
% Note:     Conversion from spherical to Cartesian coordinates, given
%           angles in radians, is:
%
%           x = r .* sin(phi) .* cos(theta);
%           y = r .* sin(phi) .* sin(theta);
%           z = r .* cos(phi);
%
%           +ve x-axis from origin through T4 (right ear)
%           +ve y-axis from origin through Nasion (at theta = +90degrees)
%
%           Results from this routine are not the same as those from
%           the inbuilt sph2cart matlab function. Phi here is 
%           elevation with respect to z-axis.  The matlab function 
%           sph2cart uses elevation with respect to xy plane.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:55 $

% Licence:  GNU GPL, no express or implied warranties
% History:  07/2001, Darren.Weber_at_radiology.ucsf.edu
%           08/2003, Darren.Weber_at_radiology.ucsf.edu
%                    changed input option 'angle' to 'degrees' so
%                    the logic below is clearer.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('degrees','var'), degrees = 0; end

if degrees,
   % convert theta and phi from degrees to radians
   theta = theta * (pi/180);
   phi   = phi   * (pi/180);
end

X = r .* sin(phi) .* cos(theta);
Y = r .* sin(phi) .* sin(theta);
Z = r .* cos(phi);

return
