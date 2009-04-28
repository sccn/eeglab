function [f] = elec_fit_ellipse_optim(r, X, Y, Z, xo, yo, zo)

% elec_fit_ellipse_optim - Optimization for elec_fit_ellipse.m
%
% Called from elec_fit_ellipse.m
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:55 $

% Licence:  GNU GPL, no implied or express warranties
% History:  02/2002, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% r is a 3x1 vector of radius values for each x,y,z axis component of ellipse
%
% equation of ellipsoid with center (xo,yo,zo) and radius for each axis (x,y,z) = (a,b,c):
% (( x - xo )^2 / a^2) + (( y - yo )^2 / b^2) + (( z - zo )^2 / c^2) = 1
%
% This function below creates a scalar value to
% return to the fminsearch function in elec_fit_sphere.

E = ( (X-xo).^2 )/r(1).^2 + ((Y-yo).^2)/r(2).^2 + ((Z-zo).^2)/r(3).^2  - 1;

f = sum( E .* E );  % sum of squares returned

return
