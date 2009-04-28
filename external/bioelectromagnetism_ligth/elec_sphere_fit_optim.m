function [f] = elec_sphere_fit_optim(r, X, Y, Z, xo, yo, zo)

% elec_sphere_fit_optim - Optimization for elec_sphere_fit.m
%
% Called from elec_sphere_fit.m
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:55 $

% Licence:  GNU GPL, no implied or express warranties
% History:  02/2002, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% with center (Xo,Yo,Zo) and radius r, the equation of a sphere is:
%
% r^2 = (x-xo)^2  +  (y-yo)^2  +  (z-zo)^2
%
% This function below creates a scalar value to
% return to the fminsearch function in elec_sphere_fit.

S = (X-xo).^2  +  (Y-yo).^2  +  (Z-zo).^2  -  r^2;

f = sum( S.^2 );
