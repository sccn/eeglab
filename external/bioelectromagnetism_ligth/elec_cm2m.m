function [X,Y,Z] = elec_cm2m(X,Y,Z)

% elec_cm2m - Convert electrode coordinates from centimeters to meters
%
%   [X,Y,Z] = elec_cm2m(X,Y,Z)
%
%   Simply:  X = X / 100;  Y = Y / 100;  Z = Z / 100;
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:54 $

% Licence:  GNU GPL, no implied or express warranties
% History:  07/2001, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    X = X / 100;  Y = Y / 100;  Z = Z / 100;
