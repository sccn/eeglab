function [Xe,Ye,Ze] = elec_2d_3d(X,Y,Xrad,Yrad,Zrad)

% elec_2d_3d - Project 2D planar Cartesian coordinates into 3D elliptical shape
%
% Useage: [Xe,Ye,Ze] = elec_2d_3d(X,Y,Xrad,Yrad,Zrad)
%
%   Project 2D planar Cartesian coordinates into 3D elliptical shape,
%   using gnomonic projection, which draws a line from points on a plane 
%   tangential to north pole, to a point halfway between equator and south pole, 
%   where the intersection with the ellipse is the position of the electrode.
%
%   Given:      Set of 2D Cartesian coordinates on a plane (X,Y)
%	            3 radii of ellipsoid (Xrad,Yrad,Zrad)
%
%   See also elec_3d_2d.m
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:54 $

% Licence:  GNU GPL, no implied or express warranties
% History:  10/1999, Chris Harvey
%           07/2001, Darren.Weber_at_radiology.ucsf.edu
%                    - using matrix algebra rather than indexed looping
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% The size of X and Y should be the same.
% The size of Xe,Ye,Ze will be set to the same.
%
Xe = X;
Ye = Y;
Ze = Y;

B = -1/2;
C = -3/4;

for i = 1:length(X)

    A = (X(i)/Xrad)^2 + (Y(i)/Yrad)^2 + 1/4;

    t1 = (-B + sqrt(B^2 - 4*A*C)) / (2*A);

    t2 = (-B - sqrt(B^2 - 4*A*C)) / (2*A);

    t = max(t1,t2);

    Xe(i) = X(i)*t;
	Ye(i) = Y(i)*t;
	Ze(i) = (-Zrad/2)*(1-t);
end

