function [Xt,Yt] = elec_3d_2d(X,Y,Z,Zrad)

% elec_3d_2d - Project Cartesian 3D coordinates to a 2D plane.
%
% [Xt,Yt] = elec_3d_2d(X,Y,Z,Zrad)
%
% Project 3D electrode positions onto a 2D plane, using 
% gnomonic projection, which draws a line from a point 
% halfway between equator and south pole, through electrode, 
% to a plane tangential to north pole.
%
% Given:  Set of 3D Cartesian coordinates (X,Y,Z with Z > 0)
%         X,Y midpoint at 0
%         Zrad is radius in Z
%
% See also elec_2d_3d.m
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:54 $

% Licence:  GNU GPL, no implied or express warranties
% History:  10/1999, Chris Harvey
%           07/2001, Darren.Weber_at_radiology.ucsf.edu
%                    - using matrix algebra rather than indexed looping
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rad = Zrad / 2;

t = (Z + rad) * ( 1/rad );

Xt = X .* t;
Yt = Y .* t;

return
