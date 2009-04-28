function [r,x,y,z] = elec_sphere_fit(X,Y,Z,xo,yo,zo,estim)

% elec_sphere_fit - find radius of sphere and closest spherical points to XYZ
%
% Usage: [r,x,y,z] = elec_sphere_fit(X,Y,Z,xo,yo,zo)
%
% Notes:    The general formula for a sphere, with radius r is given by:
%
%           (x - xo)^2  +  (y - yo)^2  +  (z - zo)^2  =  r^2
%
%           This function takes arguments for cartesian co-ordinates
%           of X,Y,Z (assume Z > 0) and the center of the sphere (xo,yo,zo).
%           If (xo,yo,zo) are not provided, they are assumed (0,0,0).
%
%           Returns values are the radius 'r' and the (x,y,z) Cartesian
%           co-ordinates for the fitted sphere.
%
%           See also, elec_sphere_project.m
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:55 $

% Licence:  GNU GPL, no implied or express warranties
% History:  06/01   Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialise centroid, unless input parameters defined
if ~exist('xo','var') xo = 0; else, if isempty(xo) xo = 0; end, end
if ~exist('yo','var') yo = 0; else, if isempty(yo) yo = 0; end, end
if ~exist('zo','var') zo = 0; else, if isempty(zo) zo = 0; end, end

% The data is assumed a rough hemisphere.
% To fit to a sphere, we duplicate first, with negative Z.

Zmin = min(Z); Z = Z + ( 0 - Zmin ); % translate Z so min(Z) = zero
% also used below to recalc Z
X2 = [X;X];
Y2 = [Y;Y];
Z2 = [Z;-Z];

% Initialise r0 as a rough guess at the sphere radius
rX = (max(X2) - min(X2)) / 2;
rY = (max(Y2) - min(Y2)) / 2;
rZ = (max(Z2) - min(Z2)) / 2;
r0 = mean([ rX rY rZ ]);

if ~exist('estim','var') estim = 1; 
else, if isempty(estim) estim = 1; end, end

if estim ==1, 
  fprintf('\n%s%f\n', 'Initial spherical radius   = ', r0); end

% perform least squares estimate of spherical radius (r)
options = optimset('fminsearch');
[r,fval,exitflag,output] = fminsearch('elec_sphere_fit_optim',...
  r0, options, X2, Y2, Z2, xo, yo, zo);
%fprintf('\n%s%f\n', 'Iterations = ', output.iterations);
%fprintf('\n%s%d\n', 'Exit = ', exitflag);

if estim ==1,
  fprintf('%s%f\n',   'Estimated spherical radius = ', r); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Fit all X,Y,Z to sphere
%
%   It'd be better to find the intersection of the line from (xo,yo,zo) 
%   to (x,y,z) intersecting the sphere:
%
%   (sX,sY,sZ) where s = 1 / (X/r)^2 + (Y/r)^2 + (Z/r)^2
%
%   Instead, this routine generates points on the sphere 
%   and then locates the closest of these to each (x,y,z).

[Xs,Ys,Zs] = elec_sphere_points(40,80,r);

% For each given X,Y,Z point, locate the nearest Xs,Ys,Zs.

x = zeros(size(X)); y = zeros(size(Y)); z = zeros(size(Z));

for s = 1:length(X)
  
  xs = X(s);  ys = Y(s);  zs = Z(s);  % Get original values
  
  % Generate matrix of linear distances between point (X,Y,Z)
  % and all points on the sphere (Xs,Ys,Zs).  Could be arc
  % length, but doesn't matter here.
  d = sqrt( (Xs-xs).^2 + (Ys-ys).^2 + (Zs-zs).^2 );
  m = min(min(d));
  
  index = find (d <= m); i = index(1);
  
  if isfinite(Xs(i)) x(s) = Xs(i); 
  else fprintf('\n\aWarning: Xs(i) is NaN: %d', Xs(i)); end
  if isfinite(Ys(i)) y(s) = Ys(i); 
  else fprintf('\n\aWarning: Ys(i) is NaN: %d', Ys(i)); end
  if isfinite(Zs(i)) z(s) = Zs(i); 
  else fprintf('\n\aWarning: Zs(i) is NaN: %d', Zs(i)); end
  
end
z = z + Zmin;    % restore z to original height
Zs = Zs + Zmin;
