function [x,y,z] = ellipse(a,b,c,N)

% elec_ellipse_points - Generate points on an ellipse
%
% Useage: [X,Y,Z] = elec_ellipse_points(a,b,c,N)
%
%           a   X radius
%           b   Y radius
%           c   Z radius
%           N   Number of points to generate
%
%           [X,Y,Z]     cartesian points
%
% Example:  [x,y,z] = elec_ellipse_points; surf(x,y,z); shading interp
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:54 $

% Licence:  GNU GPL, no implied or express warranties
% History:  08/01 Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('a','var'), a = 20; end
if ~exist('b','var'), b = 15; end
if ~exist('c','var'), c = 15; end
if ~exist('N','var'), N = 124; end

theta = linspace(0,2*pi,N);
phi   = linspace(-pi,pi,N)';

sinphi = sin(phi);
cosphi = cos(phi);
sintheta = sin(theta);
costheta = cos(theta);
invrho = (sinphi*sintheta/a).^2 + (sinphi*costheta/b).^2 + (cosphi*ones(size(theta))/c).^2;
invrho = sqrt(invrho);

x = sinphi*costheta./invrho;
y = sinphi*sintheta./invrho;
z = cosphi*ones(size(theta))./invrho;

return
