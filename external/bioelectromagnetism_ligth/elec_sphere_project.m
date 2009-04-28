function [r,x,y,z] = elec_sphere_project(X,Y,Z,xo,yo,zo,varargin)

% elec_sphere_project - calculates spherical projections of electrode positions
%
% Usage: [r,x,y,z] = elec_sphere_project(X,Y,Z,xo,yo,zo,[estim],[plot])
%
% Notes:    The general formula for a sphere, with radius r is given by:
%
%           (x - xo)^2  +  (y - yo)^2  +  (z - zo)^2  =  r^2
%
%           This function takes arguments for cartesian co-ordinates
%           of X,Y,Z (assume Z > 0) and the center of the sphere (xo,yo,zo).
%           If (xo,yo,zo) are not provided, they are assumed (0,0,0).
%
%           Returned values are the fitted radius 'r' (constant)
%           and the (x,y,z) Cartesian co-ordinates on the projected sphere.
%
%           Options: estim - echo spherical radius estimates (default),
%                    plot  - plot the input/projected xyz on a sphere.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:55 $

% Licence:  GNU GPL, no express or implied warranties
% History:  02/2002, Darren.Weber_at_radiology.ucsf.edu
%                    adapted from elec_fit_sph
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialise centroid, unless input parameters defined
if ~exist('xo','var') xo = 0; else, if isempty(xo) xo = 0; end, end
if ~exist('yo','var') yo = 0; else, if isempty(yo) yo = 0; end, end
if ~exist('zo','var') zo = 0; else, if isempty(zo) zo = 0; end, end

if (nargin>6),        estim = varargin{1};
else                  estim = 1;
end
if (nargin>7),        plot = varargin{2};
else                  plot = 0;
end

eegversion = '$Revision: 1.1 $';
fprintf('ELEC_SPHERE_PROJECT [v %s]\n',eegversion(11:15));
fprintf('...projecting electrodes to spherical points\n');
tic;

% Initialise r0 as a rough guess at the sphere radius
rX = (max(X) - min(X)) / 2;
rY = (max(Y) - min(Y)) / 2;
rZ =  max(Z) - zo;
r0 = mean([ rX rY rZ ]);

if isequal(estim,1),
    fprintf('...initial spherical radius   = %f\n', r0); end

% perform least squares estimate of spherical radius (r)
options = optimset('fminsearch');
[r,fval,exitflag,output] = fminsearch('elec_sphere_fit_optim',...
                           r0, options, X, Y, Z, xo, yo, zo);
%fprintf('\n%s%f\n', 'Iterations = ', output.iterations);
%fprintf('\n%s%d\n', 'Exit = ', exitflag);

if isequal(estim,1),
    fprintf('...estimated spherical radius = %f\n', r); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the projection point of X,Y,Z to the fitted sphere radius r

% Convert Cartesian X,Y,Z to spherical (radians)
theta = atan2( (Y-yo), (X-xo) );
phi = atan2( sqrt( (X-xo).^2 + (Y-yo).^2 ), (Z-zo) );
% do not recalc: r = sqrt( (X-xo).^2 + (Y-yo).^2 + (Z-zo).^2);

%   Recalculate X,Y,Z for constant r, given theta & phi.
R = ones(size(phi)) * r;    
x = R .* sin(phi) .* cos(theta);
y = R .* sin(phi) .* sin(theta);
z = R .* cos(phi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the input & projected electrode positions on a sphere
if isequal(plot,1),
    figure('NumberTitle','off','Name','Electrode Placements');
    set(gca,'Projection','perspective');
    set(gca,'DataAspectRatio',[1 1 1]);

    hold on

    plot3(x,y,z,'b.');
    plot3(X,Y,Z,'ro');
    legend('input xyz','fitted sphere','Location','BestOutside');

    [Xs,Ys,Zs]=sphere;
    Xs = Xs * r;
    Ys = Ys * r;
    Zs = Zs * r;

    surf(Xs,Ys,Zs,...
        'EdgeColor','none',...
        'FaceColor',[0.7 0.7 0.7],...
        'FaceAlpha',0.4);
    view(2);
    rotate3d;

    axis tight; hold off;
end

t = toc; fprintf('...done (%6.2f sec).\n\n',t);

return
