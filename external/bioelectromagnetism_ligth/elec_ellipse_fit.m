function [r,x,y,z,Xe,Ye,Ze] = elec_ellipse_fit(X,Y,Z,xo,yo,zo,gridSize,plot)

% elec_ellipse_fit - Find radii of ellipse and elliptical points for XYZ
%
% Usage: [r,x,y,z,Xe,Ye,Ze] = elec_ellipse_fit(X,Y,Z,xo,yo,zo,gridSize,plot)
%
% Notes:    The general formula for a 3D ellipse, with center (xo,yo,zo) and
%           major axis length 2a in X, 2b in Y, and 2c in Z, is given by:
%
%           (x-xo/a)^2  +  (y-yo/b)^2  +  (z-zo/c)^2  =  1
%
%           This function takes arguments for cartesian co-ordinates
%           of X,Y,Z and returns the radius 'r(x,y,z)' and the 'x,y,z'
%           co-ordinates for the estimated elliptical surface points
%           corresponding to X,Y,Z (Z > 0 assumed).
%
%           Also returns Xe,Ye,Ze - the matrix values for an ellipsoid
%           of dimensions r(x,y,z).  The size of Xe,Ye,Ze depends on
%           'gridSize' (default 50).  Note that x,y,z are a subset of
%           Xe,Ye,Ze.
%
%           'plot' is a boolean option for plotting of the input/fitted xyz
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:54 $

% Licence:  GNU GPL no implied or express warranties
% History:  08/99   Chris Harvey
%                   function needs to swap X,Y input for some reason?
%           06/01   Darren.Weber_at_radiology.ucsf.edu
%                   replaced deprecated fmins with fminsearch function;
%                   rectified swapping of X,Y input & scaling of Z;
%                   included input options for xo,yo,zo;
%                   included output of fitted radius estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialise centroid, unless input parameters defined
if ~exist('xo','var') xo = 0; else, if isempty(xo) xo = 0; end, end
if ~exist('yo','var') yo = 0; else, if isempty(yo) yo = 0; end, end
if ~exist('zo','var') zo = 0; else, if isempty(zo) zo = 0; end, end

if ~exist('plot','var') plot = 0; else, if isempty(plot) plot = 0; end, end

eegversion = '$Revision: 1.1 $';
fprintf('ELEC_ELLIPSE_FIT [v %s]\n',eegversion(11:15));
tic;

% The data is assumed a rough hemisphere.
% To fit to a sphere, we duplicate first, with negative Z.

Zmin = min(Z); Z = Z + ( 0 - Zmin ); % translate Z so min(Z) = zero
% also used below to recalc Z
X2 = [X;X];
Y2 = [Y;Y];
Z2 = [Z;-Z];


% Initialise Ro(x,y,z) as a rough guess at axis radius in X,Y,Z
rX = (max(X2) - min(X2)) / 2;
rY = (max(Y2) - min(Y2)) / 2;
rZ = (max(Z2) - min(Z2)) / 2;
r0 = [ rX rY rZ ];

fprintf('...initial elliptical radii   (X,Y,Z) = %6.4f, %6.4f, %6.4f\n', r0);

if r0 == 0,
  msg = sprintf('...initial elliptical radii are zero!');
  error(msg);
end

% r(X,Y,Z) = radius in X,Y,Z of fitted ellipse
options = optimset('fminsearch');
[r,fval,exitflag,output] = fminsearch('elec_ellipse_fit_optim', r0, options, X2, Y2, Z2, xo, yo, zo);
ax = r(1);  by = r(2);  cz = r(3);

fprintf('...estimated elliptical radii (X,Y,Z) = %6.4f, %6.4f, %6.4f\n', r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Fit all X,Y,Z to ellipsoid
%
%   It'd be better to find the intersection of the line from (xo,yo,zo) 
%   to (x,y,z) intersecting the ellipsoid:
%
%   (sX,sY,sZ) where s = 1 / (X/a)^2 + (Y/b)^2 + (Z/c)^2
%
%   Instead, this routine generates points on the ellipsoid surface 
%   and then locates the closest of these to each (x,y,z).

if ~exist('gridSize','var')  gridSize = 50;   end

[Xe,Ye,Ze] = ellipsoid(xo,yo,zo,ax,by,cz,(gridSize-1));

Ze = max(Ze,0);	% Force into a hemi-ellipse.

% For each given X,Y,Z point, locate the nearest Xe,Ye,Ze.

x = zeros(size(X)); y = zeros(size(Y)); z = zeros(size(Z));

for e = 1:length(X)
  
  xe = X(e);  ye = Y(e);  ze = Z(e);
  
  % Generate matrix of linear distances between point (xe,ye,ze)
  % and all points on the ellipsoid (Xe,Ye,Ze).  Could be arc
  % length, but doesn't matter here.
  d = (Xe-xe).^2 + (Ye-ye).^2 + (Ze-ze).^2; d = sqrt(d);
  m = min(min(d));
  
  index = find (d <= m); i = index(1);
  
  if isfinite(Xe(i)),
    x(e) = Xe(i);
  else
    fprintf('ELEC_ELLIPSE_FIT: Warning: Xe(i) is NaN: %d\n', Xe(i));
  end
  if isfinite(Ye(i)),
    y(e) = Ye(i);
  else
    fprintf('ELEC_ELLIPSE_FIT: Warning: Ye(i) is NaN: %d\n', Ye(i));
  end
  if isfinite(Ze(i)),
    z(e) = Ze(i);
  else
    fprintf('ELEC_ELLIPSE_FIT: Warning: Ze(i) is NaN: %d\n', Ze(i)); 
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the input & projected electrode positions on a sphere
if isequal(plot,1),
    figure('NumberTitle','off','Name','Electrode Placements');
    set(gca,'Projection','perspective');
    set(gca,'DataAspectRatio',[1 1 1]);

    hold on

    plot3(X,Y,Z,'b.');
    plot3(x,y,z,'ro');
    legend('input xyz','fitted ellipse','Location','BestOutside');

    surf(Xe,Ye,Ze,...
        'EdgeColor','none',...
        'FaceColor',[0.7 0.7 0.7],...
        'FaceAlpha',0.4);
    view(2);
    rotate3d on;

    axis tight; 
    hold off;
end

z = z + Zmin;                       % restore z to original height
Ze = Ze + Zmin;

t = toc; fprintf('...done (%6.2f sec).\n\n',t);

return
