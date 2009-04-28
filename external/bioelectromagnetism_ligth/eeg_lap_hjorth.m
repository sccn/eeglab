function [Hjorth] = eeg_lap_hjorth(voltage,X,Y)

% eeg_lap_hjorth - 2D Laplacian of Potential at XY
%
% Useage: [Hjorth] = eeg_lap_hjorth(voltage [,X,Y])
%
% The Hjorth nearest neighbour, finite difference Laplacian.
% For a continuous approximation of the Laplacian, use the
% 'eeg_lap' function.
%
% This routine simply calls the del2 matlab command, which requires
% 'voltage' to be a rectangular matrix.  See 'help del2', esp:
%
% L = DEL2(U,HX,HY) when U is 2-D, uses the spacing specified by HX
%     and HY. If HX is a scalar, it gives the spacing between points in
%     the x-direction. If HX is a vector, it must be of length SIZE(U,2)
%     and specifies the x-coordinates of the points.  Similarly, if HY
%     is a scalar, it gives the spacing between points in the
%     y-direction. If HY is a vector, it must be of length SIZE(U,1) and
%     specifies the y-coordinates of the points.
%
% For example:
%
% [x,y] = meshgrid(-10:.5:10);
% z = (x.^2).*(y.^2); % simulate monotonic potential
% del2z = 4*del2(z);  % calculate Hjorth laplacian
% figure('name','potential vs laplacian','numbertitle','off','position',[500 10 512 512]);
% subplot(2,1,1); surf(x,y,z);
% title('potential'); shading interp; colorbar; rotate3d; axis tight
% subplot(2,1,2); surf(x,y,del2z), 
% title('Hjorth laplacian'); shading interp; colorbar; rotate3d; axis tight
%
% refs:  Hjorth B (1975).  An on-line transformation of EEG scalp
%          potentials into orthogonal source derivations. 
%          Electroencephalography & Clinical Neurophysiology, 39: 526-530.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:51 $

% Licence:  GNU GPL, no implied or express warranties
% History:  06/01, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('X','var') X = 1; end
if ~exist('Y','var') Y = 1; end

Hjorth = 4*del2(voltage,X,Y);

return
