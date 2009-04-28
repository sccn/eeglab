function H = elec_meshplot(X,Y,Z,interp,grid_res)

% elec_meshplot - Plot a mesh and points for 3D electrode positions.
%
% Usage: H = elec_meshplot(X,Y,Z [,interp] [,grid_resolution])
%
% Uses an interpolated mesh at 0.25 cm resolution, unless
% user specifies alternative grid resolution (optional).
%
% Interpolation can be 'linear' or 'cubic'.  Unless specified,
% cubic is the default. 
%
% Returns a handle to the figure created.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:55 $

% Licence: Gnu GPL
% Author: Darren.Weber_at_radiology.ucsf.edu
% Created:  18/05/00
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  H = figure('NumberTitle','off','Name','Electrode Mesh Plot','Position',[200 200 650 500],'PaperUnits','centimeters','PaperType','A4','Units','centimeters');

  colormap(gray);

  if ~exist( 'grid_res', 'var' )  grid_res = 0.25;  end

  xi =   (min(X) - 2):grid_res:(max(X) + 2);
  yi =  [(min(Y) - 2):grid_res:(max(Y) + 2)]';

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Use 'linear' for a linear grid or 'cubic' for a polynomial grid

  if ~exist( 'interp', 'var' )  interp = 'cubic';  end
  
  [Xi,Yi,Zi] = griddata(X,Y,Z,xi,yi,interp);

  surf(Xi,Yi,Zi);

  brighten(0.75);

  mesh(Xi,Yi,Zi);

  hold on;

  plot3(X,Y,Z,'ro');
  
  rotate3d;

  clear xi yi Xi Yi Zi X Y Z;

%  view(2);  % top view
%  view(   0,0);  % back view
%  view(-180,0);  % front view
%  view(-90, 0);  % left side view
%  view( 90, 0);  % right side view
