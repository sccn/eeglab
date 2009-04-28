% eeg_lap_test_script - explore the vector calculus of scalp current density.
% 
%   It starts with a potential field (ie, scalp voltage),
%   which is a scalar field in x,y.
%
%   We then calculate the first directional derivative, called
%   the gradient.  This is a vector field with unit magnitude and
%   direction of greatest incline of the potential field (ie,
%   all vectors are directed away from a negative potential peak 
%   and toward a positive potential peak).  
%
%   The next step is to calculate the current density vector.  This
%   depends on the conductivity (sigma) of the material from which the
%   potential field is measured, given by:
%
%   current density =  -sigma * grad(v)
%
%   Lets assume the scalp potential is measured from a scalp 
%   with constant conductivity and a good estimate of this value 
%   is 0.35 Siemens/m (0.035 Siemens/cm).  Given a gradient of the
%   potential field over meters, we have the following:
%
%   scalp current density = -0.35 * grad(v)
%
%   The next step is to calculate the divergence of the gradient 
%   vector.  This is a scalar field that indicates the amount of 
%   expansion (+ve) or contraction (-ve) in the direction of the 
%   gradient field.
%
%   div( scalp current density ) = -0.35 * div( grad(v) )
%
%   Mathematically, this step is performed with the Laplacian
%   operator (del^2).  The Laplacian is equivalent to the divergence
%   of the gradient - expressed as:
%
%   div( grad(v) )     or     del . del(v)     or     del^2(v)
%
%   However, for a Guassian scalp surface, with no current sources/sinks,
%   the divergence of the gradient is zero, because no current can
%   be conducted along a scalp normal direction from the scalp into the 
%   air.
%
%   To solve this problem we need the sum of the partial derivates of
%   scalp potential in x and y (tangential components).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all;

SHADE = 'interp'; % FLAT, INTERP, {FACETED}
VIEW = [-20,20];

rot = 0; % Zaxis label rotation
align = 'right';

fig1 = figure('Position',[10 10 1200 900]), rotate3d;
sub11 = subplot(3,2,1); sub12 = subplot(3,2,2); sub13 = subplot(3,2,3); 
sub14 = subplot(3,2,4); sub15 = subplot(3,2,5); sub16 = subplot(3,2,6);

fig2 = figure('Position',[500 100 500 800]);
sub21 = subplot(2,1,1); sub22 = subplot(2,1,2);

%   A scalar potential field in x,y is any function f(x,y), eg:

    MIN = -pi/2; MAX = pi/2; grid = pi/10;
    [x,y] = meshgrid(MIN:grid:MAX);    r = x.^2 + y.^2;
    
    %   Model a charge distributions
    Q1 = sin(x);    Q2 = sin(y);    z = (Q1 .* Q2);
    % cubic spline interpolation of charge density rho(x,y)
    zi = interp2(z,1,'cubic');    xi = linspace(MIN,MAX,size(zi,1));    yi = xi;
    
    figure(fig1); subplot(sub11); surfc(xi,yi,zi), shading (SHADE), hold on, colorbar, view(VIEW);
    title 'Charge Density (C/m^2)', zlabel('coulombs (C)','Rotation',rot,'HorizontalAlignment',align), axis tight
    
% Electric Field Strength
    
    %e0 = 8.854E-12;           % permittivity of free space
    %k = 1 / (4 * pi * e0);    % permittivity constant (8.9877e+009)
    
    [u,v,en] = surfnorm(x,y,z);      % 3D surface normals
    en = en.*z.*r;                % orient normal in direction of charge & scale by distance
    
    [ex,ey] = gradient(-z,grid,grid);   % electric field, gradient of -charge (del[-z])
    ex = ex.*r;                      % scale ex component by distance
    ey = ey.*r;                      % scale ey component by distance
    
    figure(fig1); subplot(sub12); surfc(xi,yi,zi), shading (SHADE), hold on, colorbar, view(VIEW)
    title 'Electric Field (E, newtons/coulomb * 8.9877e09)', zlabel('N/C','Rotation',rot,'HorizontalAlignment',align)
    quiver3(x,y,z,ex,ey,en,0), axis tight
    
    figure(fig2); subplot(sub21); contour(xi,yi,zi), hold on, quiver(x,y,ex,ey,0), hold off
    title 'Electric Field Vectors (N/C)'
    
%   Electric potential (volts = joules/coulomb) is proportional to minus Electric Field Strength
%   and equivalent to the charge density in free space
    
    figure(fig1); subplot(sub13); surfc(xi,yi,zi), shading (SHADE), hold on, colorbar, view(VIEW);
    title 'Electric Potential (V; 1 volt = 1 joule/coulomb)', axis tight
    zlabel('Volt','Rotation',rot,'HorizontalAlignment',align)
    
    
% Current Density = Conductivity * Electric Field Strength

    % For a material conducting medium, current density is equal to
    % electric field strength multiplied by conductivity (sigma = 1/resistivity).
    
    sigma = 0.35;    % conductivity of human scalp tissues
    
    cdx = ex.*sigma; % scale ex component by conductivity
    cdy = ey.*sigma; % scale ey component by conductivity
    cdn = en.*sigma; % scale en component by conductivity
    
    figure(fig1); subplot(sub14); surfc(xi,yi,zi), shading (SHADE), hold on, colorbar, view(VIEW)
    title 'Current Density (Jscalp = E * 0.35)', zlabel('Amp/m^2','Rotation',rot,'HorizontalAlignment',align)
    quiver3(x,y,z,cdx,cdy,cdn,0), axis tight
    
    figure(fig2); subplot(sub22); contour(xi,yi,zi), hold on, quiver(x,y,cdx,cdy,0), hold off
    title 'Current Density 2D Gradient'
    
%   The Laplacian of a scalar field, eg electric potential (voltage)
    
    lap = 4*del2(z,grid,grid);            			% Hjorth, nearest neighbour Laplacian
    lapi = interp2(lap,1,'cubic');
    
    figure(fig1); subplot(sub15); surfc(xi,yi,lapi), shading (SHADE), colorbar, axis tight, view(VIEW)
    title 'Laplacian of Electric Potential (Lap(V) = div[grad[V]])'
    zlabel('Volt/m^2','Rotation',rot,'HorizontalAlignment',align);
    
    %figure('Name','Potential, Coloured by Laplacian','Position',[650 10 600 400])
    %surfc(xi,yi,zi,lapi), shading (SHADE), colorbar, axis tight, view(VIEW)
    
%   The divergence of the current density; Div(Current Density) = - conductivity * Laplacian(voltage)
    
    scdi = -0.35 * lapi;
    
    figure(fig1); subplot(sub16); surfc(xi,yi,scdi), shading (SHADE), colorbar, axis tight, view(VIEW)
    title 'Divergence of Current Density (div(J) = -0.35 * Lap(V))', 
    zlabel('Amp/m^2','Rotation',rot,'HorizontalAlignment',align);
