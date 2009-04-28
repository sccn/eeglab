% eeg_topo_surface_script - explore surface interpolation/plots of ERP potentials

clear;

shape = 'ellipse';
shade = 'interp';

elecpath = 'D:\MyDocuments\emse_data\ptsd-pet\source modelling\c07\meshes\';
elecfile = 'c07_124fit.txt';

voltpath = 'D:\MyDocuments\emse_data\ptsd-pet\source modelling\c07\eeg\';
voltfile = 'c07oac.dat';

lappath = 'D:\MyDocuments\emse_data\ptsd-pet\lap14hz\';
lapfile  = 'c07oac_lap14hz.dat';

hjorthpath = 'D:\MyDocuments\emse_data\ptsd-pet\source modelling\c07\eeg\';
hjorthfile = 'c07oac_hjorth.dat';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load electrode co-ordinates and create 3D interpolated surface

file = strcat(elecpath,elecfile);
[elec,type,X,Y,Z,th,phi,r] = elec_load(file,'cart');
% Get electrode dataset centroid and reload
index = find(type == 99);    xo = X(index);    yo = Y(index);    zo = Z(index);
[elec,type,X,Y,Z,th,phi,r] = elec_load(file,'cart',xo,yo,zo);

% Select electrodes only
index = find(type == 69);
elec = elec(index);
X = X(index);    Y = Y(index);    Z = Z(index);
th = th(index);   phi = phi(index);    r = r(index);

% Gradations within co-ordinates
Xmin = floor(min(X));    Xmax = ceil(max(X));
Ymin = floor(min(Y));    Ymax = ceil(max(Y));

if abs(Xmin) < abs(Xmax) x = abs(Xmax); else x = abs(Xmin); end
if abs(Ymin) < abs(Ymax) y = abs(Ymax); else y = abs(Ymin); end
if (x < y) m = y; else m = x; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate the electrode co-ordinates

% grid, given cm input for electrode locations
grid = 1.0;    xi = -m:grid:m;    yi = xi';     gridSize = length(xi);

switch shape
    case 'ellipse'
        [r,Xe,Ye,Ze,Xes,Yes,Zes] = elec_ellipse_fit(X,Y,Z,xo,yo,0,gridSize); % Fit ellipse to X,Y,Z
        [Xi,Yi,Zi] = griddata(Xe,Ye,Ze,xi,yi,'cubic');    % Cubic spline interpolation
        %figure; plot3(X,Y,Z,'d'), hold on; plot3(Xe,Ye,Ze,'r.'), hold off; rotate3D
        %figure; surf(Xes,Yes,Zes); shading (shade); view([-20 20]); rotate3D
        clear Xes Yes Zes;
    case 'sphere'
        [r,Xs,Ys,Zs,Xss,Yss,Zss] = elec_sphere_fit1(X,Y,Z,xo,yo,0,gridSize); % Fit sphere to Z
        [Xi,Yi,Zi] = griddata(Xs,Ys,Zs,xi,yi,'cubic');   % Cubic spline interpolation
        %figure; plot3(X,Y,Z,'d'), hold on; plot3(Xs,Ys,Zs,'r.'), hold off; rotate3D
        %figure; surf(Xss,Yss,Zss); shading (shade); view([-20 20]); rotate3D
        clear Xss Yss Zss;
    otherwise
        [Xi,Yi,Zi] = griddata(X,Y,Z,xi,yi,'cubic');        % Cubic spline interpolation
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% load and plot data (124 rows x 680 columns)

timecol = 450;

% Original Voltage
file = strcat(voltpath,voltfile);
volt = eeg_ascii_load(file);
%volt = eeg_uv2volt(volt);
V = volt(:,timecol);
Vi = griddata(X,Y,V,xi,yi,'cubic');
%    figure; surfc(Xi,Yi,Vi); title 'Orig Potential'; shading (shade); colorbar; rotate3d
figure; surf(Xi,Yi,Zi,Vi); title 'Orig Potential'; shading (shade); colorbar; rotate3d

% May have to convert Vi to spherical coordinates for this to work.
%    figure; surf(Xes,Yes,Zes,Vi); title 'Orig Potential'; shading (shade); colorbar; rotate3d


% EMSE Laplacian of spherical spline
%    file = strcat(lappath,lapfile);
%    lap = eeg_ascii_load(file);
%    L = lap(:,timecol);                     
%    Li = griddata(X,Y,L,xi,yi,'cubic');
%    figure; surfc(Xi,Yi,Li); title 'EMSE Laplacian'; shading (shade); colorbar; rotate3d
%    figure; surf(Xi,Yi,Zi,Li); title 'EMSE Laplacian'; shading (shade); colorbar; rotate3d


% My Laplacian of spherical spline
%    Lsn = eeg_lap(V,X,Y,Z,xo,yo,zo);
%    Lsni = griddata(X,Y,Lsn,xi,yi,'cubic');
%    figure; surf(Xi,Yi,Lsni); title 'My Spline Laplacian'; shading (shade); colorbar; rotate3d
%    figure; surf(Xi,Yi,Zi,Lsni); title 'My Spline Laplacian'; shading (shade); colorbar; rotate3d


% EMSE hjorth
%    file = strcat(hjorthpath,hjorthfile);
%    hjorth = eeg_ascii_load(file);
%    hjorth = hjorth(:,1:124)';
%    H = hjorth(:,timecol);
%    Hi = griddata(X,Y,H,xi,yi,'cubic');
%    figure; surf(Xi,Yi,Zi,Hi); title 'EMSE Hjorth'; shading (shade); colorbar; rotate3d


% My Hjorth
Hji = eeg_lap_hjorth(Vi,xi,yi);
%    figure; surfc(Xi,Yi,Hji); title 'My Hjorth Laplacian'; shading (shade); colorbar; rotate3d
figure; surf(Xi,Yi,Zi,Hji); title 'My Hjorth Laplacian'; shading (shade); colorbar; rotate3d


% Scalp current density
%    SCDi = eeg_lap2scd(Li);
%    figure; surf(Xi,Yi,Zi,SCDi); title 'Spline SCD'; shading (shade); colorbar; rotate3d
