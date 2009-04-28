% elec_load_testscript - test loading and display of electrode data files

clear;

shape = '';
shade = 'interp';

elecpath = 'D:\MyDocuments\emse_data\ptsd-pet\source modelling\c07\meshes\';
elecfile = 'c07_124fit.txt';

elecpath = '';
elecfile = 'elec_simul32.txt';
elec_type = 'cart'; %'sphere1';


% load electrode co-ordinates and create 3D interpolated surface

file = strcat(elecpath,elecfile);
[elec,type,X,Y,Z,th,phi,r] = elec_load(file,'cart');

% Get electrode dataset centroid & reload
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

% grid, given cm input for electrode locations
grid = 1.0;    xi = -m:grid:m;    yi = xi';     gridSize = length(xi);

% Interpolate the electrode co-ordinates
switch shape
    case 'ellipse'
        [r,Xe,Ye,Ze,Xes,Yes,Zes] = elec_ellipse_fit1(X,Y,Z,xo,yo,0,gridSize); % Fit ellipse to X,Y,Z
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
        [Xi,Yi,Zi] = griddata(X,Y,Z,xi,yi,'linear');        % Cubic spline interpolation
end


figure; polar(th,phi,'d');

figure; plot3(X,Y,Z,'d'), hold on;
surf(Xi,Yi,Zi); axis tight; rotate3d;
