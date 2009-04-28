% eeg_topo_surfacespline_script - explore spline toolbox for topographic mapping


%   If you want to interpolate at sites other than the breaks and/or by 
%   splines other than cubic splines with simple knots, then you use 
%   the spapi command. In its simplest form, you would say:
%
%   sp = spapi(k,x,y);
%
%   in which the first argument, k, specifies the order of the interpolating 
%   spline; this is the number of coefficients in each polynomial piece, 
%   i.e., 1 more than the nominal degree of its polynomial pieces. For 
%   example, the next figure shows a linear, a quadratic, and a quartic 
%   spline interpolant to our data, as obtained by the statements:
%
%   x = (4*pi)*[0 1 rand(1,20)]; y = sin(x);
%   figure
%   sp2 = spapi(2,x,y); fnplt(sp2), hold on
%   sp3 = spapi(3,x,y); fnplt(sp3,'--')
%   sp5 = spapi(5,x,y); fnplt(sp5,'-.'), plot(x,y,'o')
%   legend('linear','quadratic','quartic','data'), hold off


%   What if the data are noisy? For example, suppose that the given values are:
%
%   noisy = y + .1*(rand(size(x))-.5);
%
%   Then you might prefer to approximate instead. For example, you might try 
%   the cubic smoothing spline, obtained by the command:
%
%   scs = csaps(x,noisy);
%   fnplt(scs), hold on, plot(x,noisy,'o')
%   legend('smoothing spline','noisy data'), hold off


%   If you don't like the level of smoothing done by csaps(x,y), you can change
%   it by specifying the smoothing parameter, p, as an optional third argument.
%   Choose this number anywhere between 0 and 1. As p changes from 0 to 1, the 
%   smoothing spline changes, correspondingly, from one extreme, the least-squares 
%   straight-line approximation to the data, to the other extreme, the "natural" 
%   cubic spline interpolant to the data. Since csaps returns the smoothing 
%   parameter actually used as an optional second output, you could now experiment, 
%   as follows:
%
%   [scs,p] = csaps(x,noisy); fnplt(scs), hold on
%   fnplt(csaps(x,noisy,p/2),'--')
%   fnplt(csaps(x,noisy,(1+p)/2),':'), plot(x,noisy,'o')
%   legend('smoothing spline','more smoothed','less smoothed','noisy data'), hold off


%   At times, you might prefer simply to get the smoothest cubic spline sp that is
%   within a specified tolerance tol of the given data in the sense that:
%
%   norm(noisy - fnval(sp,x))^2 <= tol
%
%   This spline is provided by the command:
%
%   sp = spaps(x,noisy,tol);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   If f is one of the cs, ch, or sp splines then, it can be displayed by the statement
%   fnplt(f)
%   Its value at a is given by the statement:
%   fnval(f,a);
%   Its second derivative is constructed by the statement:
%   DDf = fnder(fnder(f));
%   or by the statement:
%   DDf = fnder(f,2);
%   Its definite integral over the interval [a..b] is supplied by the statement:
%   [1 -1]*fnval(fnint(f),[b;a]);
%   and the difference between the spline in cs and the one in ch can be computed as
%   fncmb(cs,'-',sp);
%
%
%   The toolbox supports vector-valued splines. For example, if you want a spline curve 
%   through given planar points (x(i), y(i)), i = 1, ..., n, then the statements
%
%   xy = [x;y]; df = diff(xy.').'; 
%   t = cumsum([0, sqrt([1 1]*(df.*df))]); 
%   cv = csapi(t,xy);
%
%   provide such a spline curve, using chord-length parametrization and cubic spline 
%   interpolation with the not-a-knot end condition, as can be verified by the
%   statements
%
%   fnplt(cv), hold on, plot(x,y,'o'), hold off
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Vector-valued splines are also used in the approximation to gridded data, in 
%   any number of variables, using tensor-product splines. The same spline-construction 
%   commands are used, only the form of the input differs. For example, 
%   if x is an m-vector, y is an n-vector, and z is an array of size [m,n], then
%
%   cs = csapi({x,y},z);
%
%   describes a bicubic spline f satisfying f( x(i), y(j) ) = z(i,j) for
%   i = 1:m, j = 1:n. Such a multivariate spline can be vector-valued. For example,
%
%   x = 0:4; y=-2:2; s2 = 1/sqrt(2);
%   z(3,:,:) = [0 1 s2 0 -s2 -1 0].'*[1 1 1 1 1];
%   z(2,:,:) = [1 0 s2 1 s2 0 -1].'*[0 1 0 -1 0];
%   z(1,:,:) = [1 0 s2 1 s2 0 -1].'*[1 0 -1 0 1];
%   sph = csape({x,y},z,{'clamped','periodic'});
%   fnplt(sph), axis equal, axis off
%
%   gives a perfectly acceptable sphere. Its projection onto the (x,z)-plane is plotted by
%
%   fnplt(fncmb(sph,[1 0 0; 0 0 1]))
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear; close all;

x = .0001 + [-4:.2:4];
y = -3:.2:3;

[yy,xx] = meshgrid(y,x);

r = pi * sqrt(xx.^2 + yy.^2);   % sinc function
z = sin(r)./r;

bcs = csapi( {x,y}, z);     % cubic spline (cs) approximate (ap) with interpolation (i)

surf(xx,yy,z)
figure
fnplt(bcs)                  % plot cubic spline interpolation



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 0:4;
b = -2:2;

R = 4;
r = 2;

V(3,:,:) = [0 (R-r)/2 0 (r-R)/2 0].'*[1 1  1  1 1];
V(2,:,:) = [R (r+R)/2 r (r+R)/2 R].'*[0 1  0 -1 0];
V(1,:,:) = [R (r+R)/2 r (r+R)/2 R].'*[1 0 -1  0 1];

dough0 = csape( {a,b}, V, 'periodic');

figure; fnplt(dough0)

% normals to surface
nx = 43; 
xy = [ones(1,nx); linspace(2,-2,nx)];
points = fnval(dough0,xy)';

ders = fnval(fndir(dough0, eye(2)), xy);

normals = cross(ders(4:6,:), ders(1:3,:));
normals = (normals ./ repmat(sqrt(sum(normals .* normals)), 3, 1))';
pn = [points; points + normals];
hold on
for j = 1:nx
    plot3( pn( [j,j+nx], 1), pn( [j,j+nx], 2), pn( [j,j+nx], 3) ),
end
hold off

% 2D projection
figure; fnplt( fncmb( dough0, [1 0 0; 0 1 0] ))

