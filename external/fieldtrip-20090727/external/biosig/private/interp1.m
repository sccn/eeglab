## Copyright (C) 2000 Paul Kienzle
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

## usage: yi = interp1(x, y, xi [, 'method' [, 'extrap']])
##
## Interpolate the function y=f(x) at the points xi. The sample 
## points x must be strictly monotonic.  If y is a matrix with
## length(x) rows, yi will be a matrix of size rows(xi) by columns(y),
## or its transpose if xi is a row vector.
##
## Method is one of:
## 'nearest': return nearest neighbour.
## 'linear': linear interpolation from nearest neighbours
## 'pchip': piece-wise cubic hermite interpolating polynomial
## 'cubic': cubic interpolation from four nearest neighbours
## 'spline': cubic spline interpolation--smooth first and second
##           derivatives throughout the curve
## ['*' method]: same as method, but assumes x is uniformly spaced
##               only uses x(1) and x(2); usually faster, never slower
##
## Method defaults to 'linear'.
##
## If extrap is the string 'extrap', then extrapolate values beyond
## the endpoints.  If extrap is a number, replace values beyond the
## endpoints with that number.  If extrap is missing, assume NaN.
##
## Example:
##    xf=[0:0.05:10]; yf = sin(2*pi*xf/5);
##    xp=[0:10];      yp = sin(2*pi*xp/5);
##    lin=interp1(xp,yp,xf);
##    spl=interp1(xp,yp,xf,'spline');
##    cub=interp1(xp,yp,xf,'cubic');
##    near=interp1(xp,yp,xf,'nearest');
##    plot(xf,yf,';original;',xf,lin,';linear;',xf,spl,';spline;',...
##         xf,cub,';cubic;',xf,near,';nearest;',xp,yp,'*;;');
##
## See also: interp

## 2000-03-25 Paul Kienzle
##    added 'nearest' as suggested by Kai Habel
## 2000-07-17 Paul Kienzle
##    added '*' methods and matrix y
##    check for proper table lengths
## 2002-01-23 Paul Kienzle
##    fixed extrapolation

function yi = interp1(x, y, xi, method, extrap)

  if nargin<3 || nargin>5
    usage("yi = interp1(x, y, xi [, 'method' [, 'extrap']])");
  endif

  if nargin < 4, 
    method = 'linear';
  else
    method = tolower(method); 
  endif

  if nargin < 5
    extrap = NaN;
  endif

  ## reshape matrices for convenience
  x = x(:);
  if size(y,1)==1, y=y(:); endif
  transposed = (size(xi,1)==1);
  xi = xi(:);

  ## determine sizes
  nx = size(x,1);
  [ny, nc] = size(y);
  if (nx < 2 || ny < 2)
     error ("interp1: table too short");
  endif

  ## determine which values are out of range and set them to extrap,
  ## unless extrap=='extrap' in which case, extrapolate them like we
  ## should be doing in the first place.
  minx = x(1);
  if (method(1) == '*')
     dx = x(2) - x(1);
     maxx = minx + (ny-1)*dx;
  else
     maxx = x(nx);
  endif
  if strcmp(extrap,"extrap")
    range=1:size(xi,1);
    yi = zeros(size(xi,1), size(y,2));
  else
    range = find(xi >= minx & xi <= maxx);
    yi = extrap*ones(size(xi,1), size(y,2));
    if isempty(range), 
      if transposed, yi = yi.'; endif
      return; 
    endif
    xi = xi(range);
  endif

  if strcmp(method, 'nearest')
    idx = lookup(0.5*(x(1:nx-1)+x(2:nx)), xi)+1;
    yi(range,:) = y(idx,:);

  elseif strcmp(method, '*nearest')
    idx = floor((xi-minx)/dx+1.5);
    yi(range,:) = y(idx,:);

  elseif strcmp(method, 'linear')
    ## find the interval containing the test point
    idx = lookup (x(2:nx-1), xi)+1; 
				# 2:n-1 so that anything beyond the ends
				# gets dumped into an interval
    ## use the endpoints of the interval to define a line
    dy = y(2:ny,:) - y(1:ny-1,:);
    dx = x(2:nx) - x(1:nx-1);
    s = (xi - x(idx))./dx(idx);
    yi(range,:) = s(:,ones(1,nc)).*dy(idx,:) + y(idx,:);

  elseif strcmp(method, '*linear')
    ## find the interval containing the test point
    t = (xi - minx)/dx + 1;
    idx = floor(t);

    ## use the endpoints of the interval to define a line
    dy = [y(2:ny,:) - y(1:ny-1,:); zeros(1,nc)];
    s = (t - idx)./dx;
    yi(range,:) = s(:,ones(1,nc)).*dy(idx,:) + y(idx,:); 

  elseif strcmp(method, 'pchip') || strcmp(method, '*pchip')
    if (nx == 2) x = linspace(minx, maxx, ny); endif
    yi(range,:) = pchip(x, y, xi);

  elseif strcmp(method, 'cubic')
    if (nx < 4 || ny < 4)
      error ("interp1: table too short");
    endif
    idx = lookup(x(3:nx-2), xi) + 1;

    ## Construct cubic equations for each interval using divided
    ## differences (computation of c and d don't use divided differences
    ## but instead solve 2 equations for 2 unknowns). Perhaps
    ## reformulating this as a lagrange polynomial would be more efficient.
    i=1:nx-3;
    J = ones(1,nc);
    dx = diff(x);
    dx2 = x(i+1).^2 - x(i).^2;
    dx3 = x(i+1).^3 - x(i).^3;
    a=diff(y,3)./dx(i,J).^3/6;
    b=(diff(y(1:nx-1,:),2)./dx(i,J).^2 - 6*a.*x(i+1,J))/2;
    c=(diff(y(1:nx-2,:),1) - a.*dx3(:,J) - b.*dx2(:,J))./dx(i,J);
    d=y(i,:) - ((a.*x(i,J) + b).*x(i,J) + c).*x(i,J);
    yi(range,:) = ((a(idx,:).*xi(:,J) + b(idx,:)).*xi(:,J) ...
		   + c(idx,:)).*xi(:,J) + d(idx,:);

  elseif strcmp(method, '*cubic')
    if (nx < 4 || ny < 4)
      error ("interp1: table too short");
    endif

    ## From: Miloje Makivic 
    ## http://www.npac.syr.edu/projects/nasa/MILOJE/final/node36.html
    t = (xi - minx)/dx + 1;
    idx = max(min(floor(t), ny-2), 2);
    t = t - idx;
    t2 = t.*t;
    tp = 1 - 0.5*t;
    a = (1 - t2).*tp;
    b = (t2 + t).*tp;
    c = (t2 - t).*tp/3;
    d = (t2 - 1).*t/6;
    J = ones(1,nc);
    yi(range,:) = a(:,J) .* y(idx,:) + b(:,J) .* y(idx+1,:) ...
		  + c(:,J) .* y(idx-1,:) + d(:,J) .* y(idx+2,:);

  elseif strcmp(method, 'spline') || strcmp(method, '*spline')
    if (nx == 2) x = linspace(minx, maxx, ny); endif
    yi(range,:) = spline(x, y, xi);

  else
    error(["interp1 doesn't understand method '", method, "'"]);
  endif
  if transposed, yi=yi.'; endif

endfunction

%!demo
%! xf=0:0.05:10; yf = sin(2*pi*xf/5);
%! xp=0:10;      yp = sin(2*pi*xp/5);
%! lin=interp1(xp,yp,xf,'linear');
%! spl=interp1(xp,yp,xf,'spline');
%! cub=interp1(xp,yp,xf,'pchip');
%! near=interp1(xp,yp,xf,'nearest');
%! plot(xf,yf,';original;',xf,near,';nearest;',xf,lin,';linear;',...
%!      xf,cub,';pchip;',xf,spl,';spline;',xp,yp,'*;;');
%! %--------------------------------------------------------
%! % confirm that interpolated function matches the original

%!demo
%! xf=0:0.05:10; yf = sin(2*pi*xf/5);
%! xp=0:10;      yp = sin(2*pi*xp/5);
%! lin=interp1(xp,yp,xf,'*linear');
%! spl=interp1(xp,yp,xf,'*spline');
%! cub=interp1(xp,yp,xf,'*cubic');
%! near=interp1(xp,yp,xf,'*nearest');
%! plot(xf,yf,';*original;',xf,near,';*nearest;',xf,lin,';*linear;',...
%!      xf,cub,';*cubic;',xf,spl,';*spline;',xp,yp,'*;;');
%! %--------------------------------------------------------
%! % confirm that interpolated function matches the original

%!shared xp, yp, xi, style
%! xp=0:5;      yp = sin(2*pi*xp/5);
%! xi = sort([-1, max(xp)*rand(1,6), max(xp)+1]);

%!test style = 'nearest';
%!assert (interp1(xp, yp, [min(xp)-1, max(xp)+1]), [NaN, NaN]);
%!assert (interp1(xp,yp,xp,style), yp, 100*eps);
%!assert (interp1(xp,yp,xp',style), yp', 100*eps);
%!assert (interp1(xp',yp',xp',style), yp', 100*eps);
%!assert (interp1(xp',yp',xp,style), yp, 100*eps);
%!assert (isempty(interp1(xp',yp',[],style)));
%!assert (isempty(interp1(xp,yp,[],style)));
%!assert (interp1(xp,[yp',yp'],xi(:),style),...
%!	  [interp1(xp,yp,xi(:),style),interp1(xp,yp,xi(:),style)]);
%!assert (interp1(xp,[yp',yp'],xi,style),
%!	  interp1(xp,[yp',yp'],xi,["*",style]));

%!test style = 'linear';
%!assert (interp1(xp, yp, [-1, max(xp)+1]), [NaN, NaN]);
%!assert (interp1(xp,yp,xp,style), yp, 100*eps);
%!assert (interp1(xp,yp,xp',style), yp', 100*eps);
%!assert (interp1(xp',yp',xp',style), yp', 100*eps);
%!assert (interp1(xp',yp',xp,style), yp, 100*eps);
%!assert (isempty(interp1(xp',yp',[],style)));
%!assert (isempty(interp1(xp,yp,[],style)));
%!assert (interp1(xp,[yp',yp'],xi(:),style),...
%!	  [interp1(xp,yp,xi(:),style),interp1(xp,yp,xi(:),style)]);
%!assert (interp1(xp,[yp',yp'],xi,style),
%!	  interp1(xp,[yp',yp'],xi,["*",style]),100*eps);

%!test style = 'cubic';
%!assert (interp1(xp, yp, [-1, max(xp)+1]), [NaN, NaN]);
%!assert (interp1(xp,yp,xp,style), yp, 100*eps);
%!assert (interp1(xp,yp,xp',style), yp', 100*eps);
%!assert (interp1(xp',yp',xp',style), yp', 100*eps);
%!assert (interp1(xp',yp',xp,style), yp, 100*eps);
%!assert (isempty(interp1(xp',yp',[],style)));
%!assert (isempty(interp1(xp,yp,[],style)));
%!assert (interp1(xp,[yp',yp'],xi(:),style),...
%!	  [interp1(xp,yp,xi(:),style),interp1(xp,yp,xi(:),style)]);
%!assert (interp1(xp,[yp',yp'],xi,style),
%!	  interp1(xp,[yp',yp'],xi,["*",style]),1000*eps);

%!test style = 'spline';
%!assert (interp1(xp, yp, [-1, max(xp) + 1]), [NaN, NaN]);
%!assert (interp1(xp,yp,xp,style), yp, 100*eps);
%!assert (interp1(xp,yp,xp',style), yp', 100*eps);
%!assert (interp1(xp',yp',xp',style), yp', 100*eps);
%!assert (interp1(xp',yp',xp,style), yp, 100*eps);
%!assert (isempty(interp1(xp',yp',[],style)));
%!assert (isempty(interp1(xp,yp,[],style)));
%!assert (interp1(xp,[yp',yp'],xi(:),style),...
%!	  [interp1(xp,yp,xi(:),style),interp1(xp,yp,xi(:),style)]);
%!assert (interp1(xp,[yp',yp'],xi,style),
%!	  interp1(xp,[yp',yp'],xi,["*",style]),10*eps);

%!error interp1
%!error interp1(1:2,1:2,1,'bogus')

%!error interp1(1,1,1, 'nearest');
%!assert (interp1(1:2,1:2,1.4,'nearest'),1);
%!error interp1(1,1,1, 'linear');
%!assert (interp1(1:2,1:2,1.4,'linear'),1.4);
%!error interp1(1:3,1:3,1, 'cubic');
%!assert (interp1(1:4,1:4,1.4,'cubic'),1.4);
%!error interp1(1:2,1:2,1, 'spline');
%!assert (interp1(1:3,1:3,1.4,'spline'),1.4);

%!error interp1(1,1,1, '*nearest');
%!assert (interp1(1:2,1:2,1.4,'*nearest'),1);
%!error interp1(1,1,1, '*linear');
%!assert (interp1(1:2,1:2,1.4,'*linear'),1.4);
%!error interp1(1:3,1:3,1, '*cubic');
%!assert (interp1(1:4,1:4,1.4,'*cubic'),1.4);
%!error interp1(1:2,1:2,1, '*spline');
%!assert (interp1(1:3,1:3,1.4,'*spline'),1.4);
