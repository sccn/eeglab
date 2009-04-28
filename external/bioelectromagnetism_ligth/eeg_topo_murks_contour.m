function isobars

% eeg_topo_murks_contour - illustrate isocontours of a surface

clear; close all;

s = 40;
w = 1/150;
cvec = 10:60:300;

p = -s:s;
[x,y] = meshgrid(p);
f = exp(-w*(x.^2 + y.^2));
f = max(f,0);
f = sqrt(f);

surf(x,y,f); shading interp; hold on; rotate3D

V = (x + 10).^2 + (y + 10).^2 + f.^2;

for k = 1:length(cvec)
  c = cvec(k);
  Vc = (V > c);

%  Vcdil = dilate(Vc); Vcdil = double(Vcdil);

  se = [0 1 0 ; 1 1 1 ; 0 1 0];
%  se = ones(1);
  Vcdil = conv2(Vc,se,'same');
  Vcdil = Vcdil > 0;

  Vc = Vcdil - Vc;
  Vcvec = find(Vc);

  plot3(x(Vcvec),y(Vcvec),f(Vcvec),'w.')

end

hold off
