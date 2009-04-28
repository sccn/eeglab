% eeg_interp_scalp3d_test - Script to test 3D interpolation of scalp potentials

clear

p = eeg_toolbox_defaults;
p = eeg_open(p);
p = elec_open(p);
p = mesh_open(p);

p = mesh_plot(p);

close all

scalpvert = p.mesh.data.vertices{3};

Mean = mean(scalpvert);
STD  = std(scalpvert);

Min = Mean - 3*STD;
Max = Mean + 3*STD;

delta = 0.05;
Xd = Min(1):delta:Max(1);
Yd = Min(2):delta:Max(2);
Zd = Min(3):delta:Max(3);

[x0,y0,z0] = meshgrid(Xd,Yd,Zd);

X0 = [x0(:), y0(:), z0(:)];

X = p.mesh.data.vertices{4};
v = p.volt.data(1,:)';

v0 = griddatan(X,v,X0);
v0 = reshape(v0, size(x0));

V = v0(:);
Vfinite = find(isfinite(V));
X0finite = X0(Vfinite,:);

trisurf(convhulln(X0finite),X0finite(:,1),X0finite(:,2),X0finite(:,3),Vfinite,'EdgeColor','none','FaceColor','interp');

fprintf('\ndone\n');
return


%X0scalpindices = dsearchn(X0,scalpvert);
X0scalpindices = dsearchn(X0finite,scalpvert);

newscalpvert(:,1) = x0(X0scalpindices);
newscalpvert(:,2) = y0(X0scalpindices);
newscalpvert(:,3) = z0(X0scalpindices);
newvoltage = v0(X0scalpindices);

trisurf(convhulln(newscalpvert),newscalpvert(:,1),newscalpvert(:,2),newscalpvert(:,3),newvoltage,'EdgeColor','none','FaceColor','interp');

return

Vplot = mean(v);

Hp = patch(isosurface(x0,y0,z0,v0,Vplot)); 
isonormals(x0,y0,z0,v0,Hp); 
set(Hp,'FaceColor','red','EdgeColor','none'); 
view(3); 
camlight; 
lighting phong
%axis equal
title('Interpolated isosurface from scattered data')
