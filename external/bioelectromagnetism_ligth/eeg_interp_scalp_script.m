% eeg_interp_scalp_script - Script to interpolate scalp potentials over mesh
% 
% This is part of the development and testing of
% functions to implement:
%
% Oostendorp T, Oosterom A, & Huiskamp G (1989),
% Interpolation on a triangulated 3D surface.  
% Journal of Computational Physics, 80: 331-343.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:51 $

% Licence:  GNU GPL, no implied or express warranties
% History:  04/2002, Darren.Weber_at_radiology.ucsf.edu
%           - needs verification and conversion to function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear all

% run R. Oostendorp's script to load simulated data
cd d:\matlab\cvs\eeg_toolbox\lapint
%sphere_load;
%cd ..

% If this script has already been run and the data
% saved, then just load the saved data and return
if isequal(exist('eeg_interp_scalp_mesh_data.mat'),2),
    fprintf('\nLoading saved data.\n');
    load 'eeg_interp_scalp_mesh_data';
    return
end


p = eeg_toolbox_defaults('create');
% Load electrode coordinates 'Cartesian'
p = elec_open(p);
% Load potential values at each electrode
p = eeg_open(p);
% Load a tesselation dataset (scalp, 'skull', cortex)
p = mesh_open(p);

% Find the nearest scalp vertices to each
% electrode vertex
p = mesh_plot(p);
close all

% Extract tesselations from p struct
scalpvert = p.mesh.data.vertices{3};
scalpface = p.mesh.data.faces{3};
elecvert  = p.mesh.data.vertices{4};
elecface  = p.mesh.data.faces{4};

% Extract electric potential at time t for electrodes
Velec = p.volt.data;

clear p;

% This command retains the shape of scalp, while reducing
% the number of faces from ~4000 to 1000.  In doing so, it
% reduces the resolution of the scalp mesh and we need to
% recompute the nearest scalp vertices for each electrode
% vertex below
%[scalpface, scalpvert] = reducepatch(scalpface,scalpvert,1000);

% find the scalp vertex indices for each electrode vertex.
% In this example, this generates replicate indices when
% the reducepatch command above is executed
Scalpindex = dsearchn(scalpvert,elecvert)';

% Calculate the Laplacian matrix
lap = mesh_laplacian(scalpvert,scalpface);

% Calculate interpolation matrix, based on Laplacian
[Lint, keepindex, repindex] = mesh_laplacian_interp(lap,Scalpindex);


% Interpolate potential to scalp
time = 100;
if isempty(repindex),
    Vscalp = Lint * Velec(time,:)';
else
    Vscalp = Lint * Velec(time,keepindex)';
end

%save eeg_interp_scalp_mesh_data


return

% Plot scalp with voltages at all interpolated vertices

fig = figure;
Hp = patch('vertices',scalpvert,'faces',scalpface,'FaceVertexCdata',Vscalp,'facecolor','interp','edgecolor','none');
axis off tight vis3d
lighting gouraud

map = eeg_colormap('Red/Blue/White');
colormap(map)

[az,el] = view;
lit = lightangle(az,el);

set(Hp,'AmbientStrength',.7);
set(Hp,'DiffuseStrength',.3);
set(Hp,'SpecularStrength',0);
set(Hp,'SpecularColorReflectance',0);

%colorbar

% Animation of the timeseries & save movie
set(gcf,'BackingStore','off');
set(Hp,'EraseMode','normal');
figure(fig);
% Define movie region as figure size
rect = get(fig,'Position');
rect(1:2) = [0 0];

clear M

start = 120;
finish = 160;
for t=start:finish,
    Vscalp = Lint * Velec(t,:)';
    set(Hp,'FaceVertexCdata',Vscalp);
    drawnow;
    M(t-(start-1)) = getframe(gcf,rect);
end
% Play the movie in another figure
playtimes = 1;
fps = 15;
figure; axis off; movie(gcf,M,playtimes,fps,rect);
movie2avi(M,'scalp_interpolation.avi','quality',100);



return
