function[p] = eeg_plot_surf(p)

% eeg_plot_surf - Plot a patch surface with voltage colourmap
% 
%[p] = eeg_plot_surf(p)
% 
% fields of p struct used:
%   
%   p.rangeMethod   - initialised by eeg_toolbox_defaults
%   p.colorMap.map  - initialised by eeg_toolbox_defaults
%   
%   p.volt.sampleTime
%   p.mesh.data.vertices{p.mesh.current}
%   p.mesh.data.faces{p.mesh.current}
%   p.mesh.data.timeseries{p.mesh.current}
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:52 $

% Licence:  GNU GPL, no express or implied warranties
% History:  10/2002, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nEEG_PLOT_SURF...\n'); tic;


% confirm that a patch surface is available
if isfield(p,'mesh'),
  if isfield(p.mesh,'data'),
    if isfield(p.mesh.data,'timeseries'),
      if isempty(p.mesh.data.Cdata{p.mesh.current}),
        msg = sprintf('...p.mesh.data.Cdata{%d} is empty\n',p.mesh.current);
        error(msg);
      end
    end
  end
else
  error('...p.mesh.data is empty - load mesh first\n');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine/validate the min, max surface color data range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select mesh timeseries values at samplePoint
switch p.mesh.data.meshtype{p.mesh.current},
  case {'scalp','elec'},
    samplePoint = p.volt.samplePoint;
    sampleTime  = p.volt.sampleTime;
  otherwise
    samplePoint = p.mesh.samplePoint;
    sampleTime  = p.mesh.sampleTime;
end

% get number of vertices
nvert = size(p.mesh.data.vertices{p.mesh.current},1);
% Assume more vertices than time points
[s1,s2] = size(p.mesh.data.Cdata{p.mesh.current});
if isequal(nvert,s1), % vertices in rows, timepoints in columns
  meshCdata = p.mesh.data.Cdata{p.mesh.current}(:,samplePoint);
else,                 % vertices in columns, timeseries in rows
  meshCdata = p.mesh.data.Cdata{p.mesh.current}(samplePoint,:)';
end

% Use absolute values for cortical surfaces
switch p.mesh.data.meshtype{p.mesh.current},
  case {'scalp','elec'},
  otherwise
    meshCdata = abs(meshCdata);
end


switch p.rangeMethod
  case 'minmaxall', % Min/Max,all points
    fprintf('...estimating color data range, min/max all time points.\n');
    p.maximumIntensity = max(max(p.mesh.data.Cdata{p.mesh.current}));
    p.minimumIntensity = min(min(p.mesh.data.Cdata{p.mesh.current}));
  case 'minmaxone', % Min/Max, single point
    fprintf('...estimating color data range, min/max single time point.\n');
    % get number of vertices
    p.maximumIntensity = max(max(meshCdata));
    p.minimumIntensity = min(min(meshCdata));
  case 'minmaxabs', % Min/Max, Absolute
    fprintf('...estimating color data range, abs min/max single time point.\n');
    absmax = max(max(abs(meshCdata)));
    p.maximumIntensity =  absmax;
    p.minimumIntensity = -absmax;
  otherwise
    % check that specified intensity range is defined
    fprintf('...checking predefined color data range.\n');
    if isempty(p.maximumIntensity),
      fprintf('...estimating color data range, min/max single time point.\n');
      p.maximumIntensity = max(max(meshCdata)); end
    if isempty(p.minimumIntensity),
      fprintf('...estimating color data range, min/max single time point.\n');
      p.minimumIntensity = min(min(meshCdata)); end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the surface patch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('...plotting the surface.\n');


% Set the window title name
if sampleTime,
  if p.volt.data,
    data = sprintf('%s %s',p.volt.file,p.mesh.data.meshtype{p.mesh.current});
    name = sprintf('Surface: %s @ %8.2f msec',data, sampleTime);
  else
    data = sprintf('%s',p.mesh.data.meshtype{p.mesh.current});
    name = sprintf('Surface: %s @ %8.2f msec',data, sampleTime);
  end
else
  if p.volt.data,
    data = sprintf('%s %s',p.volt.file,p.mesh.data.meshtype{p.mesh.current});
    name = sprintf('Surface: %s',data);
  else
    data = sprintf('%s',p.mesh.data.meshtype{p.mesh.current});
    name = sprintf('Surface: %s',data);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot it

% create the figure window
H.gui = figure('NumberTitle','off','Name',name,'color',[0 0 0]);

% Create the patch
p.mesh.patch = patch('vertices',p.mesh.data.vertices{p.mesh.current},...
  'faces',p.mesh.data.faces{p.mesh.current},...
  'facevertexCdata',meshCdata,...
  'facecolor','interp','EdgeColor','none',...
  'FaceLighting','phong');

%data = interp3(Vscalp,3,'cubic');
%isonormals(Vscalp,p.mesh.patch)

set(gca,'Projection','perspective')
set(gca,'DataAspectRatio',[1 1 1]);
axis off tight vis3d

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lighting

lighting phong
set(H.gui,'RendererMode','auto')

% Create a light above the z axis
%H.light = light('style','infinite','position',[0 0 10000]);

H.light(1) = camlight(-20, 30,'infinite');
H.light(2) = camlight( 20, 30,'infinite');
H.light(3) = camlight(-20,-30,'infinite');
for i = 1:length(H.light),
  % mute the intensity of the lights
  set(H.light(i),'color',[.8 1 1]/length(H.light)/1.2);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Colormap

meshtype = lower(p.mesh.data.meshtype{p.mesh.current});
switch meshtype,
  case 'scalp',
    set(gca,'AmbientLightColor',[.9 .8 .7]);  % skin color
    colormap(p.colorMap.map);
  case {'cortex','pial','white','smoothwm','orig'},
    % Ambient, Diffusion, Specular, spec exp, spec reflect (see material)
    %MaterialBack =   [0.5 0.7 0.3 10 1]; % qualities of the background image
    %MaterialSource = [0.8 0.7 0.8 20 1]; % qualities for the data
    
    if exist('bluehot'),
      p.colorMap.map = grayish(bluehot(128),.33); % brainstorm cmap
      p.colorMap.map = grayish(hot(128),.33); % brainstorm cmap
      p.colorMap.map = hot(128);
      colormap(p.colorMap.map);
    else
      set(gca,'AmbientLightColor',[.6 .6 .6]);  % grey color
      colormap(p.colorMap.map);
    end
    
  otherwise
    set(gca,'AmbientLightColor',[.9 .9 .9]);  % white light
end

caxis([p.minimumIntensity p.maximumIntensity]);
colorbar


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface Properties

p.reflect{1} = 0.9;
p.reflect{2} = 0.2;
p.reflect{3} = 0.0;
p.reflect{4} = 500;
p.reflect{5} = 0;
set(p.mesh.patch,'AmbientStrength',p.reflect{1});
set(p.mesh.patch,'DiffuseStrength',p.reflect{2});
set(p.mesh.patch,'SpecularStrength',p.reflect{3});
set(p.mesh.patch,'SpecularExponent',p.reflect{4});
set(p.mesh.patch,'SpecularColorReflectance',p.reflect{5});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface Contours

if (p.contour.plot3D == 1),
 [p] = eeg_plot_surf_contours(p);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional controls

H.p = p;

set(H.gui,'userdata',H);
set(gca,'Visible','off');

if isequal(exist('mouse_rotate'),2),
  mouse_rotate('init',H);
else
  rotate3d on;
end
if isequal(exist('gui_topo_animate'),2),
  gui_topo_animate('init',p);
end


t = toc;
fprintf('...done (%6.2f sec).\n\n',t);

return
