function [p] = mesh_plot(p),

% mesh_plot - plot a mesh from the p structure
%
% Usage: [p] = mesh_plot(p)
%
% The p structure is described in eeg_toolbox_defaults.
%
% Here, the p structure should contain the following fields:
%
% p.mesh.current -  which mesh to plot (integer)
%                   cell of p.mesh.data, see MESH_CHECK
% p.mesh.data    -  contains the vertices and faces in
%                   p.mesh.data.vertices{p.mesh.current}
%                   p.mesh.data.faces{p.mesh.current}
%                   If p.mesh.data is empty, a gui is opened
%                   for loading a mesh data file(s).
%
% This function can plot electrodes on a given scalp mesh
% when p.mesh.data.meshtype contains a 'scalp' mesh
% and p.elec.data is not empty.  However, this function
% does not provide any coregistration/alignment of these
% two datasets and at present it converts the electrode
% data by swapping x/y and converting from centimeters to
% meters (see code for details on why).  It plots given 
% electrodes in blue circles and the nearest vertex 
% points in red dots.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:57 $

% Licence:  GNU GPL, no express or implied warranties
% History:  02/2002, Darren.Weber_at_radiology.ucsf.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('p','var'),[p] = eeg_toolbox_defaults; end

if isfield(p,'mesh'),
  if isfield(p.mesh,'data'),
    if isempty(p.mesh.data),
     [p] = gui_mesh_open(p);
      return;
    end
  else
   [p] = gui_mesh_open(p);
  end
else
 [p] = gui_mesh_open(p);
  return;
end

if isfield(p,'mesh'),
  if isfield(p.mesh,'current'),
    if isempty(p.mesh.current),
      p.mesh.current = 1;
    end
  else
    p.mesh.current = 1;
  end
else
  msg = sprintf('No mesh field in p struct.\n');
  error(msg);
end


% -- 'Coregister' electrodes, if required

if strmatch('scalp',p.mesh.data.meshtype{p.mesh.current}) > 0,
  if ~isempty(p.elec.data),
   [p] = mesh_fit_elec(p);
    p.mesh.current = mesh_check(p,'scalp');
  end
end


% -- Get the vertices and faces

vertices = p.mesh.data.vertices{p.mesh.current};
faces    = p.mesh.data.faces{p.mesh.current};

if isempty(vertices),
  fprintf('...p.mesh.current points to an empty cell\n');
  return
end

% -- Get the figure or create it

if isfield(p.mesh.plot,'fig'),
  figure(p.mesh.plot.fig);
else
  p.mesh.plot.fig = figure;
end

set(gcf,'color',[0 0 0]);

set(gca,'Projection','perspective')
%set(gca,'Projection','orthographic')
set(gca,'DataAspectRatio',[1 1 1]);


% -- Get the handles to the current patch

if isfield(p.mesh.plot,'Hpatch'),
  if p.mesh.plot.overlay,
    % leave the current plot as is
    hold on;
  else
    % clear all plot handles
    hold off;
    
    for i = 1:length(p.mesh.plot.Hpatch),
      if ishandle(p.mesh.plot.Hpatch{i}),
        delete(p.mesh.plot.Hpatch{i});
      end
    end
    p.mesh.plot.Hpatch = cell(0);
  end
end
p.mesh.plot.Hpatch{p.mesh.current} = patch('Vertices',vertices,'Faces',faces);





axis off tight vis3d


set(gcf,'Renderer','opengl'); lighting phong

% Create the light source
if isfield(p.mesh.plot,'Hlight'),
  for i = 1:length(p.mesh.plot.Hlight),
    if ishandle(p.mesh.plot.Hlight(i)),
      delete(p.mesh.plot.Hlight(i));
    end
  end
  p.mesh.plot.Hlight = [];
end

H.light(1) = camlight(  0,  0,'local'); % top
H.light(2) = camlight(-75,  0); % left
H.light(3) = camlight( 75,  0); % right
%H.light(4) = camlight(  0,-90,'infinite'); % back
%H.light(5) = camlight(  0, 90,'infinite'); % front
%H.light(6) = camlight(  0,180,'infinite'); % bottom

% for i = 1:length(H.light),
%   % mute the intensity of the lights
%   color = [.8 1 1] / length(H.light) / 1.2;
%   set(H.light(i),'color',color);
%   %get(H.light(i),'color')
% end
p.mesh.plot.Hlight = H.light;


% MATERIAL([ka kd ks n sc]) sets the ambient/diffuse/specular strength,
%    specular exponent and specular color reflectance of the objects.
%reflect = material('dull');
if ~isfield(p.mesh.plot,'reflect')
  p.mesh.plot.ambient = 0.5;
  p.mesh.plot.diffuse = 0.4;
  p.mesh.plot.specular = 0.05;
  p.mesh.plot.specexp = 5;
  p.mesh.plot.speccolor = 1;
end
if ~isfield(p.mesh.plot,'edgecolor'),
  % plot the cortex without edgecolor
  switch p.mesh.data.meshtype{p.mesh.current},
    case 'cortex',      p.mesh.plot.edgecolor = 'none';
    case 'pial',        p.mesh.plot.edgecolor = 'none';
    case 'smoothwm',    p.mesh.plot.edgecolor = 'none';
    case 'white',       p.mesh.plot.edgecolor = 'none';
    otherwise,          p.mesh.plot.edgecolor = '.7 .7 .7';
  end
end

if ~isfield(p.mesh.plot,'facecolor'),
  p.mesh.plot.facecolor = '.9 .8 .7';
end
if ~isfield(p.mesh.plot,'facealpha'),
  p.mesh.plot.facealpha = '1';
end

if findstr(p.mesh.plot.edgecolor,'none'),
  set(p.mesh.plot.Hpatch{p.mesh.current},'EdgeColor','none');
else
  set(p.mesh.plot.Hpatch{p.mesh.current},'EdgeColor',str2num(p.mesh.plot.edgecolor));
end
if findstr(p.mesh.plot.facecolor,'none'),
  set(p.mesh.plot.Hpatch{p.mesh.current},'FaceColor','none');
else
  set(p.mesh.plot.Hpatch{p.mesh.current},'FaceColor',str2num(p.mesh.plot.facecolor));
end
if p.mesh.plot.facealpha,
  set(p.mesh.plot.Hpatch{p.mesh.current},'FaceAlpha',str2num(p.mesh.plot.facealpha));
else
  set(p.mesh.plot.Hpatch{p.mesh.current},'FaceAlpha',1);
end
set(p.mesh.plot.Hpatch{p.mesh.current},'AmbientStrength',p.mesh.plot.ambient);
set(p.mesh.plot.Hpatch{p.mesh.current},'DiffuseStrength',p.mesh.plot.diffuse);
set(p.mesh.plot.Hpatch{p.mesh.current},'SpecularStrength',p.mesh.plot.specular);
set(p.mesh.plot.Hpatch{p.mesh.current},'SpecularExponent',p.mesh.plot.specexp);
set(p.mesh.plot.Hpatch{p.mesh.current},'SpecularColorReflectance',p.mesh.plot.speccolor);


% -- Plot the electrodes, if scalp and electrodes available
if strmatch('scalp',p.mesh.data.meshtype{p.mesh.current}) > 0,
  if ~isempty(p.elec.data),
   [p] = elec_plot(p);
  end
end



% -- setup the viewer rotation command

if isequal(exist('mouse_rotate'),2),
  userdata = get(p.mesh.plot.fig,'userdata');
  if isempty(userdata),
    mouse_rotate('init')
  else
    mouse_rotate('init',get(p.mesh.plot.fig,'userdata'));
  end
else
  rotate3d;
end


if nargout < 1,[p] = []; end

return







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[p] = elec_plot(p),

elecN = mesh_check(p,'elec');

if isempty(elecN),
  warning('MESH_FIT_ELEC: No ''elec'' mesh in p struct.\n');
end

% Get the scalp vertices and faces
% that correspond to the electrodes
elecvert  = p.mesh.data.vertices{elecN};

x = elecvert(:,1);
y = elecvert(:,2);
z = elecvert(:,3);

hold on
scatter3(x,y,z,70,'b','filled');

%     % add electrode labels
%     for e = 1:length(x),
%         % relocate xyz outside of mesh
%         c = [0 0 0];
%         r = sqrt(x(e).^2 + y(e).^2 + z(e).^2); r = r + (r .* 0.05);
%         v = sphere_project([x(e),y(e),z(e)],r,c);
%         H(e)  = text(v(1),v(2),v(3),p.elec.data.label{e},...
%                      'HorizontalAlignment','center',...
%                      'VerticalAlignment','middle');
%     end


%     % plot regions of electrodes in different colors
%     if isfield(p.elec.data,'regions'),
%         if ~isempty(p.elec.data.regions),
%             elecregions = p.elec.data.regions;
%             C = zeros(length(x),1); % generate color array for electrodes
%             colors = rand(length(x),1);
%             for r = 1:length(elecregions),
%                 % assign random colors to electrode indices
%                 C(elecregions(r).elec) = colors(r);
%                 %C(elecregions(r).elec) = r;
%             end
%             map = colormap(colorcube);
%             %map = colormap(hsv);
%             %map = colormap(lines);
%             scatter3(x,y,z,80,C,'b','filled');
%             
% %             % attempt to plot an interpolation of these regions
% %             p.colorMap.map = map;
% %             p.volt.data = repmat(C,1,2);
% %             p.volt.timeArray = [1 2];
% %             p.volt.points    = 2
% %             p.volt.sampleHz  = 1000;
% %             p.volt.sampleRate = 1;
% %             p.volt.sampleTime = 1;
% %             p.volt.channels   = length(x);
% %             p.volt.epochStart = 1;
% %             p.volt.epochEnd   = 2;
% %             p.elec.plotSurf     = 0;
% %             p.mesh.plotSurf     = 1;
% %             p.clickTimePoint    = 0;
% %             p.volt.samplePoint         = 1;
% %             eeg_contours_engine(p)
%         end
%     end

hold off
return
