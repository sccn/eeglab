function [p] = eeg_contours_engine(p)

% eeg_contours_engine - Generates contour/topographic maps
% 
% Usage: [p] = eeg_contours_engine(p)
%
% p is the eeg_toolbox struct (see eeg_toolbox_defaults)
%
% Given:
%   data file of electrode positions (see elec_load for help)
%   data file of electrode voltage values (see eeg_ascii_load for help)
%   other parameters for selection and plot control
%
%   Parameters can be passed as a structure (p) with the function call or
%   passed in the file 'eeg_contours_default.txt'.  A standard set of
%   parameters can be created using 'eeg_contours_create_defaults'. See
%   that function for more information about the parameters created.
%
% Output:
%	various topographic maps and 2D or 3D contour plots
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:51 $

% Licence:  Gnu GPL, no express or implied warranties
% History:  08/99 Chris Harvey,	Created, combining work by Darren Weber, Murk Bottomer
%           10/99 CRH	Accepts input parameters, provides some validation
%           07/01 Darren.Weber_at_radiology.ucsf.edu
%               - modified parameter handling, now using 'p' structure
%               - now calls separate functions to create/read/write defaults
%               - modified most aspects of code to adapt to development of other
%                 functions called from this function
% Develop:  08/01 Darren.Weber_at_radiology.ucsf.edu
%               - Needs scripting facility to wind through time points in
%                 a range of timepoints, outputing contour plot(s) for each in
%                 a given image format (EPS, GIF, TIF, etc.) or a movie of
%                 the timeseries.
%           05/02 Darren.Weber_at_radiology.ucsf.edu
%               - now integrated with mesh interpolation and calls
%                 animation interface, with image saving options
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('p','var'),
  fprintf('Setting default parameters.\n');
 [p] = eeg_toolbox_defaults('create');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load electrode co-ordinates 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(p.elec.data),[p] = elec_open(p); end

%elecLabels = p.elec.data.label;
X = p.elec.data.x;
Y = p.elec.data.y;
Z = p.elec.data.z;

Xrad = p.elec.data.R(1);
Yrad = p.elec.data.R(2);
Zrad = p.elec.data.R(3);

elecNumElectrodes = length(p.elec.data.x);

% Estimate sphere or ellipse that best fits electrode co-ordinates
switch p.elec.shape
  case 'ellipse'
    Xp = p.elec.data.Xel;
    Yp = p.elec.data.Yel;
    Zp = p.elec.data.Zel;
  case 'sphere'
    Xp = p.elec.data.Xsp;
    Yp = p.elec.data.Ysp;
    Zp = p.elec.data.Zsp;
  otherwise
    Xp = p.elec.data.x;
    Yp = p.elec.data.y;
    Zp = p.elec.data.z;           
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read voltage data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(p.volt.data),[p] = eeg_open(p); end

[voltNumTimePoints, voltNumElectrodes] = size(p.volt.data);

% Validate voltage matrix orientation against electrode file

if (voltNumElectrodes ~= elecNumElectrodes)
  fprintf('\nWarning: # electrodes, data file %d, electrode file %d\n', voltNumElectrodes, elecNumElectrodes);
  if (voltNumTimePoints == elecNumElectrodes)
    % Assume data file needs rotation
    fprintf('Continuing with data file rotated\n');
    p.volt.data = p.volt.data';   [voltNumTimePoints, voltNumElectrodes] = size(p.volt.data);
  else
    fprintf('\nError: Cannot reconcile electrodes with voltage data.\n');
    return
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select a timePoint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if or(isequal(p.clickTimePoint, 1),isempty(p.volt.samplePoint))
  figure('NumberTitle','off','Name','MouseClick the TimePoint')
  plot(p.volt.data); 
  [Xt,Yt] = ginput(1);
  p.volt.samplePoint = round(Xt);
  close
end
if (p.volt.samplePoint > voltNumTimePoints)
  msg = sprintf(' Selected timepoint exceeds data points of %d\n\n', voltNumTimePoints);
  warning(msg);
  p.volt.samplePoint = voltNumTimePoints;
end
if ~isempty(p.volt.timeArray),
  p.volt.sampleTime = p.volt.timeArray(p.volt.samplePoint);
end

V = p.volt.data(p.volt.samplePoint,:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine/validate the min, max range and the step size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch p.rangeMethod
  case 'minmaxall', % Min/Max,all points
    p.maximumIntensity = max(max(p.volt.data));
    p.minimumIntensity = min(min(p.volt.data));
  case 'minmaxone', % Min/Max, single point
    p.maximumIntensity = max(max(V));
    p.minimumIntensity = min(min(V));
  case 'minmaxabs', % Min/Max, Absolute
    absmax = max(max(abs(V)));
    p.maximumIntensity =  absmax;
    p.minimumIntensity = -absmax;
  otherwise
    % check that specified intensity range is defined
    if isempty(p.maximumIntensity),
      p.maximumIntensity = max(max(V)); end
    if isempty(p.minimumIntensity),
      p.minimumIntensity = min(min(V)); end
end

switch p.contour.stepMethod
  case 0  % StepSize method
    p.contour.levels = eeg_contour_levels( p.contour.stepSize, V );
    p.contour.Nsteps = length(p.contour.levels);
  otherwise % Number of steps method
    p.contour.stepSize = abs(p.maximumIntensity - p.minimumIntensity) / p.contour.Nsteps;
    p.contour.levels = eeg_contour_levels( p.contour.stepSize, V );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw the plots
%	X,Y,Z		    = actual electrode placements
% 	Xel,Yel,Zel 	= electrode position on the ellipsoid surface
%	gridSize	    = number of grid lines to use on each plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine interpolation gridding
m = max(abs([Xp;Yp]));
m = ceil(m);
if (p.grid.method == 2), xg = -m:p.grid.res:m;             % grid resolution method
else,                   xg = linspace(-m,m,p.grid.size);  % grid size method
end
yg = xg';
[Xi,Yi,Zi] = griddata(Xp,Yp,Zp,xg,yg,p.interpMethod);    % Surface interpolation
Vi         = griddata(Xp,Yp,V ,xg,yg,p.interpMethod);    % Voltage interpolation

% Calculate voltage at refined electrode vertices
%if ~isempty(p.elec.data.Vsp),
%    p.elec.data.Psp = p.elec.data.Isp * V;
%end


% Calculate voltage at scalp mesh vertices
if p.mesh.plotSurf,
  
  % get the scalp mesh
  [p.mesh.current,meshExists] = mesh_check(p,'scalp');
  
  if isempty(p.mesh.current),
    fprintf('...no scalp mesh, using electrode surface.\n');
    p.mesh.plotSurf = 0;
    p.elec.plotSurf = 1;
  else
    % Assume that p.mesh.data has this field, as it should
    % be initialised by mesh_open, but check anyway
    if ~isfield(p.mesh.data,'Cdata'),
     [p] = mesh_scalp_interp(p);
    end
    % Confirm that p.mesh.data.timeseries{p.mesh.current} has data
    if length(p.mesh.data.Cdata) < p.mesh.current,
     [p] = mesh_scalp_interp(p);
    end
    if isempty(p.mesh.data.Cdata{p.mesh.current}),
     [p] = mesh_scalp_interp(p);
    end
    
    % Now we should have a mesh timeseries
    
    % Verify the timeseries data against the elec timeseries;
    % if generated with mesh_scalp_interp, they should be equal
    TMP = p.mesh.data.Cdata{p.mesh.current};
    Nelec = size(p.volt.data,2);
    if ~isequal(p.volt.data',TMP(1:Nelec,:)),
      % Something has changed, either a new electrode 
      % set or scalp mesh, so run the interpolation
      msg = sprintf('scalp time series not consistent with electrode time series - recalculating.\n');
      warning(msg);
      p.mesh.data.Cdata{p.mesh.current} = [];
     [p] = mesh_scalp_interp(p);
      TMP = p.mesh.data.Cdata{p.mesh.current};
      if ~isequal(p.volt.data',TMP(1:Nelec,:)),
        msg = sprintf('scalp timeseries not consistent with electrode timeseries.\n');
        error(msg);
      end
    end
    clear TMP;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% surface plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if or(p.elec.plotSurf,p.mesh.plotSurf),
  
  if ~isequal(p.volt.sampleTime,0),
    name = sprintf('Surface: %s @ %8.2f msec',p.volt.file, p.volt.sampleTime);
  else
    name = sprintf('Surface: %s',p.volt.file);
  end
  
  if p.elec.plotSurf,
    
    [p.mesh.current,meshExists] = mesh_check(p,'elec');
    
    if ~meshExists,
      p.mesh.data.meshtype{p.mesh.current}    = 'elec';
      p.mesh.data.vertices{p.mesh.current}    = [Xp,Yp,Zp];
      p.mesh.data.faces{p.mesh.current}       = convhulln([Xp,Yp,Zp]);
      p.mesh.data.lapmat{p.mesh.current}      = [];
      p.mesh.data.lapint{p.mesh.current}      = [];
      p.mesh.data.timeseries{p.mesh.current}  = p.volt.timeArray(:,1);
      p.mesh.data.Cdata{p.mesh.current}       = p.volt.data';
    end
    
    %p = eeg_plot_surf(p);
    %return
    
    H.gui = figure('NumberTitle','off','Name',name,...
      'PaperType','A4','PaperUnits','centimeters');
    
    set(gca,'Projection','perspective')
    set(gca,'DataAspectRatio',[1 1 1]);
    
    p.elec.patch = patch('vertices',[Xp,Yp,Zp],'faces',convhulln([Xp,Yp,Zp]),...
      'facevertexCdata',V,'facecolor','interp','EdgeColor','none',...
      'FaceLighting','phong');
    %patch('vertices',p.elec.data.Vsp,'faces',p.elec.data.Fsp,...
    %      'facevertexCdata',p.elec.data.Psp,'facecolor','interp','EdgeColor','none');
    
    set(gca,'Projection','perspective')
    set(gca,'DataAspectRatio',[1 1 1]);
    axis off tight vis3d
    
    lighting phong
    set(H.gui,'Renderer','zbuffer')
    
    % Create a light above the z axis
    H.light = light('style','infinite','position',[0 0 1000]);
    
    % For a near skin color
    %set(gca,'AmbientLightColor',[.9 .8 .7]);
    
    % MATERIAL([ka kd ks n sc]) sets the ambient/diffuse/specular strength,
    %    specular exponent and specular color reflectance of the objects.
    %p.reflect = material('dull');
    p.reflect{1} = 0.9;
    p.reflect{2} = 0.1;
    p.reflect{3} = 0.0;
    p.reflect{4} = 500;
    p.reflect{5} = 0;
    set(p.elec.patch,'AmbientStrength',p.reflect{1});
    set(p.elec.patch,'DiffuseStrength',p.reflect{2});
    set(p.elec.patch,'SpecularStrength',p.reflect{3});
    set(p.elec.patch,'SpecularExponent',p.reflect{4});
    set(p.elec.patch,'SpecularColorReflectance',p.reflect{5});
    
    if p.elec.plot,
      hold on, plot3(Xp,Yp,Zp,'k.');
    end
    
    set(gca,'Visible','off');
    colormap(p.colorMap.map);
    caxis([p.minimumIntensity p.maximumIntensity]);
    colorbar
    
    H.p = p;
    H.p.mesh.plotSurf = 0;
    set(H.gui,'userdata',H);
    
    SaveGraphics(H.gui,'_3Delec_',p);
    if isequal(exist('mouse_rotate'),2),
      mouse_rotate;
    else
      rotate3D;
    end
    if isequal(exist('gui_topo_animate'),2),
      gui_topo_animate('init',p);
    end
  end
  
  
  % Plot a mesh surface
  
  if p.mesh.plotSurf,
    
    [p.mesh.current,meshExists] = mesh_check(p,'scalp');
    
   [p] = eeg_plot_surf(p);
    
    return
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contour in 2-D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CONTOUR(Z,N) and CONTOUR(X,Y,Z,N) draw N contour lines, 
%overriding the automatic value.
%CONTOUR(Z,V) and CONTOUR(X,Y,Z,V) draw LENGTH(V) contour lines 
%at the values specified in vector V. 

if isequal(p.contour.raw2D,1),
  
  if ~isequal(p.volt.sampleTime,0),
    name = sprintf('2-D Contour: %s @ %8.2f msec',p.volt.file, p.volt.sampleTime);
  else
    name = sprintf('2-D Contour: %s',p.volt.file);
  end
  fig = figure('NumberTitle','off','Name',name); colormap(p.colorMap.map);
  set(gca,'Projection','perspective')
  %set(gca,'Projection','orthographic')
  set(gca,'DataAspectRatio',[1 1 1]);
  
  %trisurf(delaunay(Xp,Yp),Xp,Yp,Zp,V,'EdgeColor','none','FaceColor','interp');
  %hold on, plot3(Xp,Yp,Zp,'.');
  %contour(Xi,Yi,Vi, p.contour.Nsteps) %, colorbar, hold on
  
  [c,ch,cf] = contourf(Xi,Yi,Vi,p.contour.levels);
  absmax = max(max(abs(Vi)));caxis([-absmax absmax]);
  %[C,CH,CF] = CONTOURF(...) also returns a column vector CH of handles
  %to PATCH objects and the contour matrix CF for the filled areas. The
  %UserData property of each object contains the height value for each
  %contour.
  
  % FaceColor is the contour patch color
  % EdgeColor is the contour line color
  for i=1:size(ch),
    height = get(ch(i),'UserData');
    
    if     (height > 0),
      set(ch(i),'LineStyle', '-'); % options are: -, --, :, -., none
      set(ch(i),'EdgeColor',[1 1 1]) % white contour lines, black = [0 0 0] (ie, [R G B])
    elseif (height < 0),
      set(ch(i),'LineStyle', ':'); % options are: -, --, :, -., none
    else
      set(ch(i),'LineStyle', '-'); % options are: -, --, :, -., none
      set(ch(i),'LineWidth',2);
    end
  end
  
  set(gca,'Visible','off'); %hold on; DrawHead(Xrad,'black'); hold off;
  colorbar
  
  SaveGraphics(fig,'_2D_',p);
  if isequal(exist('mouse_rotate'),2),
    mouse_rotate;
  else
    rotate3D;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Projected Contour in 2-D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isequal(p.contour.plot2D,1)
  
  % Project 3D electrode locations onto 2D plane.
  [X2p,Y2p] = elec_3d_2d(Xp,Yp,Zp,Zrad);
  m = max(abs([X2p;Y2p]));
  m = ceil(m);
  if (p.grid.method == 2)
    % Use grid resolution method
    X2g = -m:p.grid.res:m;
  else
    % Use grid size method
    X2g = linspace( -m, m, p.grid.size );
  end
  Y2g = X2g';
  [X2i Y2i V2i] = griddata(X2p,Y2p,V,X2g,Y2g,p.interpMethod);
  
  if ~isequal(p.volt.sampleTime,0),
    name = sprintf('Projected 2D Contours: %s @ %8.2f msec',p.volt.file, p.volt.sampleTime);
  else
    name = sprintf('Projected 2D Contours: %s',p.volt.file);
  end
  fig = figure('NumberTitle','off','Name',name);
  set(gca,'Projection','perspective')
  %set(gca,'Projection','orthographic')
  set(gca,'DataAspectRatio',[1 1 1]);
  
  colormap(p.colorMap.map);
  contourf(X2i,Y2i,V2i, p.contour.Nsteps);
  absmax = max(max(abs(V2i)));caxis([-absmax absmax]);colorbar
  xlabel('X'); ylabel('Y');
  
  SaveGraphics(fig,'_proj2D_',p);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Projected Contour in 3-D
% Draw the contour lines,
%   Get the contours in 2-D.
%   For each point in each contour line, know x,y
%   Calculate z for the x,y from ellipse formula
%   Plot3 point x,y,z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (p.contour.plot3D == 1)
  
  
  % Project 3-D electrode locations onto 2-D plane.
  [X2p,Y2p] = elec_3d_2d(Xp,Yp,Zp,Zrad);
  m = max(abs([X2p;Y2p]));
  m = ceil(m);
  if isequal(p.grid.method,2)
    % Use grid resolution method
    X2g = -m:p.grid.res:m;
  else
    % Use grid size method
    X2g = linspace( -m, m, p.grid.size );
  end
  Y2g = X2g';
  
  [X2i Y2i V3i] = griddata(Xp, Yp, V,X2g,Y2g,p.interpMethod);
  
  % Caculate the contour matrix
  C = contourc(X2g, Y2g, V3i, p.contour.Nsteps);
  
  if ~isequal(p.volt.sampleTime,0),
    name = sprintf('3D Contours: %s @ %8.2f msec',p.volt.file, p.volt.sampleTime);
  else
    name = sprintf('3D Contours: %s',p.volt.file);
  end
  
  %fig = figure('NumberTitle','off','Name',name);
  %set(gca,'Projection','perspective')
  %set(gca,'DataAspectRatio',[1 1 1]);
  %[C,Cpatch] = contourf(X2g, Y2g, V3i, p.contour.Nsteps);
  %get(Cpatch(1))
  %P.cdata    = get(Cpatch,'CData');
  %P.vertices = get(Cpatch,'Vertices');
  %P.faces    = get(Cpatch,'Faces');
  %close(fig);
  %fig = figure('NumberTitle','off','Name',name); hold on
  %i = 0;
  %while i < (length(Cpatch) - 1),
  %%vertices = P.vertices{i};
  %i = i + 1;
  %vertices = [P.vertices{i}; P.vertices{i+1}];
  %Xc = vertices(:,1);
  %Yc = vertices(:,2);
  %[X3,Y3,Z3] = elec_2d_3d(Xc,Yc,Xrad,Yrad,Zrad);
  %tri = delaunay(X3,Y3);
  %Htri = trisurf(tri,X3,Y3,Z3,P.cdata{i},'EdgeColor','none','FaceColor','flat');
  
  %%vertices = [X3 Y3 Z3];
  %%Hpatch = patch('Vertices',vertices,'Faces',P.faces{i},'CData',P.cdata{i},'FaceColor','flat');
  %contourValues(i) = P.cdata{i};
  %end
  
  %view(0,90), axis tight, rotate3d
  %set(gca,'Visible','off');
  %colormap(p.colorMap.map);
  %absmax = max(abs(contourValues));
  %caxis([-absmax absmax]);colorbar;
  
  
  
  
  fig = figure('NumberTitle','off','Name',name); hold on
  set(gca,'Projection','perspective')
  %set(gca,'Projection','orthographic')
  set(gca,'DataAspectRatio',[1 1 1]);
  
  % When we view the contour lines, those on the other side of
  % the head can be seen, which makes it confusing.
  
  % The array C is 2-row, x on 1st, y on 2nd.
  % Each contour set has a preceeding column, where
  % 1st is value of contour and 2nd is number of points.
  %
  % We need to get each contour separately, and calculate
  % Z for each X,Y.
  i = 1; j = 0;
  limit = size(C,2);
  while (i < limit),
    
    j = j + 1; contourValues(j) = C(1,i);
    
    numberPoints = C(2,i);
    
    endContourSet = i + numberPoints;
    
    Xc = [ C(1,i+1:endContourSet) ]';
    Yc = [ C(2,i+1:endContourSet) ]';
    
    %cdata = contourValues(j) + 0*Xc;  % Make cdata the same size as xdata
    
    % Given each X,Y we project back from 2D to 3D.
    [X3,Y3,Z3] = elec_2d_3d(Xc,Yc,Xrad,Yrad,Zrad);
    
    if (contourValues(j) < 0)
      % Values < 0 with dashed black line
      if isequal(p.colorMap.style,'Gray')
        line(X3,Y3,Z3,'LineStyle',':','LineWidth',0.5,'Color','black');
      else
        line(X3,Y3,Z3,'LineStyle',':','LineWidth',0.5,'Color','blue');
      end
    elseif (contourValues(j) > 0)
      % Values > 0 with solid black line
      if isequal(p.colorMap.style,'Gray')
        line(X3,Y3,Z3,'LineStyle','-','LineWidth',0.5,'Color','black');
      else
        line(X3,Y3,Z3,'LineStyle','-','LineWidth',0.5,'Color','red');
      end
    else
      % Values = 0 with thicker solid black line
      line(X3,Y3,Z3,'LineStyle','-','LineWidth',0.75,'Color','black');
    end
    i = endContourSet + 1;
  end
  
  view(0,90), axis tight, rotate3d
  set(gca,'Visible','off');
  %colormap(p.colorMap.map); H = colorbar;
  %absmax = max(abs(contourValues));
  
  SaveGraphics(fig,'_3D_',p);
end

return    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DrawHead(Hradius, Hcolor)

switch Hcolor
  case 'black'
    HCOLOR = [0 0 0];
  case 'white'
    HCOLOR = [1 1 1];
  otherwise
    HCOLOR = [0 0 0];
end
HLINEWIDTH = 30;

rmax = 0.95 * Hradius;

% Plot Head

theta1 = linspace(0,2*pi,50);       % 360 degree rotation
thetad = (theta1(2) - theta1(1))/2; % provide staggered rotation
theta = [theta1, theta1 + thetad];

plot(cos(theta).*rmax,sin(theta).*rmax,'color',HCOLOR,'LineWidth',HLINEWIDTH);

% Plot Nose
width = rmax * 0.1;
tip   = rmax * 1.075;
base  = rmax * 1.005;

plot([width;0;-width],[rmax;tip;rmax],'Color',HCOLOR,'LineWidth',HLINEWIDTH);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveGraphics(F,type,p)
if ~isequal(p.saveGraphics,'No Save Plots'),
  
  [path,file,ext] = fileparts(strcat(p.volt.path, filesep, p.volt.file));
  file = strcat(file, type, num2str(p.volt.samplePoint),'.',p.saveGraphics);
  file = fullfile(path,file);
  saveas(F,file);
end
return
