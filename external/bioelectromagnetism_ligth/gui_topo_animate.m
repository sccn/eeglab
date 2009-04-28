function [p] = gui_topo_animate(command,p)

% gui_topo_animate - GUI controls for animating topography
%
% Usage: [p] = gui_topo_animate(command,[p])
%
% command   - either 'init' or 'animate'
%
% p is the eeg_toolbox struct
%
% This function adds controls to the current figure.
% It assumes this figure has been created with a topographic
% map patch and the userdata contains the p struct for the
% data of the patch.
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:56 $

% Licence:  GNU GPL, no express or implied warranties
% History:  05/2002, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('command','var'), command = 'init'; end

switch command,
  
  case 'init',
    
    if exist('p','var'),
      if ~isempty(p),
        ANIM = INIT(p);
      end
    else
      ANIM = INIT;
    end
    
  case 'animate',
    
    ANIM = get(gcbf,'Userdata');
    if isempty(ANIM),
      ANIM = get(gcf,'Userdata');
    end
    
    % Obtain the current view from mouse_rotate, which
    % can be used in eeg_save_graphics
    views = get(ANIM.View,'String');
    ANIM.p.topoView = lower(char(views(get(ANIM.View,'Value'))));
    
    
    
    % ---  First get patch handle (Hp)  ---
    
    axsibs = get(ANIM.axis,'Children');
    for i=1:length(axsibs),
      type = get(axsibs(i),'Type');
      if isequal(type,'patch'),
        ANIM.Hp = axsibs(i);
      end
    end
    
    
    % Animation setup
    ANIM = toggle_visible(ANIM);
    set(ANIM.Hp,'EraseMode','normal');
    figure(ANIM.gui);
    
    
    % ---  Get the animation parameters (msec)  ---
    
    start_msec  = str2num(get(ANIM.Start,'string'));
    finish_msec = str2num(get(ANIM.Finish,'string'));
    step_msec   = str2num(get(ANIM.Step,'string'));
    
    
    % --- get data from p struct ---
    
    [ sample, timeArray ] = extract_data(ANIM);
    
    
    % ---  Check step parameter  ---
    
    step = step_msec/sample.rate;
    stepREM = step - fix(step);
    if stepREM,
      nointerp = 0;
      % have not yet implemented interpolation
      warning('GUI_TOPO_ANIMATE: Cannot interpolate for steps ~= sample rate\n');
      step = 1;
      set(ANIM.Step,'string',sprintf('%7.2f',sample.rate));
    else
      nointerp = 1;
    end
    
    
    
    % ---  Convert start/finish to data rows  ---
    
    [ index, point ] = NearestXYArrayPoint( timeArray, start_msec, 'exact' );
    if isempty(index),
      [ index, point ] = NearestXYArrayPoint( timeArray, start_msec );
    end
    start_index = index;
    start_point = point;
    
    [ index, point ] = NearestXYArrayPoint( timeArray, finish_msec, 'exact' );
    if isempty(index),
      [ index, point ] = NearestXYArrayPoint( timeArray, finish_msec );
    end
    finish_index = index;
    finish_point = point;
    
    set(ANIM.Start,'value',start_point);
    set(ANIM.Start,'string',sprintf('%7.2f',start_point));
    set(ANIM.Finish,'value',finish_point);
    set(ANIM.Finish,'string',sprintf('%7.2f',finish_point));
    
    
    % ---  Verify start/finish values  ---
    
    if isequal(start_index,finish_index),
      
      % Given start == finish, update topo for that time
      
      % set voltage scale limits
      ANIM = voltage_scale(ANIM,start_index,finish_index);
      caxis([ANIM.p.minimumIntensity ANIM.p.maximumIntensity]);
      colorbar
      
      t = start_index;
      
      ANIM = update_topo(ANIM,t,step);
      
      ANIM = toggle_visible(ANIM);
     [p] = ANIM.p;
      return;
      
    elseif start_index > finish_index,
      ANIM = toggle_visible(ANIM);
      error('GUI_TOPO_ANIMATE: Finish < Start!\n');
    end
    
    
    
    
    
    
    % ---  Run the animation  ---
    
    % Update voltage scale limits
    ANIM = voltage_scale(ANIM,start_index,finish_index);
    caxis([ANIM.p.minimumIntensity ANIM.p.maximumIntensity]);
    colorbar
    
    % do we want a movie?
    doMovie = ANIM.Movie;
    if doMovie,
      movieFrameIndex = 0;
      % Define movie region as figure size
      rect = get(ANIM.gui,'Position');
      rect(1:2) = [0 0];
    end
    
    for t=start_index : step : finish_index,
      
      ANIM = update_topo(ANIM,t,step);
      
      % Capture movie frame
      if doMovie,
        movieFrameIndex = movieFrameIndex + 1;
        movieMatrix(:,movieFrameIndex) = getframe(ANIM.gui,rect);
      end
    end
    ANIM = toggle_visible(ANIM);
    
    
    if doMovie,
      
      % Play movie in another figure
      playtimes = 10;
      framesPerSecond = 10;
      
      movieFigure = figure('NumberTitle','off','Name','Movie Playback (10 FPS, 10 LOOPS)',...
        'PaperType','A4','PaperUnits','centimeters');
      
      axis off;
      
      movie(movieFigure,movieMatrix,playtimes,framesPerSecond,rect);
      
      selection = 'replay';
      while isequal(selection,'replay'),
        
        [selectedButton,dlgShown]=uigetpref('eeg_toolbox','saveTopoMovie',...
          'Movie Controls',...
          'Replay, save or close this movie?',...
          {'replay','save','close';'Replay','Save','Close'},...
          'CheckboxString','save this preference',...
          'DefaultButton','Close');
        
        % Repeat, Save or Close movie
        if dlgShown,
          selection = selectedButton;
        else
          selection = getpref('eeg_toolbox','saveTopoMovie');
        end
        
        switch selection,
          case 'replay',
            movie(movieFigure,movieMatrix,playtimes,framesPerSecond,rect);
          case 'save',
            % save to AVI file
            aviFileName = uiputfile('*.avi','Save Movie As .avi');
            if aviFileName,
              movie2avi(movieMatrix,aviFileName,'quality',100);
              filename = strcat(pwd,filesep,aviFileName);
              fprintf('Saved movie to %s\n',filename);
            end
            close(movieFigure);
          otherwise,
            close(movieFigure);
        end
      end
      
    end
    
   [p] = ANIM.p;
    
    
  otherwise,
    
end

return






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ sample, timeArray, Cdata ] = extract_data(ANIM),

% Unlike most of the eeg_toolbox, this function requires
% the potential data to be electrodes in rows, voltage in columns;
% because of the way the matlab patch command uses Cdata.
%
% This is especially so for using this data as Cdata in the patch
% function and facilitates storage and handling of very large
% patch Cdata, such as cortical patches

if ANIM.p.mesh.current & ~isempty(ANIM.p.mesh.data),
  
  if ANIM.p.mesh.current > length(ANIM.p.mesh.data.meshtype),
    % assume we are plotting the p.volt.data with the electrodes
    sample.rate  = ANIM.p.volt.sampleMsec;
    sample.point = ANIM.p.volt.samplePoint;
    sample.time  = ANIM.p.volt.sampleTime;
    timeArray    = ANIM.p.volt.timeArray(:,1);
    Cdata        = ANIM.p.volt.data';
    return
  end
  
  meshtype = lower(ANIM.p.mesh.data.meshtype{ANIM.p.mesh.current});
  switch meshtype,
    case 'elec',
      % assume we are plotting the p.volt.data with the electrodes
      sample.rate  = ANIM.p.volt.sampleMsec;
      sample.point = ANIM.p.volt.samplePoint;
      sample.time  = ANIM.p.volt.sampleTime;
      timeArray    = ANIM.p.volt.timeArray(:,1);
      Cdata        = ANIM.p.volt.data';
    case 'scalp',
      % assume we are plotting the p.volt.data on the scalp
      sample.rate  = ANIM.p.volt.sampleMsec;
      sample.point = ANIM.p.volt.samplePoint;
      sample.time  = ANIM.p.volt.sampleTime;
      timeArray    = ANIM.p.volt.timeArray(:,1);
      % check for vertices in rows (default) or columns
      nvert = size(ANIM.p.mesh.data.vertices{ANIM.p.mesh.current},1);
      [s1,s2] = size(ANIM.p.mesh.data.Cdata{ANIM.p.mesh.current});
      if isequal(nvert,s1), % vertices in rows, timepoints in columns
        Cdata = ANIM.p.mesh.data.Cdata{ANIM.p.mesh.current};
      else,                 % vertices in columns, timeseries in rows
        % try to rotate it, but it may exceed memory for detailed surfaces
        Cdata = ANIM.p.mesh.data.Cdata{ANIM.p.mesh.current}';
      end
    otherwise
      % assume the data is in p.mesh.data.Cdata, with a 
      % timeseries in p.mesh.data.timeseries
      sample.rate  = ANIM.p.mesh.data.timeseries{ANIM.p.mesh.current}(2) - ANIM.p.mesh.data.timeseries{ANIM.p.mesh.current}(1);
      sample.point = ANIM.p.mesh.samplePoint;
      sample.time  = ANIM.p.mesh.sampleTime;
      timeArray    = ANIM.p.mesh.data.timeseries{ANIM.p.mesh.current};
      % get number of vertices
      nvert = size(ANIM.p.mesh.data.vertices{ANIM.p.mesh.current},1);
      % Assume more vertices than time points
      [s1,s2] = size(ANIM.p.mesh.data.Cdata{ANIM.p.mesh.current});
      if isequal(nvert,s1), % vertices in rows, timepoints in columns
        Cdata = ANIM.p.mesh.data.Cdata{ANIM.p.mesh.current};
      else,                 % vertices in columns, timeseries in rows
        % try to rotate it, but it may exceed memory
        Cdata = ANIM.p.mesh.data.Cdata{ANIM.p.mesh.current}';
      end
  end
else
  % assume we are plotting the p.volt.data with the electrodes
  sample.rate  = ANIM.p.volt.sampleMsec;
  sample.point = ANIM.p.volt.samplePoint;
  sample.time  = ANIM.p.volt.sampleTime;
  timeArray    = ANIM.p.volt.timeArray(:,1);
  Cdata        = ANIM.p.volt.data';
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ ANIM ] = update_topo(ANIM,t,step),

[ sample, timeArray, Cdata ] = extract_data(ANIM);


% Do we need to average?
if step > 1,
  Cdata = Cdata(:,t-(step-1):t);
  V = mean(Cdata,2);
else
  V = Cdata(:,t);
end
clear Tmp;

% Use absolute values for cortical surfaces
switch ANIM.p.mesh.data.meshtype{ANIM.p.mesh.current},
  case {'scalp','elec'},
  otherwise
    V = abs(V);
end

set(ANIM.Hp,'FaceVertexCdata',V);

%ANIM.p.timePoint = t;
switch ANIM.p.mesh.data.meshtype{ANIM.p.mesh.current},
  case {'scalp','elec'},
    ANIM.p.volt.samplePoint = t;
  otherwise
    ANIM.p.mesh.samplePoint = t;
end

if (ANIM.p.contour.plot3D == 1),
  ANIM.p = eeg_plot_surf_contours(ANIM.p);
end

drawnow;

% Update the msec text
ANIM.p.volt.sampleTime = timeArray(t);
timeString = sprintf('%8.2f msec',ANIM.p.volt.sampleTime);
set(ANIM.msec,'String',timeString);

% Update the figure name (title)
ANIM.p.volt.sampleTime = timeArray(t);
name = sprintf('Surface: %s @ %8.2f msec',ANIM.p.volt.file, ANIM.p.volt.sampleTime);
set(ANIM.gui,'Name',name);

% Save the new figure to a graphics file
if get(ANIM.Save,'Value'),
  types = get(ANIM.Type,'String');
  type = char(types(get(ANIM.Type,'Value')));
  eeg_save_graphics(ANIM.gui,type,ANIM.p,0);
end

set(ANIM.gui,'Userdata',ANIM);

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ ANIM ] = voltage_scale(ANIM,start,finish),

[ sample, timeArray, Cdata ] = extract_data(ANIM);

yset = 0;
if isfield(ANIM,'Yset'),
  if get(ANIM.Yset,'value'),
    ANIM.p.minimumIntensity = get(ANIM.Ymin,'value');
    ANIM.p.maximumIntensity = get(ANIM.Ymax,'value');
    yset = 1;
  end
end
if ~yset,
  
  % For very large Cdata matrices, it is more memory
  % efficient to use max(Cdata) that first use abs(Cdata)
  Cmax = max(max(Cdata(:,start:finish)));
  Cmin = min(min(Cdata(:,start:finish)));
  absmax = max([Cmax abs(Cmin)]);
  
  ANIM.p.minimumIntensity = -absmax;
  ANIM.p.maximumIntensity =  absmax;
  if isfield(ANIM,'Ymin'),
    set(ANIM.Ymin,'string',sprintf('%7.2f',-absmax));
    set(ANIM.Ymin,'value',-absmax);
  end
  if isfield(ANIM,'Ymax'),
    set(ANIM.Ymax,'string',sprintf('%7.2f', absmax));
    set(ANIM.Ymax,'value', absmax);
  end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ H ] = toggle_visible(H),

toggle = get(H.Ymin,'visible');
if isequal(toggle,'on'),
  set(H.gui,'BackingStore','off','Pointer','watch');
  set(H.Yset,'visible','off');
  set(H.Ymin,'visible','off');
  set(H.Ymax,'visible','off');
  set(H.Start,'visible','off');
  set(H.Finish,'visible','off');
  set(H.Step,'visible','off');
  set(H.BMovie,'visible','off');
  set(H.Animate,'visible','off');
  set(H.Save,'visible','off');
  set(H.Type,'visible','off');
  set(H.AClose,'visible','off');
  % In case mouse rotate is on
  if ishandle(H.RClose),
    set(H.RClose,'visible','off');
    set(H.View,'visible','off');
    set(H.Info,'visible','off');
    set(H.Az,'visible','off');
    set(H.El,'visible','off');
  end
  
else
  set(H.gui,'BackingStore','on','Pointer','arrow');
  set(H.Yset,'visible','on');
  set(H.Ymin,'visible','on');
  set(H.Ymax,'visible','on');
  set(H.Start,'visible','on');
  set(H.Finish,'visible','on');
  set(H.Step,'visible','on');
  set(H.Animate,'visible','on');
  set(H.Save,'visible','on');
  set(H.Type,'visible','on');
  set(H.BMovie,'visible','on');
  set(H.AClose,'visible','on');
  % In case mouse rotate is on
  if ishandle(H.RClose),
    set(H.RClose,'visible','on');
    set(H.View,'visible','on');
    set(H.Info,'visible','on');
    set(H.Az,'visible','on');
    set(H.El,'visible','on');
  end
  
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ index, point ] = NearestXYArrayPoint( data_array, point, type )

if ~exist('type','var') type = ''; end

% In this function, input data_array is an array, not a matrix.
% This function returns the data point in the array
% that has the closest value to the value given (point).  In
% the context of 'gui_topo_animate' the point is a time.

if      point >= max(data_array)
  point  = max(data_array);
  index  = find(data_array == point);
  return;
elseif  point <= min(data_array)
  point  = min(data_array);
  index  = find(data_array == point);
  return;
end

data_sorted = sort(data_array);

greater = find(data_sorted > point);
greater_index = greater(1);

lesser = find(data_sorted < point);
lesser_index = lesser(length(lesser));

greater_dif = data_sorted(greater_index) - point;
lesser_dif  = point - data_sorted(lesser_index);

if     strcmp(type,'exact'),  index = find(data_array == point);
elseif strcmp(type,'nextx'),  index = greater_index;
elseif strcmp(type,'prevx'),  index = lesser_index;
else
  if (greater_dif < lesser_dif)
    index = find(data_array == data_sorted(greater_index));
  else
    index = find(data_array == data_sorted(lesser_index));
  end
end
point = data_array(index);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H] = INIT(p),

H = get(gcbf,'userdata');

% only one per figure
if isfield(H,'Start'),
  if ~isempty(H.Start),
    return;
  end
end

H = get(gcf,'userdata');

% enable right click context menu for animation
if isempty(get(H.gui,'uicontextmenu')),
  menu = uicontextmenu;
  uimenu(menu,'Label','Animate','Callback','gui_topo_animate; ');
  set(H.gui,'uicontextmenu',menu);
else
  menu = get(H.gui,'uicontextmenu');
  sibs = get(menu,'children');
  nomenu = 1;
  for i=1:length(sibs),
    label = get(sibs(i),'Label');
    if strmatch(label,'Animate'),
      nomenu = 0;
    end
  end
  if nomenu,
    uimenu(menu,'Label','Animate','Callback','gui_topo_animate; ');
    set(H.gui,'uicontextmenu',menu);
  end
end
try,
set(H.axis,'uicontextmenu',menu);
axsibs = get(H.axis,'Children');
for i=1:length(axsibs),
  type = get(axsibs(i),'Type');
  if isequal(type,'patch'),
    if isempty(get(axsibs(i),'uicontextmenu')),
      set(axsibs(i),'uicontextmenu',menu);
    end
  end
end
catch, end;
% Match background figure colour
bgcolor = get(H.gui,'Color');
% Try to adapt the foreground colour a little
black = find(bgcolor <= .6);
fgcolor = [0 0 0]; %black text
if length(black)>2, fgcolor = [1 1 1]; end

Font.FontName   = 'Helvetica';
Font.FontUnits  = 'Pixels';
Font.FontSize   = 12;
Font.FontWeight = 'normal';
Font.FontAngle  = 'normal';

if exist('p','var'), H.p = p; end

% -- Set Voltage Range
switch H.p.mesh.data.meshtype{H.p.mesh.current},
  case {'scalp','elec'},
    samplePoint = H.p.volt.samplePoint;
  otherwise
    samplePoint = H.p.mesh.samplePoint;
end
H = voltage_scale(H,samplePoint,samplePoint);

H.msec = uicontrol('Parent',H.gui,'Style','text','Units','Normalized',Font,...
  'Position',[.40 .00 .2 .05],...
  'HorizontalAlignment','left',...
  'BackGroundColor',bgcolor,...
  'ForeGroundColor',fgcolor,...
  'String','msec',...
  'TooltipString','Time of topographic map.');

H.Yset = uicontrol('Parent',H.gui,'Style','checkbox','Units','Normalized',Font,...
  'Position',[.0 .90 .1 .05],...
  'HorizontalAlignment','left',...
  'BackGroundColor',bgcolor,...
  'ForeGroundColor',fgcolor,...
  'String','Y set',...
  'TooltipString','Set specific Y min/max, otherwise automatic.');

H.Ymin = uicontrol('Parent',H.gui,'Style','edit','Units','Normalized',Font,...
  'Position',[.0 .85 .10 .05],...
  'HorizontalAlign','right',...
  'Tag','YMIN',...
  'TooltipString','Set Y min', ...
  'String',sprintf('%7.1f',H.p.minimumIntensity),...
  'Callback',strcat('H = get(gcbf,''userdata'');',...
  'ymin = str2num(get(H.Ymin,''string''));',...
  'H.p.minimumIntensity = ymin;',...
  'set(H.Ymin,''string'',sprintf(''%7.2f'',ymin));',...
  'set(H.Ymin,''value'',ymin);',...
  'set(gcbf,''userdata'',H); clear H ymin;'));

H.Ymax = uicontrol('Parent',H.gui,'Style','edit','Units','Normalized',Font,...
  'Position',[.0 .80 .10 .05],...
  'HorizontalAlign','right',...
  'Tag','YMAX',...
  'TooltipString','Set Y max', ...
  'String',sprintf('%7.1f',H.p.maximumIntensity),...
  'Callback',strcat('H = get(gcbf,''userdata'');',...
  'ymax = str2num(get(H.Ymax,''string''));',...
  'H.p.maximumIntensity = ymax;',...
  'set(H.Ymax,''string'',sprintf(''%7.2f'',ymax));',...
  'set(H.Ymax,''value'',ymax);',...
  'set(gcbf,''userdata'',H); clear H ymax;'));

% -- Set Time Range

H.Start = uicontrol('Parent',H.gui,'Style','edit','Units','Normalized',Font,...
  'Position',[.0 .70 .1 .05],...
  'HorizontalAlignment','right',...
  'BackGroundColor',bgcolor,...
  'ForeGroundColor',fgcolor,...
  'String','','TooltipString','Start Time (msec)',...
  'Callback',strcat('H = get(gcbf,''userdata''); ',...
  'val = str2num(get(H.Start,''string'')); ',...
  'set(H.Start,''value'',val); ',...
  'set(H.Start,''string'',sprintf(''%7.2f'',val)); ',...
  'clear H val;'));

H.Finish  = uicontrol('Parent',H.gui,'Style','edit','Units','Normalized',Font,...
  'Position',[.0 .65 .1 .05],...
  'HorizontalAlignment','right',...
  'BackGroundColor',bgcolor,...
  'ForeGroundColor',fgcolor,...
  'String','','TooltipString','Finish Time (msec)',...
  'Callback',strcat('H = get(gcbf,''userdata''); ',...
  'val = str2num(get(H.Finish,''string'')); ',...
  'set(H.Finish,''value'',val); ',...
  'set(H.Finish,''string'',sprintf(''%7.2f'',val)); ',...
  'clear H val;'));

H.Step    = uicontrol('Parent',H.gui,'Style','edit','Units','Normalized',Font,...
  'Position',[.0 .60 .1 .05],...
  'HorizontalAlignment','right',...
  'BackGroundColor',bgcolor,...
  'ForeGroundColor',fgcolor,...
  'String','','TooltipString','Time Step (msec)',...
  'Callback',strcat('H = get(gcbf,''userdata''); ',...
  'val = str2num(get(H.Step,''string'')); ',...
  'set(H.Step,''value'',val); ',...
  'set(H.Step,''string'',sprintf(''%7.2f'',val)); ',...
  'clear H val;'));


[ sample, timeArray, Cdata ] = extract_data(H); % local function


set(H.msec,  'String',sprintf('%7.2f msec',sample.time));
set(H.Start, 'String',sprintf('%7.2f',sample.time));
set(H.Finish,'String',sprintf('%7.2f',sample.time));
set(H.Step,  'String',sprintf('%7.2f',sample.rate));


% -- Set Image and Movie Save Options

H.Save = uicontrol('Parent',H.gui,'Style','checkbox','Units','Normalized',Font,...
  'Position',[.0 .55 .1 .05],...
  'HorizontalAlignment','left',...
  'BackGroundColor',bgcolor,...
  'ForeGroundColor',fgcolor,...
  'String','Save',...
  'TooltipString','Save each step to graphics file (via eeg_save_graphics).');

H.Type = uicontrol('Parent',H.gui,'Style','popup','Units','Normalized',Font,...
  'Position',[.0 .50 .1 .05],...
  'HorizontalAlignment','left',...
  'BackGroundColor',bgcolor,...
  'ForeGroundColor',fgcolor,...
  'String',{'png','jpeg','tiff','eps'},...
  'Value',1,...
  'TooltipString','Save format (warning: TIFF and EPS are large files).');

% -- Set Action Controls

H.Animate = uicontrol('Parent',H.gui,'Style','pushbutton','Units','Normalized',Font,...
  'Position',[.0 .40 .1 .045],...
  'HorizontalAlignment','center',...
  'BackGroundColor',[0 .6 0],...
  'ForeGroundColor',[1 1 1],...
  'String','Animate',...
  'Callback','gui_topo_animate(''animate'');');

H.Movie = 0;
H.BMovie = uicontrol('Parent',H.gui,'Style','pushbutton','Units','Normalized',Font,...
  'Position',[.0 .35 .1 .045],...
  'HorizontalAlignment','center',...
  'BackGroundColor',[0 .6 0],...
  'ForeGroundColor',[1 1 1],...
  'String','Movie',...
  'Callback',strcat('H = get(gcbf,''userdata'');',...
  'H.Movie = 1; ',...
  'set(gcbf,''userdata'',H); clear H;',...
  'gui_topo_animate(''animate'');'));

H.AClose = uicontrol('Parent',H.gui,'Style','pushbutton','Units','Normalized',Font,...
  'Position',[.0 .30 .1 .045],...
  'HorizontalAlignment','center',...
  'BackGroundColor',[.6 0 0],...
  'ForeGroundColor',[1 1 1],...
  'String','Close','Value',0,...
  'TooltipString','Use right click context menu after ''Close'' to get animation back.',...
  'Callback',strcat('H = get(gcbf,''userdata'');',...
  'delete(H.Yset);    delete(H.Ymin);   delete(H.Ymax); ',...
  'delete(H.msec);    delete(H.Start);  delete(H.Finish); delete(H.Step); ',...
  'delete(H.Animate); delete(H.Save);   delete(H.Type); ',...
  'delete(H.BMovie);  delete(H.AClose); ',...
  'fields = fieldnames(H); ',...
  'trashfields(1) = strmatch(''msec'',  fields,''exact''); ',...
  'trashfields(1) = strmatch(''Start'', fields,''exact''); ',...
  'trashfields(2) = strmatch(''Finish'',fields,''exact''); ',...
  'trashfields(3) = strmatch(''Step'',  fields,''exact''); ',...
  'trashfields(4) = strmatch(''Animate'', fields,''exact''); ',...
  'trashfields(5) = strmatch(''Save'',  fields,''exact''); ',...
  'trashfields(6) = strmatch(''Type'',  fields,''exact''); ',...
  'trashfields(7) = strmatch(''AClose'',fields,''exact''); ',...
  'trashfields(8) = strmatch(''Ymin'',  fields,''exact''); ',...
  'trashfields(9) = strmatch(''Ymax'',  fields,''exact''); ',...
  'trashfields(10) = strmatch(''Yset'', fields,''exact''); ',...
  'trashfields(11) = strmatch(''BMovie'',fields,''exact''); ',...
  'for i=1:12, ',...
  '   H = setfield(H,char(fields(trashfields(i))),[]); ',...
  'end; ',...
  'set(H.gui,''userdata'',H); clear H trashfields fields; '));


set(H.gui,'userdata',H);

return
