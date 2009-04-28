function [cnt] = eeg_view_cnt(filename,command,parent)

% eeg_view_cnt - Plot scan CNT data
% 
% eeg_view_cnt(filename,command,parent)
% 
%   filename    -   a string filename
%   command     -   a string, 'init', 'hslider', 'vslider'
%   parent      -   optional, for GUI parent
% 
% Note: Works only with Scan 4.1+ data files
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:53 $

% Licence:  GNU GPL, no express or implied warranties
% History:  2002, Sean.Fitzgibbon@flinders.edu.au
%           06/2002, Darren.Weber_at_radiology.ucsf.edu
%                    adapted to eeg_toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('filename','var'),
  msg = sprintf('EEG_VIEW_CNT: No input filename\n');
  error(msg);
end
if ~exist('command','var'), command = 'init'; end


switch command,
  
  case 'init',
    
    if exist('parent','var'),
      VIEWCNT = INIT(filename,parent);
    else
      VIEWCNT = INIT(filename);
    end
    
  case 'plot',
    
    VIEWCNT = get(gcbf,'userdata');
    
    % ------- Baseline parameters -------
    
    base = str2num(get(VIEWCNT.handles.baseline,'String'));
    
    % adjust baseline to multiple of sample rate
    minbase = (1/VIEWCNT.data.cnt.srate)*1000;
    base = floor(base / minbase) * minbase;
    % make sure baseline is within range of display
    if (base < minbase), base = minbase; end;
    maxbase = (VIEWCNT.data.dispNPoints/VIEWCNT.data.cnt.srate)*1000;
    if (base > maxbase), base = maxbase; end;
    % update baseline number displayed
    set(VIEWCNT.handles.baseline,'String',sprintf('%6.2f',base));
    VIEWCNT.data.base = 1:base;
    
    %  ------- Seconds displayed per page -------
    
    VIEWCNT.data.dispSeconds = str2num(get(VIEWCNT.handles.seconds,'String'));
    
    % must display at least 2 time points
    if (VIEWCNT.data.dispSeconds <= 0),
      minsec = (1000/VIEWCNT.data.cnt.srate) / 1000;
      VIEWCNT.data.dispSeconds = 2 * minsec;
    end;
    
    if (VIEWCNT.data.dispSeconds > VIEWCNT.data.points / VIEWCNT.data.cnt.srate ),
      % Request to display more seconds than currently loaded
      totalseconds = VIEWCNT.data.cnt.numSamples / VIEWCNT.data.cnt.srate;
      if(VIEWCNT.data.dispSeconds > totalseconds ),
        % Load ALL data
        VIEWCNT.data.dispSeconds = totalseconds;
        VIEWCNT.data.dispNPoints = floor(VIEWCNT.data.dispSeconds * VIEWCNT.data.cnt.srate);
        
        range = 'all';
      else,
        % Load more data
        VIEWCNT.data.dispNPoints = floor(VIEWCNT.data.dispSeconds * VIEWCNT.data.cnt.srate);
        VIEWCNT.data.dispEPoint = VIEWCNT.data.dispSPoint - 1 + VIEWCNT.data.dispNPoints;
        
        range(1) = VIEWCNT.data.dispSPoint;
        range(2) = VIEWCNT.data.dispEPoint;
      end;
      filename = [VIEWCNT.data.cnt.path,VIEWCNT.data.cnt.name];
      VIEWCNT.data.cnt = eeg_load_scan4_cnt(filename,'all',range);
      VIEWCNT.data.points = size(VIEWCNT.data.cnt.volt,1);
    else
      VIEWCNT.data.dispNPoints = floor(VIEWCNT.data.dispSeconds * VIEWCNT.data.cnt.srate);
      VIEWCNT.data.dispEPoint = VIEWCNT.data.dispSPoint - 1 + VIEWCNT.data.dispNPoints;
    end
    
    set(VIEWCNT.handles.seconds,'String',sprintf('%8.4f',VIEWCNT.data.dispSeconds));
    
    
    
    
    %  ------- Slider parameters -------
    
    slider = get(VIEWCNT.handles.Slider,'Value');
    
    sec = slider/VIEWCNT.data.cnt.srate;
    set(VIEWCNT.handles.Tslider,'String',sprintf('%8.4f sec',sec));
    
    point = floor(slider);
    if (point > VIEWCNT.data.dispEPoint),
      % Slider request for data beyond display range
      pointdif = point - VIEWCNT.data.dispEPoint;
      range(1) = VIEWCNT.data.dispSPoint -1 + pointdif;
      range(2) = range(1) + VIEWCNT.data.dispNPoints;
      
      if (range(2) > VIEWCNT.data.cnt.numSamples),
        % Ensure display range within data range
        range(1) = VIEWCNT.data.cnt.numSamples - VIEWCNT.data.dispNPoints;
        range(2) = VIEWCNT.data.cnt.numSamples;
      end;
      
      % Update display range and check whether or not to reload data
      if (VIEWCNT.data.dispSPoint < range(1)),
        filename = [VIEWCNT.data.cnt.path,VIEWCNT.data.cnt.name];
        VIEWCNT.data.cnt = eeg_load_scan4_cnt(filename,'all',range);
        VIEWCNT.data.dispSPoint = range(1);
      elseif (VIEWCNT.data.dispEPoint < range(2)),
        filename = [VIEWCNT.data.cnt.path,VIEWCNT.data.cnt.name];
        VIEWCNT.data.cnt = eeg_load_scan4_cnt(filename,'all',range);
        VIEWCNT.data.dispEPoint = range(2);
      end
      VIEWCNT.data.dispSPoint = range(1);
      VIEWCNT.data.dispEPoint = range(2);
      
    elseif (point < VIEWCNT.data.dispSPoint),
      % Slider request for data beyond display range
      pointdif = VIEWCNT.data.dispSPoint - point;
      range(1) = VIEWCNT.data.dispSPoint -1 - pointdif;
      if range(1) < 1, range(1) = 1; end;
      range(2) = range(1) + VIEWCNT.data.dispNPoints;
      
      if (range(2) > VIEWCNT.data.cnt.numSamples),
        % Ensure display range within data range
        range(1) = VIEWCNT.data.cnt.numSamples - VIEWCNT.data.dispNPoints;
        range(2) = VIEWCNT.data.cnt.numSamples;
      end;
      
      % Update display range and check whether or not to reload data
      if (VIEWCNT.data.dispSPoint > range(1)),
        filename = [VIEWCNT.data.cnt.path,VIEWCNT.data.cnt.name];
        VIEWCNT.data.cnt = eeg_load_scan4_cnt(filename,'all',range);
        VIEWCNT.data.dispSPoint = range(1);
      elseif (VIEWCNT.data.dispEPoint > range(2)),
        filename = [VIEWCNT.data.cnt.path,VIEWCNT.data.cnt.name];
        VIEWCNT.data.cnt = eeg_load_scan4_cnt(filename,'all',range);
        VIEWCNT.data.dispEPoint = range(2);
      end
      VIEWCNT.data.dispSPoint = range(1);
      VIEWCNT.data.dispEPoint = range(2);
    end;
    
    
    stime = (VIEWCNT.data.dispSPoint /VIEWCNT.data.cnt.srate) * 1000;
    etime = (VIEWCNT.data.dispEPoint /VIEWCNT.data.cnt.srate) * 1000;
    
    VIEWCNT.data.time = [stime:VIEWCNT.data.tstep:etime]';
    VIEWCNT.data.time = VIEWCNT.data.time(1:VIEWCNT.data.dispNPoints);
    
    volt = VIEWCNT.data.cnt.volt(1:VIEWCNT.data.dispNPoints,:);
    
    
    
    VIEWCNT.data.basemean = mean(VIEWCNT.data.cnt.volt(VIEWCNT.data.base,:));
    basemean = repmat(VIEWCNT.data.basemean,VIEWCNT.data.dispNPoints,1);
    
    % Now calculate voltage data and plot it
    volt = volt - basemean;
    plot(VIEWCNT.data.time,volt);
    set(VIEWCNT.axes,'XLim',[VIEWCNT.data.time(1) VIEWCNT.data.time(end)]);
    eeg_plot_metric;
    set(VIEWCNT.gui,'userdata',VIEWCNT);
    
  otherwise
    
    close gcbf;
  end
  
  cnt = VIEWCNT.data.cnt;
  
  return
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [VIEWCNT] = INIT(filename,parent),
  
  % GUI General Parameters
  
  GUIwidth  = 150;
  GUIheight =  40;
  
  GUI = figure('Name','CNT View [alpha 1.0]','Tag','CNTVIEW',...
    'NumberTitle','off',...
    'units','characters',...
    'MenuBar','none','Position',[1 1 GUIwidth GUIheight]);
  movegui(GUI,'center');
  
  VIEWCNT.gui = GUI;
  
  VIEWCNT.axes = axes('Parent',GUI,'YDir','reverse');
  
  Font.FontName   = 'Helvetica';
  Font.FontUnits  = 'Pixels';
  Font.FontSize   = 12;
  Font.FontWeight = 'normal';
  Font.FontAngle  = 'normal';
  
  % ---- Display Parameters
  
  data.dispSeconds =  1; % seconds per screen
  data.chanPage =  4; % channels displayed per page (not used yet)
  
  % ---- Control Parameters
  
  % Load first 2 points of CNT file to check parameters
  data.cnt = eeg_load_scan4_cnt(filename,'all',[1 2]);
  
  data.dispNPoints = data.dispSeconds * data.cnt.srate;
  data.dispSPoint = 1;
  data.dispEPoint = data.dispNPoints;
  
  range(1) = data.dispSPoint;
  range(2) = data.dispEPoint;
  data.cnt = eeg_load_scan4_cnt(filename,'all',range);
  
  data.points = size(data.cnt.volt,1);
  data.nChan  = size(data.cnt.volt,2);
  
  data.totalsec  = data.cnt.numSamples / data.cnt.srate;
  data.totalmsec = data.totalsec * 1000;
  
  stime = (data.dispSPoint / data.cnt.srate) * 1000;
  etime = (data.dispEPoint / data.cnt.srate) * 1000;
  
  data.tstep = (1/data.cnt.srate) * 1000;
  
  data.time = [stime:data.tstep:etime]';
  
  set(VIEWCNT.axes,'XLim',[data.time(1) data.time(end)]);
  
  % Baseline data on first 10 msec
  data.base = data.dispSPoint:data.dispSPoint + (data.cnt.srate / 1000) * 10;
  data.basemean = mean(data.cnt.volt(data.base,:));
  basemean = repmat(data.basemean,data.dispNPoints,1);
  volt = data.cnt.volt(data.dispSPoint:data.dispEPoint,:);
  volt = volt - basemean;
  
  plot(data.time,volt);
  
  eeg_plot_metric;
  
  %dataSD   = std(data);
  %dataVar  = var(data);
  %dataMean = repmat(dataMean,points,1);
  %dataVar  = repmat(dataVar ,points,1);
  %data = (data - dataMean) ./ dataVar;
  %increment = [1:nChan];
  %data = data + repmat(increment,points,1);
  if isfield(data.cnt,'labels'),
    data.labels = data.cnt.labels;
  else
    data.labels = [];
  end
  if (length(data.labels) ~= data.nChan),
    data.labels = [];
  end
  %set(VIEWCNT.axes,'YTick',increment,'YTickLabel',labels);
  
  axpos = get(VIEWCNT.axes,'position');
  
  % Baseline data interface
  G.baseline = uicontrol(GUI,'style','edit','units','normalized',Font,...
    'Position',[.92 .8 .075 .05],...
    'TooltipString','Baseline (msec)',...
    'String',10,...
    'min',1,'max',(data.dispNPoints/data.cnt.srate)*1000,...
    'Callback','eeg_view_cnt('''',''plot'');');
  
  G.seconds = uicontrol(GUI,'style','edit','units','normalized',Font,...
    'Position',[.92 .7 .075 .05],...
    'TooltipString','Seconds / Page',...
    'String',sprintf('%8.4f',data.dispSeconds), ...
    'min',1,'max',data.points/data.cnt.srate,...
    'Callback','eeg_view_cnt('''',''plot'');');
  
  tpos = [ axpos(1,1)-.12 .01 .12        .05];
  spos = [ axpos(1,1)     .01 axpos(1,3) .05];
  
  hbuttonstep = .01; % move 1% of datafile
  hsliderstep = .05; % move 5% of datafile
  horizSliderStep = [ hbuttonstep hsliderstep ];
  
  G.Tslider = uicontrol(GUI,'style','text','units','normalized',Font,...
    'Position',tpos,...
    'TooltipString','% of CNT file',...
    'String',sprintf('%6.2f sec',1/data.cnt.srate),...
    'HorizontalAlignment', 'center');
  G.Slider = uicontrol(GUI,'style','slider','units','normalized',Font,...
    'Position',spos,...
    'sliderstep',horizSliderStep, ...
    'min',1,'max',data.cnt.numSamples,'value',1,...
    'Callback','eeg_view_cnt('''',''plot'');');
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Font.FontWeight = 'bold';
  
  G.exit = uicontrol(GUI,'style','pushbutton','units','normalized',Font,...
    'Position',[.925 .01 .07 .05],...
    'TooltipString','Close',...
    'String','EXIT',...
    'BackgroundColor',[0.75 0.0 0.0],...
    'ForegroundColor', [1 1 1], 'HorizontalAlignment', 'center',...
    'Callback','close gcbf;');
  
  VIEWCNT.data = data;
  VIEWCNT.handles = G;
  set(VIEWCNT.gui,'userdata',VIEWCNT);
  
  
  return
