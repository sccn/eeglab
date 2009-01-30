function viewedf(Action, ActOption)
% viewedf
% Display EDF (European Data Format) files
% This program requires Matlab Version >= 5.2 and should work with Matlab 6.x
%
% Command line functionality (note: all commands are case sensitive):
% viewedf : start program
% viewedf('OpenEDFFile', FILENAME) : open FILENAME
% viewedf('Goto', RECORD) : move to RECORD
% viewedf('GotoSecond', SECOND) : move to a certain SECOND
% viewedf('Prev') : move back one screen
% viewedf('Next') : move forward one screen
% viewedf('PrevFast') : move back 5 screens
% viewedf('NextFast') : move forward 5 screens
% viewedf('Close') : close viewer
%
% Keyboard shortcuts:
% '+' : move to next screen
% '-' : move to previous screen
% 'U' : scale all channels up
% 'D' : scale all channels down
%
% Version 3.04Alpha, 10/22/01
% (c) Herbert Ramoser (herbert.ramoser@arcs.ac.at)
%
% This Software is subject to the GNU public license.

% Changes:
% 01/27/98 Ramoser: Plugins may replace EDF data, optional parameter string
%   is passed to the plugin
% 01/30/98 Ramoser: plotting of EDF.Valid information
% 02/02/98 Ramoser: print option added
% 02/18/98 Ramoser: Stairplot-feature added, bug for hypnograms removed
% 06/16/98 Ramoser: bugs of Matlab 5.2 fixed (return, inputdlg)
%          display multiple EDF files
%          fixed order of channels
%          export of matrices to workspace
%          increase in plotting speed
% 09/18/98 Ramoser: new EDF open and read functions included
%          upside down display of plots
%          use of cell arrays to store EDF data
%          data passed to plugins has changed
% 11/01/98 Woertz: reset display after calls of plugin menus
% 11/23/98 Ramoser: give only one file to plugin
%          handling of plugins when EDF files are closed
% 11/30/98 Ramoser: GDF functionality added
% 12/02/98 Ramoser: optionally display plot Y tick-labels
% 12/03/98 Ramoser: fixed dialog boxes under X (drawnow added)
%          upside-down plotting bug removed
% 12/07/98 Ramoser: scaling buttons added
% 12/09/98 Ramoser: physical dimension added to Y-label
% 02/18/99 Ramoser: problems with multiple files, multiple record lengths
%          and scrollbar removed
% 03/04/99 Ramoser: Assign matrix - plugin bug removed
% 03/26/99 Ramoser: GDF 0.12 added, set properties of all channels
%          simultaneusly ('Apply to all' button)
% 03/31/99 Ramoser: time tracking line added
% 04/01/99 Ramoser: print problems resolved
% 04/29/99 Ramoser: do not reset display properties after a call to
%          plugin-menu
% 05/04/99 Ramoser: time tracking line bugs removed
% 05/12/99 Ramoser: file positioning information added (EDFHead.AS.startrec
%          & numrec)
% 03/17/00 Ramoser: add command line functionality, add 'Goto second'
% 10/22/01 Ramoser: some updates for Matlab 6
% 01 Sep 2008 Schloegl: use sopen, sclose 

FastPageIncrement = 5;
PageIncrement = 1;
warning('off');

switch nargin,
 case 0,
  if ~isempty(findobj('Tag', 'ViewEDFFigure'))
    error('Only one copy of VIEWEDF may be run');
  end
  LocalInitWindow;
 case 1,
  switch Action
   case 'OpenEDFFile',
    LocalEDFOpen;
   case 'CloseEDFFile',
    LocalEDFClose;
   case 'Repaint'
    LocalRepaint;
   case 'AddPlugin'
    LocalAddPlugin;
   case 'RemovePlugin'
    LocalRemovePlugin;
   case 'PluginMenu'
    LocalPluginMenu;
   case 'PrevFast'
    LocalChangeRecord(-FastPageIncrement,0);
   case 'Prev'
    LocalChangeRecord(-PageIncrement,0);
   case 'NextFast'
    LocalChangeRecord(FastPageIncrement,0);
   case 'Next'
    LocalChangeRecord(PageIncrement,0);
   case 'HScroll'
    LocalHScroll;
   case 'ToggleUpdatePlugin'
    LocalToggleUpdatePlugin;
   case 'ToggleShowRange'
    LocalToggleShowRange;
   case 'ToggleShowCursor'
    LocalToggleShowCursor;
   case 'RecordProp'
    LocalRecordProp;
   case 'Goto'
    LocalGotoRecord;
   case 'GotoSecond'
    LocalGotoSecond;
   case 'Records'
    LocalNumRecords;
   case 'FileInfo'
    LocalFileInfo;
   case 'Channels'
    LocalSelectChannels;
   case 'About'
    LocalAbout;
   case 'Help'
    LocalHelp;
   case 'KeyPress'
    LocalKeyPress;
   case 'Print'
    LocalPrint;
   case 'AssignMat'
    LocalAssignMat;
   case 'ScaleUp'
    LocalRescale('up');
   case 'ScaleDown'
    LocalRescale('down');
   case 'Close'
    LocalCloseViewEDF;
   otherwise
    error('VIEWEDF must be called without parameters');
  end
 case 2,
  switch Action
   case 'OpenEDFFile',
    LocalEDFOpen(ActOption);
   case 'FileInfo'
    LocalFileInfo(ActOption);
   case 'Channels'
    LocalSelectChannels(ActOption);
   case 'MoveCursor'
    LocalMoveCursor(ActOption);
   case 'Goto'
    LocalGotoRecord(ActOption);
   case 'GotoSecond'
    LocalGotoSecond(ActOption);
   otherwise
    error('Unknown argument!');
  end
 otherwise
  error('Unknown argument!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% LocalInitWindow
% initializes the window (display buttons, ...)
function LocalInitWindow()
fig=figure( ...
    'NumberTitle', 'off', ...
    'CloseRequestFcn', 'viewedf Close', ...
    'IntegerHandle', 'off', ...
    'KeyPressFcn', 'viewedf KeyPress', ...
    'MenuBar', 'none', ...
    'Units', 'Normalized', ...
    'Color', [0.9, 0.9, 0.9], ...
    'ResizeFcn', 'viewedf Repaint', ...
    'PaperPositionMode', 'Auto', ...
    'PaperType', 'A4', ...
    'PaperUnits', 'Normalized', ...
    'Tag', 'ViewEDFFigure', ...
    'Name', 'EDF File Viewer - (C) 1998-2001 DPMI, 2008 HCI, Graz University of Technology');

Data=get(fig,'UserData');

% Menus
mh = uimenu('Label', '&File');
uimenu(mh, ...
    'Label', '&Open EDF', ...
    'Callback', 'viewedf OpenEDFFile');
uimenu(mh, ...
    'Label', '&Close EDF', ...
    'Callback', 'viewedf CloseEDFFile');
uimenu(mh, ...
    'Label', '&Add Plugin', ...
    'Callback', 'viewedf AddPlugin');
uimenu(mh, ...
    'Label', '&Remove Plugin', ...
    'Callback', 'viewedf RemovePlugin');
uimenu(mh, ...
    'Label', '&Print', ...
    'Separator', 'on', ...
    'Callback', 'viewedf Print');
uimenu(mh, ...
    'Label', 'Assign &Matrix', ...
    'Callback', 'viewedf AssignMat');
uimenu(mh, ...
    'Label', 'E&xit', ...
    'Separator', 'on', ...
    'Callback', 'viewedf Close');

mh = uimenu('Label', '&Display');
uimenu(mh, ...
    'Label', 'File &info', ...
    'Callback', 'viewedf FileInfo');
uimenu(mh, ...
    'Label', '&Goto Record', ...
    'Callback', 'viewedf Goto');
uimenu(mh, ...
    'Label', '&Goto Second', ...
    'Callback', 'viewedf GotoSecond');
uimenu(mh, ...
    'Label', '&Records on Screen', ...
    'Callback', 'viewedf Records');
uimenu(mh, ...
    'Label', '&Channels', ...
    'Callback', 'viewedf Channels');

mh = uimenu('Label', '&Options');
uimenu(mh, ...
    'Label', '&Update Plugin', ...
    'Checked', 'on', ...
    'Callback', 'viewedf ToggleUpdatePlugin');
uimenu(mh, ...
    'Label', 'Show &Range', ...
    'Checked', 'off', ...
    'Callback', 'viewedf ToggleShowRange');
uimenu(mh, ...
    'Label', 'Show &Cursor', ...
    'Checked', 'off', ...
    'Callback', 'viewedf ToggleShowCursor');
Data.Display.PluginMenu.Main = uimenu(mh, ...
    'Label', '&Plugin');

mh = uimenu('Label', '&?');
uimenu(mh, ...
    'Label', '&Help', ...
    'CallBack', 'viewedf Help');
uimenu(mh, ...
    'Label', '&About', ...
    'CallBack', 'viewedf About');

% Strings
textx = 0.02;
textheight = LocalGetFontHeight;
texty = 1 - textheight / 4;
Data.Display.Strings.CurrTime = uicontrol(fig, ... % Time string
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'String', '', ...
    'HorizontalAlignment', 'left', ...
    'FontUnits', 'normalized', ...
    'Position', [textx, texty - textheight, 0.2, ...
      textheight]);
Data.Display.Strings.TotTime = uicontrol(fig, ... % Total time string
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'String', 'Total length : ', ...
    'HorizontalAlignment', 'left', ...
    'Position', [textx, texty - textheight, 0.2, ...
      textheight]);
Data.Display.Strings.DispTime = uicontrol(fig, ... % Displayed time stretch string
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'String', 'Displayed : ', ...
    'HorizontalAlignment', 'left', ...
    'Position', [textx + 0.22, texty - textheight, 0.2, ...
      textheight]);
Data.Display.Strings.FileName = uicontrol(fig, ... % File string
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'String', 'File : ', ...
    'HorizontalAlignment', 'left', ...
    'Position', [textx + 0.44, texty - textheight, 0.2, ...
      textheight]);
Data.Display.Strings.NumString = 4;

% Build data structure
% settings for the display (where to plot things, ...)
Data.Display.Axes = [];       % nothing to display
Data.Display.YSpacing = 0.0025;   % spacing between plots
Data.EDF = [];
Data.Plugin = [];
Data.UpdatePlugin = 1;
Data.ShowRange = 0;
Data.ShowCursor = 0;
Data.Display.Cursor = [];
set(fig,'UserData',Data);
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% LocalResetDisplay
% reset display variables
function LocalResetDisplay()
Data = get(findobj('Tag', 'ViewEDFFigure'), 'UserData');
if length(Data.EDF)==0
  return;
end
%EDF data
MaxDur = max(LocalGetEDFInfo('Dur', Data.EDF));
for i=1:length(Data.EDF)
  Data.Display.EDF(i).ShowRecords = [0, round(MaxDur / Data.EDF(i).Head.Dur)];
  [Data.EDF(i).Record, Data.EDF(i).Head] = ...
      LocalEDFRead(Data.EDF(i).Head, Data.Display.EDF(i).ShowRecords);
  Data.Display.EDF(i).ShowSignals = 1:Data.EDF(i).Head.NS;
  Data.Display.EDF(i).DisplayMin = Data.EDF(i).Head.PhysMin;
  Data.Display.EDF(i).DisplayMax = Data.EDF(i).Head.PhysMax;
  Data.Display.EDF(i).StairPlot = zeros(1, Data.EDF(i).Head.NS);
  Data.Display.EDF(i).UpDownPlot = Data.EDF(i).Head.PhysMin>Data.EDF(i).Head.PhysMax; %zeros(1, Data.EDF(i).Head.NS);
end

% Plugin data
for i = 1:length(Data.Plugin)
  [Data.Plugin(i).EDF, Data.Plugin(i).UserData] = ...
      feval(Data.Plugin(i).Name, Data.EDF(1), Data.Plugin(i).UserData, ...
      'Reset');
  Data.Display.Plugin(i).ShowRecords = [0, round(MaxDur / Data.Plugin(i).Head.Dur)];
  Data.Display.Plugin(i).ShowSignals = 1:Data.Plugin(i).Head.NS;
  Data.Display.Plugin(i).DisplayMin = Data.Plugin(i).EDF.Head.PhysMin;
  Data.Display.Plugin(i).DisplayMax  = ...
      Data.Plugin(i).EDF.Head.PhysMax;
  Data.Display.Plugin(i).StairPlot = zeros(1, Data.Plugin(i).Head.NS);
  Data.Display.Plugin(i).UpDownPlot = Data.EDF(i).Head.PhysMin>Data.EDF(i).Head.PhysMax; %zeros(1, Data.Plugin(i).Head.NS);
end

% update all strings
tsec = max(LocalGetEDFInfo('FileDur', Data.EDF));
tmin = floor(tsec/60);
tsec = rem(tsec,60);
th = floor(tmin/60);
tmin = rem(tmin,60);
set(Data.Display.Strings.TotTime, 'String', ...
  sprintf('Total length : %02d:%02d:%02d', th,tmin,tsec));
fname = 'File : ';
% add '+' for additional files
for i = 1:(length(Data.EDF) - 1)
  fname = [fname Data.EDF(i).Head.FILE.Name '.' Data.EDF(i).Head.FILE.Ext ...
        ' + ']; 
end
fname = [fname Data.EDF(end).Head.FILE.Name '.' Data.EDF(end).Head.FILE.Ext];
set(Data.Display.Strings.FileName, 'String', fname);
set(findobj('Tag', 'ViewEDFFigure'), 'UserData', Data);
LocalRepaint(0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% LocalRepaint
% initializes the window (display buttons, ...)
function LocalRepaint(KeepOldPlots)
figure(findobj('Tag', 'ViewEDFFigure'));
Data = get(gcf, 'UserData');
fhght = LocalGetFontHeight;
if length(Data.EDF)==0
  return;
end
if nargin < 1
  KeepOldPlots = 0;
end
KeepOldPlots = (KeepOldPlots ~= 0);

% set display size (area left for plots)
Data.Display.DrawRect = [0.005 1.5*fhght 0.85 1-1.5*fhght]; 

% get width of figure
unit = get(gcf, 'Units');
set(gcf, 'Units', 'Pixels');
xpix = get(gcf, 'Position');
set(gcf, 'Units', unit);
xpix = xpix(3);

% calculate total number of plots 
TotPlots = 0;
for i=1:length(Data.EDF)
  TotPlots = TotPlots + length(Data.Display.EDF(i).ShowSignals);
end
for i=1:length(Data.Plugin)
  TotPlots = TotPlots + length(Data.Display.Plugin(i).ShowSignals);
end

% update plugin-data
LocalWatchOn;
if Data.UpdatePlugin | ~KeepOldPlots
  for i = 1:length(Data.Plugin)
    [Data.Plugin(i).EDF, Data.Plugin(i).UserData] = ...
        feval(Data.Plugin(i).Name, Data.EDF(Data.Plugin(i).EDFFile), ...
        Data.Plugin(i).UserData);
  end
end
LocalWatchOff;

Opt.YSpace = Data.Display.YSpacing;
Opt.TxtHeight = fhght;
Opt.TextWidth = 1 - Data.Display.DrawRect(3) - Data.Display.DrawRect(1) - ...
    0.01; 
Opt.TextX = Data.Display.DrawRect(1) + Data.Display.DrawRect(3) + 0.005;
Opt.ButHeight = Opt.TxtHeight * 1.2;
Opt.ValidHeight = 0.005;
Opt.FigX = Data.Display.DrawRect(1);
Opt.FigWidth = Data.Display.DrawRect(3);
if Data.ShowRange
  Opt.LabelWidth = min([0.1, 4*LocalGetFontWidth]);
else
  Opt.LabelWidth = 0;
end
Opt.FigY = Data.Display.DrawRect(2) + Data.Display.DrawRect(4) - ...
    1.5*Opt.TxtHeight - Opt.ValidHeight;
Opt.FigHeight = (Data.Display.DrawRect(4) - 1.5*Opt.TxtHeight - ...
    Opt.ValidHeight)  / TotPlots - Opt.YSpace;
Opt.TxtHeight = min([(Opt.FigHeight-Opt.YSpace)/2, Opt.TxtHeight]);
Opt.ButHeight = min([(Opt.FigHeight-Opt.YSpace)/2, Opt.ButHeight]);
Opt.XPixels = 4*xpix;
Opt.ShowRange = Data.ShowRange;

% clear old buttons and labels
if ~KeepOldPlots & ~isempty(Data.Display.Axes) 
  delete(Data.Display.HScrollBar);
  if ~isempty(Data.Display.Cursor); 
    delete(Data.Display.Cursor.Line(:));
    delete(Data.Display.Cursor.Menu.Text);
    delete(Data.Display.Cursor.Menu.Menu);
    Data.Display.Cursor = [];
  end
  delete(Data.Display.Axes(:).RecButton);
  delete(Data.Display.Axes(:).ScaleUpButton);
  delete(Data.Display.Axes(:).ScaleDownButton);
  delete(Data.Display.Axes(:).RecLabel);
  delete(Data.Display.Axes(:).MaxLabel);
  delete(Data.Display.Axes(:).MinLabel);
  delete(Data.Display.Axes(:).DimLabel);
  delete(Data.Display.Axes(:).PlotLine);
  delete(Data.Display.Axes(:).Plot);
  Data.Display.Axes = [];
end

% draw HScrollBar
FileDur = LocalGetEDFInfo('FileDur',Data.EDF);
[MaxVal, MFInd] = max(FileDur);
Dur = LocalGetEDFInfo('Dur', Data.EDF);
[MaxDur, MDInd] = max(Dur);
Len = Data.Display.EDF(MFInd).ShowRecords(1);
if ~KeepOldPlots
  Data.Display.HScrollBar = uicontrol(gcf, ...
      'Style', 'Slider', ...
      'Units', 'Normalized', ...
      'Position', [Data.Display.DrawRect(1), 0, ...
        Data.Display.DrawRect(3), LocalGetFontHeight], ...
      'SliderStep', [MaxDur/MaxVal 5*MaxDur/MaxVal] * ...
      Data.Display.EDF(MDInd).ShowRecords(2), ...
      'Min', 0, ...
      'Max', 1, ...
      'Value', Data.Display.EDF(MFInd).ShowRecords(1)/MaxVal*Dur(MFInd), ...
      'CallBack', 'viewedf HScroll');
else
  set(Data.Display.HScrollBar, ...
      'Value', Data.Display.EDF(MFInd).ShowRecords(1)/MaxVal*Dur(MFInd)); 
end

%create new buttons, plots and labels
cnt = 0;
for j=1:length(Data.EDF)
  for i=1:length(Data.Display.EDF(j).ShowSignals)
    cnt = cnt+1;
    showsig = Data.Display.EDF(j).ShowSignals(i);
    if ~KeepOldPlots  
      % make new plots, buttons, labels
      Opt.YPos = Opt.FigY - Opt.YSpace*(cnt-1) - Opt.FigHeight * cnt;
      Opt.DisplayMin = Data.Display.EDF(j).DisplayMin(showsig);
      Opt.DisplayMax = Data.Display.EDF(j).DisplayMax(showsig);
      Opt.UserData = { 'EDF', showsig, j};
      Opt.Label = Data.EDF(j).Head.Label(showsig,:);
      %Opt.YLabel = deblank(Data.EDF(j).Head.PhysDim(showsig,:));
      Opt.YLabel = deblank(Data.EDF(j).Head.PhysDim{showsig});
      yd = Data.EDF(j).Record{showsig};
      if (length(yd) ~= 1)
        yd = [yd(:); NaN];
      end
      
      Temp = ...
          LocalPlotNewData(yd, ...
          Data.Display.EDF(j).UpDownPlot(showsig), ...
          Data.Display.EDF(j).StairPlot(showsig), Opt);
      if cnt == 1
          Data.Display.Axes = Temp;
      else
          Data.Display.Axes(cnt) = Temp;
      end
    else
      % set new ydata
      yd = Data.EDF(j).Record{showsig} * ...
          (Data.Display.EDF(j).UpDownPlot(showsig)-0.5) * -2;
      if (length(yd) == 1)
        yd = yd([1 1]);
      else
        yd = yd(1:ceil(length(yd)/Opt.XPixels):length(yd));
        yd = [yd(:); NaN];
        if Data.Display.EDF(j).StairPlot(showsig)
          yd = yd(floor(1:0.5:length(yd)));
        end
      end
      set(Data.Display.Axes(cnt).PlotLine, 'YData', yd);
    end
  end
end

% plot plugin-data
for j = 1:length(Data.Plugin)
  for i=1:length(Data.Display.Plugin(j).ShowSignals)
    cnt = cnt+1;
    showsig = Data.Display.Plugin(j).ShowSignals(i);
    if ~KeepOldPlots
      % make new plots, buttons, labels
      Opt.YPos = Opt.FigY - Opt.YSpace*(cnt-1) - Opt.FigHeight * cnt;
      Opt.DisplayMin = Data.Display.Plugin(j).DisplayMin(showsig);
      Opt.DisplayMax = Data.Display.Plugin(j).DisplayMax(showsig);
      Opt.UserData = { 'PLUGIN', showsig, j};
      Opt.Label = Data.Plugin(j).EDF.Head.Label(showsig,:);
      %Opt.YLabel = deblank(Data.Plugin(j).EDF.Head.PhysDim(showsig,:));
      Opt.YLabel = deblank(Data.Plugin(j).EDF.Head.PhysDim{showsig});
      yd = Data.Plugin(j).EDF.Record{showsig};
      if (length(yd) ~= 1)
        yd = [yd(:); NaN];
      end
      Data.Display.Axes(cnt) = ...
          LocalPlotNewData(yd, ...
          Data.Display.Plugin(j).UpDownPlot(showsig), ...
          Data.Display.Plugin(j).StairPlot(showsig), Opt);
    else
      yd = Data.Plugin(j).EDF.Record{showsig} * ...
          (Data.Display.Plugin(j).UpDownPlot(showsig)-0.5) * -2;
      if (length(yd) == 1)
        yd = yd([1 1]);
      else
        yd = yd(1:ceil(length(yd)/Opt.XPixels):length(yd));
        yd = [yd(:); NaN];        
        if Data.Display.Plugin(j).StairPlot(showsig)
          yd = yd(floor(1:0.5:length(yd)));
        end
      end
      set(Data.Display.Axes(cnt).PlotLine, 'YData', yd);
    end
  end
end

% change size of all strings
set(Data.Display.Strings.CurrTime, ...
    'Position', [Opt.TextX, 0, Opt.TextWidth, fhght]);
pos = get(Data.Display.Strings.TotTime, 'Position');
set(Data.Display.Strings.TotTime, ...
    'Position', [pos(1), 1 - 1.25 * fhght, pos(3), fhght]);
pos = get(Data.Display.Strings.DispTime, 'Position');
set(Data.Display.Strings.DispTime, ...
    'Position', [pos(1), 1 - 1.25 * fhght, pos(3), fhght]);
pos = get(Data.Display.Strings.FileName, 'Position');
ext = get(Data.Display.Strings.FileName, 'Extent'); 
set(Data.Display.Strings.FileName, ...
    'Position', [pos(1), 1 - 1.25 * fhght, ext(3) + 0.01, fhght]);

% change size of HScrollBar
pos = get(Data.Display.HScrollBar, 'Position');
set(Data.Display.HScrollBar, ...
    'Position', [pos(1) pos(2) pos(3) fhght]); 

% draw cursor
if Data.ShowCursor & ~KeepOldPlots
  Data.Display.Cursor = LocalDrawCursor(Data.Display);
end

% display time
tsec = Data.EDF(1).Head.Dur * Data.Display.EDF(1).ShowRecords(1);
tmin = floor(tsec/60);
tsec = rem(tsec,60);
th = floor(tmin/60);
tmin = rem(tmin,60);
set(Data.Display.Strings.CurrTime, 'String', ...
    sprintf('Time %02d:%02d:%02d', th,tmin,tsec));
tsec = Data.EDF(1).Head.Dur * Data.Display.EDF(1).ShowRecords(2);
tmin = floor(tsec/60);
tsec = rem(tsec,60);
th = floor(tmin/60);
tmin = rem(tmin,60);
set(Data.Display.Strings.DispTime, 'String', ...
    sprintf('Displayed : %02d:%02d:%02d', th,tmin,tsec));
set(gcf, 'UserData', Data);
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LocalPlotNewData
% plots a single data vector
function [ploth] =  LocalPlotNewData(Data, UpDownPlot, StairPlot, Opt)
len = length(Data);
ploth.Plot = subplot('Position', ...
    [Opt.FigX + Opt.LabelWidth, Opt.YPos, Opt.FigWidth - Opt.LabelWidth, ...
      Opt.FigHeight - Opt.YSpace]);
if (len == 1)
  xd = [1 2];
  yd = Data([1 1]);
else
  yd = Data(1:ceil(len/Opt.XPixels):len);
  len = length(yd);
  if StairPlot 
    xd = ceil(1:0.5:len);
    yd = yd(floor(1:0.5:len));
  else
    xd = 1:len;
  end
end

if UpDownPlot
  yd = yd * -1;
  tmp = Opt.DisplayMin;
  Opt.DisplayMin = -Opt.DisplayMax;
  Opt.DisplayMax = -tmp;
end

ploth.PlotLine = plot(xd, yd);
set(gca, ...
    'XLim', [1 max([2, len])], ...
    'YLim', sort([Opt.DisplayMin Opt.DisplayMax]));
if ~Opt.ShowRange 
  set(gca, ...
      'XTickLabel', '', ...
      'YTickLabel', '', ...
      'XTick', [], ...
      'YTick', []);
  ploth.MaxLabel = [];
  ploth.MinLabel = [];
  ploth.DimLabel = [];
else
  % get TickLabels used by MatLab 
  set(gca, ...
      'YTick', [Opt.DisplayMin Opt.DisplayMax]);
  labelstr = get(gca, 'YTickLabel');
  labelstr = {deblank(labelstr(1,:)) deblank(labelstr(2,:))};
  
  if UpDownPlot
    for i=1:2
      if labelstr{i}(1) == '-'
        labelstr{i} = labelstr{i}(2:length(labelstr{i}));
      else
        labelstr{i} = ['-' labelstr{i}];
      end
    end
  end
    
  set(gca, 'YTickLabel', '');
  % place TickLabels 
  ploth.MinLabel = text(0, 0, ...
      labelstr{1}, ...
      'Parent', gca, ...
      'HorizontalAlignment', 'right', ...
      'VerticalAlignment', 'baseline', ...
      'Units', 'normalized');
  ploth.MaxLabel = text(0, 1, ...
      labelstr{2}, ...
      'Parent', gca, ...
      'HorizontalAlignment', 'right', ...
      'VerticalAlignment', 'cap', ...
      'Units', 'normalized');
  % place physical dimension
  ploth.DimLabel = text(0, 0.5, ...
      ['[' Opt.YLabel ']'], ...
      'Parent', gca, ...
      'HorizontalAlignment', 'right', ...
      'VerticalAlignment', 'middle', ...
      'Units', 'normalized');
end

ButWidth = Opt.TextWidth / 3;
ploth.RecButton = uicontrol(gcf, ...
    'Style', 'PushButton', ...
    'Units', 'Normalized', ...
    'String', '?', ...
    'Position', [Opt.TextX, Opt.YPos, ButWidth, Opt.ButHeight], ...
    'Callback', 'viewedf RecordProp', ...
    'UserData', Opt.UserData);
ploth.ScaleUpButton = uicontrol(gcf, ...
    'Style', 'PushButton', ...
    'Units', 'Normalized', ...
    'String', '+', ...
    'Position', [Opt.TextX+ButWidth, Opt.YPos, ButWidth, Opt.ButHeight], ...
    'Callback', 'viewedf ScaleUp', ...
    'UserData', Opt.UserData);
ploth.ScaleDownButton = uicontrol(gcf, ...
    'Style', 'PushButton', ...
    'Units', 'Normalized', ...
    'String', '-', ...
    'Position', [Opt.TextX+2*ButWidth, Opt.YPos, ButWidth, Opt.ButHeight], ...
    'Callback', 'viewedf ScaleDown', ...
    'UserData', Opt.UserData);
ploth.RecLabel = uicontrol(gcf, ...
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'String', Opt.Label, ...
    'Position', [Opt.TextX, Opt.YPos + Opt.FigHeight - Opt.YSpace - ...
      Opt.TxtHeight, Opt.TextWidth, Opt.TxtHeight]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LocalAssignMat
% Save display to Matrix
function LocalAssignMat()
Data = get(findobj('Tag', 'ViewEDFFigure'), 'UserData');
if length(Data.EDF)==0
  errordlg('No EDF-File is open!', 'Error');
  return;
end

answer = inputdlg({'Target name'}, 'Save to cell-array', 1, {''});
if (length(answer) == 0)
  return;
end
cnt = 0;
for j=1:length(Data.EDF)
  for i=1:length(Data.Display.EDF(j).ShowSignals)
    cnt = cnt+1;
    showsig = Data.Display.EDF(j).ShowSignals(i);
    res{cnt} = Data.EDF(j).Record{showsig};
  end
end
for j=1:length(Data.Plugin)
  for i=1:length(Data.Display.Plugin(j).ShowSignals)
    cnt = cnt+1;
    Data.Plugin(j)
    showsig = Data.Display.Plugin(j).ShowSignals(i);
    res{cnt} = Data.Plugin(j).EDF.Record{showsig};
  end
end
assignin('base', answer{1}, res);
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LocalGetEDFInfo
% get information of all open EDF files
function [res] = LocalGetEDFInfo(Which, EDF)

res = [];
for k=1:length(EDF),
	switch Which
	case 'NRec'
		res = [res,EDF(k).Head.NRec];
	case 'Dur'
		res = [res,EDF(k).Head.SPR / EDF(k).Head.SampleRate]; 
	case 'SPR'
		res = [res,EDF(k).Head.SPR];
	case 'FileDur'
		res = [res,EDF(k).Head.NRec * EDF(k).Head.SPR / EDF(k).Head.SampleRate]; 
	end
end; 	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LocalAddPlugin
% Add plugin function
function LocalAddPlugin()
Data = get(findobj('Tag', 'ViewEDFFigure'), 'UserData');
if length(Data.EDF)==0
  errordlg('No EDF-File is open!', 'Error');
  return;
end

res=inputdlg({ 'Enter Plugin name','Label (Optional)', ...
    'Parameters (Optional)'}, 'Select Plugin', 1, { '', '', ''});
if isempty(res)
  return;
end
plugname = res{1};
pluglabel = res{2};
plugopt = res{3};
if isempty(pluglabel)
  pluglabel = plugname;
end
if ~exist(plugname, 'file');
  errordlg(sprintf('Plugin %s can not be found in standard search path', ...
      upper(plugname)), 'Plugin Error');
  return
end;

% select EDF file if several files are open
infile = 1;
if (length(Data.EDF) > 1)
  [infile, cancelled] = LocalSelectEDFFile({'EDF-data passed to plugin', ...
        'Select'}, Data.EDF);
  if cancelled
    return
  end
end

% set all variables
ind = length(Data.Plugin) + 1;
Data.Plugin(ind).Name = plugname;
Data.Plugin(ind).Label = pluglabel;
Data.Plugin(ind).UserData = plugopt;
Data.Plugin(ind).EDFFile = infile;
LocalWatchOn;
[Data.Plugin(ind).EDF, Data.Plugin(ind).UserData] = ...
    feval(Data.Plugin(ind).Name, Data.EDF(infile), ...
    Data.Plugin(ind).UserData, 'Reset');
LocalWatchOff;
Data.Display.Plugin(ind).ShowSignals = 1:Data.Plugin(ind).EDF.Head.NS;
Data.Display.Plugin(ind).DisplayMin = Data.Plugin(ind).EDF.Head.PhysMin;
Data.Display.Plugin(ind).DisplayMax = Data.Plugin(ind).EDF.Head.PhysMax;
Data.Display.Plugin(ind).StairPlot = ...
    zeros(1, Data.Plugin(ind).EDF.Head.NS);
Data.Display.Plugin(ind).UpDownPlot = ...
    zeros(1, Data.Plugin(ind).EDF.Head.NS);
MaxDur = max(LocalGetEDFInfo('Dur', Data.EDF));
Data.Display.Plugin(ind).ShowRecords = [0 ...
      round(Data.Plugin(ind).EDF.Head.Dur / MaxDur)];

% Add menu entry to options menu
Data.Display.PluginMenu.Sub(ind) = uimenu(Data.Display.PluginMenu.Main, ...
    'Label', [pluglabel, ' (File ', upper(plugname), ')'], ...
    'UserData', ind, ...
    'Callback', 'viewedf PluginMenu');
set(findobj('Tag', 'ViewEDFFigure'), 'UserData', Data);
LocalRepaint(0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% LocalRemovePlugin
% Remove plugin function
function LocalRemovePlugin(whichplugin)
Data = get(findobj('Tag', 'ViewEDFFigure'), 'UserData');

if nargin == 0
  % determine which plugin to remove
  if length(Data.EDF)==0
    errordlg('No EDF-File is open!', 'Error');
    return;
  end
  if length(Data.Plugin) == 0
    errordlg('No Plugins are loaded!', 'Error');
    return;
  end

  dlgh = dialog(...
      'Name', 'Remove Plugin', ...
      'CloseRequestFcn', 'set(gcf,''UserData'',''Cancel'');uiresume;');
  dlgpos = get(dlgh, 'Position');
  set(dlgh, 'Position', [dlgpos(1),dlgpos(2),200,100]);

  PlugNames = {};
  % get Plugin names
  for i = 1:length(Data.Plugin)
    PlugNames = {PlugNames{:}, [Data.Plugin(i).Label, ' (File ', ...
            upper(Data.Plugin(i).Name), ')']};
  end
  
  poph = uicontrol(dlgh, ...
      'Style', 'Popup', ...
      'Units', 'Normalized', ...
      'String', PlugNames, ...
      'Value', 1, ...
      'Position', [0.05, 0.5, 0.9, 0.35]);
  % buttons
  uicontrol(dlgh, ...
      'Style', 'PushButton', ...
      'Units', 'Normalized', ...
      'String', 'Remove', ...
      'Position', [0.1, 0.05, 0.3, 0.3], ...
      'Callback', 'set(gco,''UserData'',''OK'');uiresume;');
  uicontrol(dlgh, ...
      'Style', 'PushButton', ...
      'Units', 'Normalized', ...
      'String', 'Cancel', ...
      'Position', [0.6, 0.05, 0.3, 0.3], ...
      'Callback', 'set(gco,''UserData'',''Cancel'');uiresume;');
  
  drawnow;
  uiwait(dlgh);
  whichplugin = get(poph, 'Value');
  selbutton = get(gco,'UserData');
  delete(dlgh);
  if strcmp(selbutton,'Cancel')
    return
  end
end

delete(Data.Display.PluginMenu.Sub(whichplugin));
for i = whichplugin:length(Data.Plugin)-1
  Data.Display.PluginMenu.Sub(i) = Data.Display.PluginMenu.Sub(i+1);
  set(Data.Display.PluginMenu.Sub(i), 'UserData', i);
  Data.Display.Plugin(i) = Data.Display.Plugin(i+1);
  Data.Plugin(i) = Data.Plugin(i+1);
end
ind = length(Data.Plugin) - 1;
if ind == 0
  Data.Plugin = [];
  Data.Display.Plugin = [];
  Data.Display.PluginMenu.Sub = [];
else
  Data.Plugin = Data.Plugin(1:ind);
  Data.Display.Plugin = Data.Display.Plugin(1:ind);
  Data.Display.PluginMenu.Sub = Data.Display.PluginMenu.Sub(1:ind);
end
set(findobj('Tag', 'ViewEDFFigure'), 'UserData', Data);
LocalRepaint;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LocalPluginMenu
% Call plugin for setup
function LocalPluginMenu()
ind = get(gcbo, 'UserData');
Data = get(findobj('Tag', 'ViewEDFFigure'), 'UserData');
[Data.Plugin(ind).EDF, Data.Plugin(ind).UserData] = ...
    feval(Data.Plugin(ind).Name, Data.EDF(Data.Plugin(ind).EDFFile), ...
    Data.Plugin(ind).UserData, 'Menu');
% Uncomment to reset display properties after a call of the plugin-menu
%Data.Display.Plugin(ind).ShowSignals = 1:Data.Plugin(ind).EDF.Head.NS;
%Data.Display.Plugin(ind).DisplayMin = Data.Plugin(ind).EDF.Head.PhysMin;
%Data.Display.Plugin(ind).DisplayMax = Data.Plugin(ind).EDF.Head.PhysMax;
%Data.Display.Plugin(ind).StairPlot = ...
%    zeros(1, Data.Plugin(ind).EDF.Head.NS);
%Data.Display.Plugin(ind).UpDownPlot = ...
%    zeros(1, Data.Plugin(ind).EDF.Head.NS);
set(findobj('Tag', 'ViewEDFFigure'), 'UserData',Data);
LocalRepaint(0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% LocalChangeRecord
% Change displayed page
function LocalChangeRecord(increment, changetype)
Data = get(findobj('Tag', 'ViewEDFFigure'), 'UserData');
if length(Data.EDF) == 0
  errordlg('No EDF-File is open!', 'Error');
  return;
end
Dur = LocalGetEDFInfo('Dur',Data.EDF);
[MDur, MDInd] = max(Dur);
FLen = max(LocalGetEDFInfo('FileDur', Data.EDF));
[MFLen, MFInd] = max(FLen);

RelInc = round(ones(size(Dur)) * max(Dur) ./ Dur);
Inc = increment * RelInc;

for i = 1:length(Data.EDF)
  switch changetype
    case 0 % incremental
      Data.Display.EDF(i).ShowRecords(1) = Data.Display.EDF(i).ShowRecords(1) ...
          + Inc(i) * Data.Display.EDF(MDInd).ShowRecords(2); 
    case 1 % absolute
      Data.Display.EDF(i).ShowRecords(1) = Inc(i);
  end
end

% check whether we are moving to far
temp = Data.Display.EDF(MDInd).ShowRecords(1);
if temp < 0
  for i = 1:length(Data.EDF)
    Data.Display.EDF(i).ShowRecords(1) = Data.Display.EDF(i).ShowRecords(1) ...
        - temp * RelInc(i);
  end
end

temp = sum(Data.Display.EDF(MFInd).ShowRecords(1:2)) - ...
    Data.EDF(MFInd).Head.NRec; 
if temp > 0
  temp = ceil(temp / RelInc(MFInd));
  for i = 1:length(Data.EDF)
    Data.Display.EDF(i).ShowRecords(1) = Data.Display.EDF(i).ShowRecords(1) - ...
        temp * RelInc(i); 
  end
end
  
for i=1:length(Data.EDF)
  temp = sum(Data.Display.EDF(i).ShowRecords(1:2)) - ...
      Data.EDF(i).Head.NRec;
  whichrec(1) = Data.Display.EDF(i).ShowRecords(1);
  if temp > 0
    % do not read to much data
    whichrec(2) = Data.Display.EDF(i).ShowRecords(2) - temp;
    if (whichrec(2) < 0)
      whichrec(2) = 0;
    end
    [temp, Data.EDF(i).Head] = LocalEDFRead(Data.EDF(i).Head, whichrec);
    for j = 1:Data.EDF(i).Head.NS
      Data.EDF(i).Record{j} = zeros(Data.Display.EDF(i).ShowRecords(2) * ...
          Data.EDF(i).Head.SPR(j));
      if ~isempty(temp)
        Data.EDF(i).Record{j}(1:length(temp{j})) = temp{j};
      end
    end

  else
    [Data.EDF(i).Record, Data.EDF(i).Head] = ...
        LocalEDFRead(Data.EDF(i).Head, Data.Display.EDF(i).ShowRecords); 
  end
end
set(findobj('Tag', 'ViewEDFFigure'), 'UserData', Data);
LocalRepaint(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% LocalHSroll
% Goto record
function LocalHScroll
Data = get(findobj('Tag', 'ViewEDFFigure'), 'UserData');
MaxVal = max(LocalGetEDFInfo('FileDur',Data.EDF));
MaxDur = max(LocalGetEDFInfo('Dur', Data.EDF));
pos = get(gcbo, 'Value');
LocalChangeRecord(round(pos*MaxVal/MaxDur), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% LocalGotoRecord
% Goto record
function LocalGotoRecord(Pos)
Data = get(findobj('Tag', 'ViewEDFFigure'), 'UserData');
if length(Data.EDF)==0
  errordlg('No EDF-File is open!', 'Error');
  return;
end

[temp, ind] = max(LocalGetEDFInfo('Dur', Data.EDF));

if (nargin == 0)
  answer = inputdlg({'Select record'}, 'Change start record', 1, ...
                    {int2str(Data.Display.EDF(ind).ShowRecords(1))});
  if (length(answer) ~= 0)
    answer = str2num(answer{1});
    if ~isempty(answer)
      LocalChangeRecord(answer, 1);
    end
  end
else
  Pos = round(Pos);
  if Pos >= 0
    LocalChangeRecord(Pos, 1);
  end
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% LocalGotoSecond
% Goto second
function LocalGotoSecond(Pos)
Data = get(findobj('Tag', 'ViewEDFFigure'), 'UserData');
if length(Data.EDF)==0
  errordlg('No EDF-File is open!', 'Error');
  return;
end

[temp, ind] = max(LocalGetEDFInfo('Dur', Data.EDF));

if (nargin == 0)
  answer = inputdlg({'Select second'}, 'Change start time', 1, ...
    {int2str(round(Data.Display.EDF(ind).ShowRecords(1) / ...
      Data.EDF(ind).Head.Dur))});
  if (length(answer) ~= 0)
    answer = str2num(answer{1});
    if ~isempty(answer)
      LocalChangeRecord(answer / Data.EDF(ind).Head.Dur, 1);
    end
  end
else
  Pos = round(Pos / Data.EDF(ind).Head.Dur);
  if Pos >= 0
    LocalChangeRecord(Pos, 1);
  end
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% LocalRecordProp
% display channel properties
function LocalRecordProp()
Data = get(findobj('Tag', 'ViewEDFFigure'), 'UserData');
if length(Data.EDF)==0
  errordlg('No EDF-File is open!', 'Error');
  return;
end

whichbut = get(gcbo, 'UserData');
Rec = whichbut{2};
Num = whichbut{3};
switch whichbut{1}
  case 'EDF'
    RecHead = Data.EDF(Num).Head;
    RecDisp = Data.Display.EDF(Num);
    DataType = 'EDF File';
    TypeString = 'Filename :';
    OptString = [Data.EDF(Num).Head.FILE.Name '.' Data.EDF(Num).Head.FILE.Ext];
  case 'PLUGIN'
    RecHead = Data.Plugin(Num).EDF.Head;
    RecDisp = Data.Display.Plugin(Num);
    DataType = 'Plugin';
    TypeString = 'Name :';
    OptString = Data.Plugin(Num).Label;
end
dlgh = dialog(...
    'Name', 'Channel Information', ...
    'CloseRequestFcn', 'set(gcf,''UserData'',''Cancel'');uiresume;');
dlgpos = get(dlgh, 'Position');
set(dlgh, 'Position', [dlgpos(1),dlgpos(2),350,300]);
Fnthgt = LocalGetFontHeight;
yinc = 0.085;
% general information
linenum=1;
uicontrol(dlgh, ...
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold', ...
    'String', [ 'Record read from ' DataType], ...
    'Position', [0.02, 1-linenum*yinc, 0.96, Fnthgt ]);
% filename
linenum=linenum+1;
uicontrol(dlgh, ...
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'HorizontalAlignment', 'left', ...
    'String', TypeString, ...
    'Position', [0.02, 1-linenum*yinc, 0.96, Fnthgt]);
uicontrol(dlgh, ...
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold', ...
    'String', OptString, ...
    'Position', [0.35, 1-linenum*yinc, 0.63, Fnthgt]);
% label
linenum=linenum+1;
uicontrol(dlgh, ...
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'HorizontalAlignment', 'left', ...
    'String', 'Label : ', ...
    'Position', [0.02, 1-linenum*yinc, 0.96, Fnthgt]);
uicontrol(dlgh, ...
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold', ...
    'String', RecHead.Label(Rec,:), ...
    'Position', [0.35, 1-linenum*yinc, 0.63, Fnthgt]);
% Transducer
linenum=linenum+1;
uicontrol(dlgh, ...
  'Style', 'Text', ...
  'Units', 'Normalized', ...
  'HorizontalAlignment', 'left', ...
  'String', 'Transducer :', ...
  'Position', [0.02, 1-linenum*yinc, 0.96, Fnthgt]);
uicontrol(dlgh, ...
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold', ...
    'String', RecHead.Transducer(Rec, :), ...
    'Position', [0.35, 1-linenum*yinc, 0.63, Fnthgt]);
% Prefilter
linenum=linenum+1;
uicontrol(dlgh, ...
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'HorizontalAlignment', 'left', ...
    'String', 'Prefilter :', ...
    'Position', [0.02, 1-linenum*yinc, 0.96, Fnthgt]);
uicontrol(dlgh, ...
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold', ...
    'String', RecHead.PreFilt(Rec, :), ...
    'Position', [0.35, 1-linenum*yinc, 0.63, Fnthgt]);
% samples per records
linenum=linenum+1;
uicontrol(dlgh, ...
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'HorizontalAlignment', 'left', ...
    'String', 'Samples per record  :', ...
    'Position', [0.02, 1-linenum*yinc, 0.96, Fnthgt]);
uicontrol(dlgh, ...
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold', ...
    'String', sprintf('%d', RecHead.SPR), ...
    'Position', [0.35, 1-linenum*yinc, 0.63, Fnthgt]);
% physical min
linenum=linenum+1;
uicontrol(dlgh, ...
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'HorizontalAlignment', 'left', ...
    'String', 'Physical min :', ...
    'Position', [0.02, 1-linenum*yinc, 0.45, Fnthgt]);
uicontrol(dlgh, ...
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold', ...
    'String', sprintf('%f', RecHead.PhysMin(Rec)), ...
    'Position', [0.25, 1-linenum*yinc, 0.25, Fnthgt]);
% physical max
uicontrol(dlgh, ...
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'HorizontalAlignment', 'left', ...
    'String', 'Physical max :', ...
    'Position', [0.52, 1-linenum*yinc, 0.45, Fnthgt]);
uicontrol(dlgh, ...
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold', ...
    'String', sprintf('%f', RecHead.PhysMax(Rec)), ...
    'Position', [0.77, 1-linenum*yinc, 0.25, Fnthgt]);
% digital min
linenum=linenum+1;
uicontrol(dlgh, ...
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'HorizontalAlignment', 'left', ...
    'String', 'Digital min :', ...
    'Position', [0.02, 1-linenum*yinc, 0.45, Fnthgt]);
uicontrol(dlgh, ...
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold', ...
    'String', sprintf('%d', RecHead.DigMin(Rec)), ...
    'Position', [0.25, 1-linenum*yinc, 0.25, Fnthgt]);
% digital max
uicontrol(dlgh, ...
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'HorizontalAlignment', 'left', ...
    'String', 'Digital max :', ...
    'Position', [0.52, 1-linenum*yinc, 0.45, Fnthgt]);
uicontrol(dlgh, ...
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold', ...
    'String', sprintf('%d', RecHead.DigMax(Rec)), ...
    'Position', [0.77, 1-linenum*yinc, 0.25, Fnthgt]);
% display min
linenum=linenum+1;
uicontrol(dlgh, ...
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'HorizontalAlignment', 'left', ...
    'String', 'Display min :', ...
    'Position', [0.02, 1-linenum*yinc, 0.2, Fnthgt]);
dminh = uicontrol(dlgh, ...
    'Style', 'Edit', ...
    'Units', 'Normalized', ...
    'HorizontalAlignment', 'left', ...
    'BackGroundColor', [1 1 1], ...
    'String', sprintf('%f', RecDisp.DisplayMin(Rec)), ...
    'Position', [0.25, 1-linenum*yinc, 0.20, Fnthgt+0.02]);
LocalResizeUI(dminh, [NaN 1.1 0 0]);
% display max
uicontrol(dlgh, ...
    'Style', 'Text', ...
    'Units', 'Normalized', ...
    'HorizontalAlignment', 'left', ...
    'String', 'Display max :', ...
    'Position', [0.52, 1-linenum*yinc, 0.2, Fnthgt]);
dmaxh = uicontrol(dlgh, ...
    'Style', 'Edit', ...
    'Units', 'Normalized', ...
    'HorizontalAlignment', 'left', ...
    'BackGroundColor', [1 1 1], ...
    'String', sprintf('%f', RecDisp.DisplayMax(Rec)), ...
    'Position', [0.77, 1-linenum*yinc, 0.20, Fnthgt+0.02]);
LocalResizeUI(dmaxh, [NaN 1.1 0 0]);
% Plot-type
linenum=linenum+1;
ploth = uicontrol(dlgh, ...
    'Style', 'CheckBox', ...
    'Units', 'Normalized', ...
    'String', '  Stairstep-Plot', ...
    'Value', RecDisp.StairPlot(Rec), ...
    'Position', [0.02, 1-linenum*yinc, 0.7, Fnthgt]);
LocalResizeUI(ploth, [1 1 0.05 0]);
% Invert plot
updownh = uicontrol(dlgh, ...
    'Style', 'CheckBox', ...
    'Units', 'Normalized', ...
    'String', '  Plot upside down', ...
    'Value', RecDisp.UpDownPlot(Rec), ...
    'Position', [0.52, 1-linenum*yinc, 0.7, Fnthgt]);
LocalResizeUI(updownh, [1 1 0.05 0]);
% buttons
linenum=linenum+1;
uicontrol(dlgh, ...
    'Style', 'PushButton', ...
    'Units', 'Normalized', ...
    'String', 'OK', ...
    'Position', [0.05, 0.02, 0.3, 0.08], ...
    'Callback', 'set(gco,''UserData'',''OK'');uiresume;');
uicontrol(dlgh, ...
    'Style', 'PushButton', ...
    'Units', 'Normalized', ...
    'String', 'Apply to all', ...
    'Position', [0.4, 0.02, 0.2, 0.08], ...
    'Callback', 'set(gco,''UserData'',''ApplyAll'');uiresume;');
uicontrol(dlgh, ...
    'Style', 'PushButton', ...
    'Units', 'Normalized', ...
    'String', 'Cancel', ...
    'Position', [0.65, 0.02, 0.3, 0.08], ...
    'Callback', 'set(gco,''UserData'',''Cancel'');uiresume;');

drawnow;
uiwait(dlgh);
changed = 0;
if ~strcmp(get(gco,'UserData'),'Cancel')
  dmin = str2num(get(dminh, 'String'));
  dmax = str2num(get(dmaxh, 'String'));
  stplot = get(ploth, 'Value');
  udplot = get(updownh, 'Value');
  if (length(dmin) == 1) & (length(dmax) == 1) & (dmin < dmax)
    changed = 1;
    if strcmp(get(gco,'UserData'),'OK')
      switch(whichbut{1})
       case 'EDF'
        Data.Display.EDF(Num).DisplayMin(Rec) = dmin;
        Data.Display.EDF(Num).DisplayMax(Rec) = dmax;
        Data.Display.EDF(Num).StairPlot(Rec) = stplot;
        Data.Display.EDF(Num).UpDownPlot(Rec) = udplot;
       case 'PLUGIN'
        Data.Display.Plugin(Num).DisplayMin(Rec) = dmin;
        Data.Display.Plugin(Num).DisplayMax(Rec) = dmax;
        Data.Display.Plugin(Num).StairPlot(Rec) = stplot;
        Data.Display.Plugin(Num).UpDownPlot(Rec) = udplot;
      end
    else
      % set all files
      for fnum = 1:length(Data.Display.EDF)
        for rec = 1:length(Data.Display.EDF(fnum).DisplayMin)
          Data.Display.EDF(fnum).DisplayMin(rec) = dmin;
          Data.Display.EDF(fnum).DisplayMax(rec) = dmax;
          Data.Display.EDF(fnum).StairPlot(rec) = stplot;
          Data.Display.EDF(fnum).UpDownPlot(rec) = udplot;
        end
      end
      if isfield(Data.Display, 'Plugin')
        for fnum = 1:length(Data.Display.Plugin)
          for rec = 1:length(Data.Display.Plugin(fnum).DisplayMin)
            Data.Display.Plugin(fnum).DisplayMin(rec) = dmin;
            Data.Display.Plugin(fnum).DisplayMax(rec) = dmax;
            Data.Display.Plugin(fnum).StairPlot(rec) = stplot;
            Data.Display.Plugin(fnum).UpDownPlot(rec) = udplot;
          end
        end
      end
    end
  end
end
delete(dlgh);
if changed
  set(findobj('Tag', 'ViewEDFFigure'), 'UserData', Data);
  LocalRepaint;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% LocalRescale
% change range of a plot
function LocalRescale(Mode)
Data = get(findobj('Tag', 'ViewEDFFigure'), 'UserData');
whichbut = get(gcbo, 'UserData');
Rec = whichbut{2};
Num = whichbut{3};

switch(whichbut{1})
  case 'EDF'
    DisData = Data.Display.EDF(Num);
    [DisData.DisplayMin(Rec), DisData.DisplayMax(Rec)] = ...
        LocalCalcNewRange(DisData.DisplayMin(Rec), DisData.DisplayMax(Rec), ...
        Mode);
    Data.Display.EDF(Num) = DisData;
  case 'PLUGIN'
    DisData = Data.Display.Plugin(Num);
    [DisData.DisplayMin(Rec), DisData.DisplayMax(Rec)] = ...
        LocalCalcNewRange(DisData.DisplayMin(Rec), DisData.DisplayMax(Rec), ...
        Mode);
    Data.Display.Plugin(Num) = DisData;
end
set(findobj('Tag', 'ViewEDFFigure'), 'UserData', Data);
LocalRepaint;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% LocalRescaleAll
% change range of all plots
function LocalRescaleAll(Mode)
Data = get(findobj('Tag', 'ViewEDFFigure'), 'UserData');
if length(Data.EDF)==0
  return;
end
% EDF files
for i=1:length(Data.EDF)
  DisData = Data.Display.EDF(i);
  [DisData.DisplayMin, DisData.DisplayMax] = ...
      LocalCalcNewRange(DisData.DisplayMin, DisData.DisplayMax, Mode);
  Data.Display.EDF(i) = DisData;
end
% plugins
for i=1:length(Data.Plugin)
  DisData = Data.Display.Plugin(i);
  [DisData.DisplayMin, DisData.DisplayMax] = ...
      LocalCalcNewRange(DisData.DisplayMin, DisData.DisplayMax, Mode);
  Data.Display.Plugin(i) = DisData;
end
set(findobj('Tag', 'ViewEDFFigure'), 'UserData', Data);
LocalRepaint;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% LocalCalcNewRange
% calculate new range for plots
function [ResMin, ResMax] = LocalCalcNewRange(PMin, PMax, Mode)
dmean = (PMin+PMax) / 2;
switch Mode
  case 'up'
    ResMin = (PMin - dmean)*0.5 + dmean;
    ResMax = (PMax - dmean)*0.5 + dmean;
  case 'down'
    ResMin = (PMin - dmean)*2 + dmean;
    ResMax = (PMax - dmean)*2 + dmean;
  otherwise
    ResMin = PMin;
    ResMax = PMax;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% LocalToggleUpdatePlugin
% Toggle flag for updating plugin-results 
function LocalToggleUpdatePlugin()
Data = get(findobj('Tag', 'ViewEDFFigure'), 'UserData');
Data.UpdatePlugin = ~Data.UpdatePlugin;
if Data.UpdatePlugin
  set(gcbo, 'Checked', 'on');
else
  set(gcbo, 'Checked', 'off');
end
set(findobj('Tag', 'ViewEDFFigure'), 'UserData', Data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% LocalToggleShowRange
% Toggle flag for showing plot ranges
function LocalToggleShowRange()
Data = get(findobj('Tag', 'ViewEDFFigure'), 'UserData');
Data.ShowRange = ~Data.ShowRange;
if Data.ShowRange
  set(gcbo, 'Checked', 'on');
else
  set(gcbo, 'Checked', 'off');
end
set(findobj('Tag', 'ViewEDFFigure'), 'UserData', Data);
LocalRepaint;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% LocalToggleShowCursor
% Toggle flag for showing cursor 
function LocalToggleShowCursor()
figure(findobj('Tag', 'ViewEDFFigure'));
Data = get(gcf, 'UserData');
Data.ShowCursor = ~Data.ShowCursor;
if Data.ShowCursor
  set(gcbo, 'Checked', 'on');
  set(gcf, 'WindowButtonDownFcn', 'viewedf MoveCursor Down');
  Data.Display.Cursor = LocalDrawCursor(Data.Display);
else
  delete(Data.Display.Cursor.Line(:));
  delete(Data.Display.Cursor.Menu.Text);
  delete(Data.Display.Cursor.Menu.Menu);
  Data.Display.Cursor = [];
  set(gcbo, 'Checked', 'off');
  set(gcf, 'WindowButtonMotionFcn', '', ...
           'WindowButtonDownFcn', '', ...
           'WindowButtonUpFcn', '');
end
set(gcf, 'UserData', Data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% LocalDrawCursor
% Show cursor
function Cursor = LocalDrawCursor(DisplayData, NewPos, TimeStr)
if nargin ~= 3
  % find good placement for lines
  Cursor.Pos = 0.5;
  Cursor.Menu.Menu = uicontextmenu;
  Cursor.Menu.Text  = uimenu(Cursor.Menu.Menu, ...
                             'Label', '', ...
                             'Callback', '');
  % generate lines
  for i=1:length(DisplayData.Axes)
    Cursor.AxesLim(i, :) = get(DisplayData.Axes(i).Plot, 'XLim');
    XPos = Cursor.AxesLim(i,1) + (Cursor.AxesLim(i,2) - Cursor.AxesLim(i,1)) ...
           * Cursor.Pos;
    Cursor.Line(i) = line('Parent', DisplayData.Axes(i).Plot, ...
                          'XData', [XPos XPos], ...
                          'YData', get(DisplayData.Axes(i).Plot, 'YLim'), ...
                          'EraseMode', 'xor', ...
                          'UIContextMenu', Cursor.Menu.Menu, ...
                          'Color', [1 0.2 0.2]);
    Cursor.AxesPos(i, :) = get(DisplayData.Axes(i).Plot, 'Position');
  end
else
  Cursor = DisplayData.Cursor;
  Cursor.Pos = NewPos;
  for i=1:length(DisplayData.Axes)
    XPos = Cursor.AxesLim(i,1) + (Cursor.AxesLim(i,2) - Cursor.AxesLim(i,1)) ...
           * Cursor.Pos;
    set(Cursor.Line(i), 'XData', [XPos XPos]);
  end
  set(Cursor.Menu.Text, 'Label', TimeStr);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% LocalMoveCursor
% Move cursor
function LocalMoveCursor(Option)
Data = get(gcf, 'UserData');
CurrPt = get(gcf, 'CurrentPoint');
% check if click was in axis
AxPos = Data.Display.Cursor.AxesPos;
temp = CurrPt(1) > AxPos(:,1) & CurrPt(1) < (AxPos(:,1) + AxPos(:,3)) & ...
       CurrPt(2) > AxPos(:,2) & CurrPt(2) < (AxPos(:,2) + AxPos(:,4));
if ~any(temp)
  return;
end

% find axis under pointer and calculate position
ax = find(temp);
Pos = get(Data.Display.Axes(ax).Plot, 'CurrentPoint');
Pos = Pos(2,1);
Pos = (Pos - Data.Display.Cursor.AxesLim(ax,1)) / (Data.Display.Cursor.AxesLim(ax,2) ...
                                                  - Data.Display.Cursor.AxesLim(ax,1));

tsec = Data.EDF(1).Head.Dur * Data.Display.EDF(1).ShowRecords(1) + ...
       Data.EDF(1).Head.Dur * Data.Display.EDF(1).ShowRecords(2) * Pos;
tmin = floor(tsec/60);
tsec = rem(tsec,60);
th = floor(tmin/60);
tmin = rem(tmin, 60);
TimeStr = sprintf('%02d:%02d:%02.1f', th,tmin,tsec);

switch Option
 case 'Down'
  Data.Display.Cursor = LocalDrawCursor(Data.Display, Pos, TimeStr);
  set(gcf, 'WindowButtonMotionFcn', 'viewedf MoveCursor Move', ...
           'WindowButtonUpFcn', 'viewedf MoveCursor Up');
 case 'Move'
  Data.Display.Cursor = LocalDrawCursor(Data.Display, Pos, TimeStr);
 case 'Up'
  set(gcf, 'WindowButtonMotionFcn', '');
end
set(gcf, 'UserData', Data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LocalResizeUI
% resize UI control according to label size
function uihandle = LocalResizeUI(uihandle, options)
if nargin < 2
  options = [1 1 0 0];
end
ext = get(uihandle, 'Extent');
pos = get(uihandle, 'Position');
if ~isnan(options(1))
  pos(3) = ext(3)*options(1) + options(3);
end
if ~isnan(options(2))
  pos(4) = ext(4)*options(2) + options(4);
end
set(uihandle, 'Position', pos);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% LocalNumRecords
% select the number of records to be shown simultaniously
function LocalNumRecords()
Data = get(findobj('Tag', 'ViewEDFFigure'), 'UserData');
if length(Data.EDF) == 0
  errordlg('No EDF-File is open!', 'Error');
  return;
end

Dur = LocalGetEDFInfo('Dur', Data.EDF);
[MaxDur, MDInd] = max(Dur);

answer = inputdlg({'Select number of records'}, ...
    'Number of records on screen', 1, ...
    {int2str(Data.Display.EDF(MDInd).ShowRecords(2))});
if (length(answer) ~= 0)
  answer = str2num(answer{1});
  if ~isempty(answer)
    for i=1:length(Data.EDF)
      Data.Display.EDF(i).ShowRecords(2) = answer * MaxDur / Dur(i);  
      [Data.EDF(i).Record, Data.EDF(i).Head] = ...
          LocalEDFRead(Data.EDF(i).Head, Data.Display.EDF(i).ShowRecords); 
    end
    set(findobj('Tag', 'ViewEDFFigure'), 'UserData', Data);
    LocalRepaint(0);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% LocalFileInfo
% display information about current file
function LocalFileInfo(Update)
if nargin == 0
  Data=get(findobj('Tag', 'ViewEDFFigure'), 'UserData');
  if length(Data.EDF)==0
    errordlg('No EDF-File is open!', 'Error');
    return;
  end
  dlgh = dialog(...
      'Name', 'EDF-File information', ...
      'CloseRequestFcn', 'set(gcf,''UserData'',''Cancel'');uiresume;');
  dlgpos = get(dlgh, 'Position');
  set(dlgh, 'Position', [dlgpos(1),dlgpos(2),450,250]);
  fnthght = LocalGetFontHeight;

  Local.EDFNames = {};
  % get filenames
  for i = 1:length(Data.EDF)
    Local.EDFNames = {Local.EDFNames{:}, Data.EDF(i).Head.FileName};
    Local.EDFInfo{i,1} = Data.EDF(i).Head.VERSION;
    Local.EDFInfo{i,2} = Data.EDF(i).Head.PID;
    Local.EDFInfo{i,3} = Data.EDF(i).Head.RID;
    Local.EDFInfo{i,4} = sprintf('%02d/%02d/%02d', Data.EDF(i).Head.T0([3 2 1]));
    Local.EDFInfo{i,5} = sprintf('%02d:%02d:%02d', Data.EDF(i).Head.T0([4 5 6]));
    Local.EDFInfo{i,6} = sprintf('%d', Data.EDF(i).Head.NS);
    Local.EDFInfo{i,7} = sprintf('%d', Data.EDF(i).Head.NRec);
  end
  
  Local.poph = uicontrol(dlgh, ...
      'Style', 'Popup', ...
      'Units', 'Normalized', ...
      'String', Local.EDFNames, ...
      'Value', 1, ...
      'Position', [0.05, 0.9, 0.9, 0.1], ...
      'Callback', 'viewedf FileInfo Update');
  LocalResizeUI(Local.poph, [NaN 1.1 0 0]);
  %HeaderVersion
  num = 1;
  yinc = (1 - 4*fnthght) / 8;
  ystart = 0.9;
  uicontrol(dlgh, ...
      'Style', 'Text', ...
      'HorizontalAlignment', 'left', ...
      'Units', 'Normalized', ...
      'String', 'HeaderVersion : ', ...
      'Position', [0.05, ystart - num*yinc, 0.3, fnthght]);
  Local.txth(num) = uicontrol(dlgh, ...
      'Style', 'Text', ...
      'HorizontalAlignment', 'left', ...
      'Units', 'Normalized', ...
      'FontWeight', 'bold', ...
      'String', '', ...
      'Position', [0.45, ystart - num*yinc, 0.5, fnthght]);
  %Patient ID
  num = num + 1;
  uicontrol(dlgh, ...
      'Style', 'Text', ...
      'HorizontalAlignment', 'left', ...
      'Units', 'Normalized', ...
      'String', 'Patient ID : ', ...
      'Position', [0.05, ystart - num*yinc, 0.3, fnthght]);
  Local.txth(num) = uicontrol(dlgh, ...
      'Style', 'Text', ...
      'HorizontalAlignment', 'left', ...
      'Units', 'Normalized', ...
      'FontWeight', 'bold', ...
      'String', '', ...
      'Position', [0.45, ystart - num*yinc, 0.5, fnthght]);
  %Recording ID
  num = num + 1;
  uicontrol(dlgh, ...
      'Style', 'Text', ...
      'HorizontalAlignment', 'left', ...
      'Units', 'Normalized', ...
      'String', 'Recording ID : ', ...
      'Position', [0.05, ystart - num*yinc, 0.3, fnthght]);
  Local.txth(num) = uicontrol(dlgh, ...
      'Style', 'Text', ...
      'HorizontalAlignment', 'left', ...
      'Units', 'Normalized', ...
      'FontWeight', 'bold', ...
      'String', '', ...
      'Position', [0.45, ystart - num*yinc, 0.5, fnthght]);
  % startdate
  num = num + 1;
  uicontrol(dlgh, ...
      'Style', 'Text', ...
      'HorizontalAlignment', 'left', ...
      'Units', 'Normalized', ...
      'String', 'Start date : ', ...
      'Position', [0.05, ystart - num*yinc, 0.3, fnthght]);
  Local.txth(num) = uicontrol(dlgh, ...
      'Style', 'Text', ...
      'HorizontalAlignment', 'left', ...
      'Units', 'Normalized', ...
      'FontWeight', 'bold', ...
      'String', '', ...
      'Position', [0.45, ystart - num*yinc, 0.5, fnthght]);
  % start time
  num = num + 1;
  uicontrol(dlgh, ...
      'Style', 'Text', ...
      'HorizontalAlignment', 'left', ...
      'Units', 'Normalized', ...
      'String', 'Start time : ', ...
      'Position', [0.05, ystart - num*yinc, 0.3, fnthght]);
  Local.txth(num) = uicontrol(dlgh, ...
      'Style', 'Text', ...
      'HorizontalAlignment', 'left', ...
      'Units', 'Normalized', ...
      'FontWeight', 'bold', ...
      'String', '', ...
      'Position', [0.45, ystart - num*yinc, 0.5, fnthght]);
  % Num channels
  num = num + 1;
  uicontrol(dlgh, ...
      'Style', 'Text', ...
      'HorizontalAlignment', 'left', ...
      'Units', 'Normalized', ...
      'String', 'Number of channels : ', ...
      'Position', [0.05, ystart - num*yinc, 0.3, fnthght]);
  Local.txth(num) = uicontrol(dlgh, ...
      'Style', 'Text', ...
      'HorizontalAlignment', 'left', ...
      'Units', 'Normalized', ...
      'FontWeight', 'bold', ...
      'String', '', ...
      'Position', [0.45, ystart - num*yinc, 0.5, fnthght]);
  % num records
  num = num + 1;
  uicontrol(dlgh, ...
      'Style', 'Text', ...
      'HorizontalAlignment', 'left', ...
      'Units', 'Normalized', ...
      'String', 'Number of records : ', ...
      'Position', [0.05, ystart - num*yinc, 0.3, fnthght]);
  Local.txth(num) = uicontrol(dlgh, ...
      'Style', 'Text', ...
      'HorizontalAlignment', 'left', ...
      'Units', 'Normalized', ...
      'FontWeight', 'bold', ...
      'String', '', ...
      'Position', [0.45, ystart - num*yinc, 0.5, fnthght]);
  % buttons
  uicontrol(dlgh, ...
      'Style', 'PushButton', ...
      'Units', 'Normalized', ...
      'String', 'OK', ...
      'Position', [0.3, 0.02, 0.4, 0.1], ...
      'Callback', 'uiresume;');
  Data=[];
  set(dlgh, 'UserData', Local);
  LocalFileInfo('Update');
  
  drawnow;
  uiwait(dlgh);
  delete(dlgh);
else
  Data=get(gcf,'UserData');
  which = get(Data.poph, 'Value');
  for i = 1:7
    set(Data.txth(i), 'String', Data.EDFInfo{which,i});
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LocalSelectChannels
% select channels to be displayed on the screen
function LocalSelectChannels(Parameter)
if nargin == 0
  Data=get(findobj('Tag', 'ViewEDFFigure'), 'UserData');
  if length(Data.EDF)==0
    errordlg('No EDF-File is open!', 'Error');
    return;
  end
  dlgh = dialog(...
      'Name', 'EDF-File information', ...
      'CloseRequestFcn', 'uiresume;');
  
  Local.Names = {};
  Local.Selected = {};
  Local.Label = {};
  % get filenames
  for i = 1:length(Data.EDF)
    Local.Names = {Local.Names{:}, ['EDF: ', Data.EDF(i).Head.FileName]};
    temp = zeros(Data.EDF(i).Head.NS,1);
    temp(Data.Display.EDF(i).ShowSignals) = 1;
    Local.Selected = {Local.Selected{:}, temp};
    Local.Label = {Local.Label{:}, Data.EDF(i).Head.Label};
  end
  %get plugin names
  for i = 1:length(Data.Plugin)
    Local.Names = {Local.Names{:}, ['Plugin: ', Data.Plugin(i).Label, ...
            ' (File ', upper(Data.Plugin(i).Name), ')']};
    temp = zeros(Data.Plugin(i).EDF.Head.NS,1);
    temp(Data.Display.Plugin(i).ShowSignals) = 1;
    Local.Selected = {Local.Selected{:}, temp};
    Local.Label = {Local.Label{:}, Data.Plugin(i).EDF.Head.Label};
  end

  uicontrol(dlgh, ...
      'Style', 'Frame', ...
      'Units', 'Normalized', ...
      'Position', [0.02, 0.12, 0.96, 0.76]);
  Local.Poph = uicontrol(dlgh, ...
      'Style', 'Popup', ...
      'Units', 'Normalized', ...
      'String', Local.Names, ...
      'Value', 1, ...
      'Position', [0.05, 0.9, 0.9, 0.1], ...
      'Callback', 'viewedf Channels Update');
  LocalResizeUI(Local.Poph, [NaN 1.1 0 0]);
  % buttons
  uicontrol(dlgh, ...
      'Style', 'PushButton', ...
      'Units', 'Normalized', ...
      'String', 'Select', ...
      'Position', [0.1, 0.14, 0.3, 0.08], ...
      'Callback', 'viewedf Channels Select');
  uicontrol(dlgh, ...
      'Style', 'PushButton', ...
      'Units', 'Normalized', ...
      'String', 'Unselect', ...
      'Position', [0.6, 0.14, 0.3, 0.08], ...
      'Callback', 'viewedf Channels Unselect');
  uicontrol(dlgh, ...
      'Style', 'PushButton', ...
      'Units', 'Normalized', ...
      'String', 'OK', ...
      'Position', [0.1, 0.02, 0.3, 0.08], ...
      'Callback', 'set(gco,''UserData'',''OK'');uiresume;');
  uicontrol(dlgh, ...
      'Style', 'PushButton', ...
      'Units', 'Normalized', ...
      'String', 'Cancel', ...
      'Position', [0.6, 0.02, 0.3, 0.08], ...
      'Callback', 'set(gco,''UserData'',''Cancel'');uiresume;');
  Local.Checkh = [];
  set(dlgh, 'UserData', Local);
  LocalSelectChannels('Update'); % draw information
  
  drawnow;
  uiwait(dlgh);
  changed = 0;
  Local = get(dlgh, 'UserData');
  temp = get(Local.Checkh, 'Value');
  Local.Selected{Local.OldWhich} = cat(1,temp{:});
  if strcmp(get(gco,'UserData'),'OK')
    changed = 1;
    for i = 1:length(Data.EDF)
      Data.Display.EDF(i).ShowSignals = find(Local.Selected{i});
    end
    for i = 1:length(Data.Plugin)
      Data.Display.Plugin(i).ShowSignals = find(Local.Selected{i+length(Data.EDF)});
    end
  end
  delete(dlgh);
  if changed
    set(findobj('Tag', 'ViewEDFFigure'), 'UserData', Data);
    LocalRepaint(0);
  end
else
  Data=get(gcf,'UserData');
  which = get(Data.Poph, 'Value');
  switch (Parameter)
    case 'Update'
      if ~isempty(Data.Checkh)
        temp = get(Data.Checkh, 'Value');
        if iscell(temp)
          Data.Selected{Data.OldWhich} = cat(1,temp{:});
        else
          Data.Selected{Data.OldWhich} = temp;
        end  
        delete(Data.Checkh);
        Data.Checkh = [];
      end
      Data.OldWhich = which;
    
      cbx = 0.1;
      cby = 0.82;
      cbxspace = 0.4;
      cbyspace = 0.6 / ceil(length(Data.Selected{which}) / 2) + 0.005; 
      cbwidth = 0.35;
      cbheight = LocalGetFontHeight;
      for i = 1:length(Data.Selected{which})
        Data.Checkh(i) = uicontrol(gcf, ...
            'Style', 'CheckBox', ...
            'Units', 'Normalized', ...
            'String', Data.Label{which}(i, :), ...
            'Value', Data.Selected{which}(i), ...
            'Position', [cbx + rem(i-1,2)*cbxspace, cby - floor((i-1)/2)*cbyspace, ...
              cbwidth, cbheight]);
      end
    case 'Select'
      Data.Selected{which}(:) = 1;
      set(Data.Checkh(:), 'Value', 1);
    case 'Unselect'
      Data.Selected{which}(:) = 0;
      set(Data.Checkh(:), 'Value', 0);
  end
  set(gcf, 'UserData', Data);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LocalKeyPress
% handle keystrokes
function LocalKeyPress()
switch upper(get(gcf, 'CurrentCharacter'))
  case '+'
    feval('viewedf', 'Next');
  case '-'
    feval('viewedf', 'Prev');
  case 'U'
    LocalRescaleAll('up');
  case 'D'
    LocalRescaleAll('down');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LocalAbout
% display program info
function LocalAbout()
helpdlg(sprintf([ 'EDF (European-Data-Format) file-viewer.\n', ...
      'Version 3.04Alpha\n\n' ...
      '(c) 1998-2001 Herbert Ramoser\n', ...
      '     herbert.ramoser@arcs.ac.at\n\n' ...
      'Comments or suggestions may be sent to the author.\n\n', ...
      'This Software is subject to the GNU public license.']), ...
    'About VIEWEDF');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LocalHelp
% display help page
function LocalHelp(Par)
name = which('viewedf');
name = [ 'file://', name(1:max(find(name == 'm'))-1), 'html'];
web(name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LocalPrint
% display print dialog
function LocalPrint()
printdlg('-crossplatform', findobj('Tag', 'ViewEDFFigure'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LocalWatchOn
% display watch pointer
function LocalWatchOn()
set(findobj('Tag', 'ViewEDFFigure'), 'Pointer', 'watch');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LocalWatchOff
% display arrow-pointer
function LocalWatchOff()
set(findobj('Tag', 'ViewEDFFigure'), 'Pointer', 'arrow');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LocalGetFontHeight
% get Fontheight
function hght=LocalGetFontHeight()
tempH = uicontrol(...
    'Style', 'Text', ...
    'String', 'Gg', ...
    'Units', 'Normalized', ...
    'FontUnits', 'Normalized', ...
    'Position', [0, 0, 1, 1], ...
    'Visible', 'off');
hght = get(tempH, 'FontSize') * 1.25; % 1.25 makes things look better
delete(tempH);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LocalGetFontWidth
% get Fontwidth
function width=LocalGetFontWidth()
tempH = uicontrol(...
    'Style', 'Text', ...
    'String', 'X', ...
    'Units', 'Normalized', ...
    'FontUnits', 'Normalized', ...
    'Position', [0, 0, 1, 1], ...
    'Visible', 'off');
width = get(tempH, 'Extent'); 
width = width(3);
delete(tempH);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LocalEDFOpen
% open a EDF (or GDF) file
function LocalEDFOpen(Filename)
% show EDF header errors
ShowHeadErr = 0;

if nargin == 0
  % get filename	
  [edfname,edfpath]=uigetfile('*.*','Open EDF File');
  Filename = [edfpath,edfname];
  if edfname == 0
    return;
  end
end 
   
if 0, 
%%% becomes obsolete 
% variables to find things in the header
H1idx=[8 80 80 8 8 8 44 8 8 4];
H2idx=[16 80 8 8 8 8 8 80 8 32];
GDFTYP_BYTE=[1 1 1 2 2 4 4 8 8 4 8 0 0 0 0 0 4 8]';
GDFTYPES=[0 1 2 3 4 5 6 7 16 17];

fid=fopen(Filename, 'r', 'ieee-le'); 
if (fid < 0) 
  errordlg('Error reading file', 'File Error');
  return
end;
Data = get(findobj('Tag', 'ViewEDFFigure'), 'UserData');
numedf = length(Data.EDF) + 1;

% read all data
EDFHead.FILE.FID = fid;
EDFHead.FILE.OPEN = 1;
EDFHead.FileName = Filename;

PPos = min([max(find(Filename == '.')) length(Filename) + 1]);
SPos = max([0 find(Filename == filesep)]);
EDFHead.FILE.Ext = Filename(PPos+1:length(Filename));
EDFHead.FILE.Name = Filename(SPos+1:PPos-1);
EDFHead.FILE.Path = edfpath;

H1 = setstr(fread(EDFHead.FILE.FID,184,'uchar')');     
EDFHead.VERSION = H1(1:8);         % 8 Byte  Versionsnummer 
IsGDF = strcmp(EDFHead.VERSION(1:3), 'GDF');
if (~strcmp(EDFHead.VERSION, '0       ') & ~IsGDF)
  errordlg('Unknown file version', 'File error');
  return;
end
EDFHead.PID = deblank(H1(9:88));   % 80 Byte local patient identification
EDFHead.RID = deblank(H1(89:168)); % 80 Byte local recording identification
if IsGDF                           % handle different file formats
  EDFHead.T0 = [str2num(H1(168 + [1:4])) ...
        str2num(H1(168 + [5 6])) ...
        str2num(H1(168 + [7 8])) ...
        str2num(H1(168 + [9 10])) ...
        str2num(H1(168 + [11 12])) ...
        str2num(H1(168 + [13:16]))];
  if str2num(EDFHead.VERSION(4:8)) < 0.12
    tmp = setstr(fread(EDFHead.FILE.FID, 8, 'uchar')');  % header-length
    EDFHead.HeadLen = str2num(tmp);
  else
    EDFHead.HeadLen = fread(EDFHead.FILE.FID, 1, 'int64');
  end
  tmp = fread(EDFHead.FILE.FID, 44, 'uchar'); % 44 bytes reserved
  EDFHead.NRec = fread(EDFHead.FILE.FID, 1, 'int64');
  if strcmp(EDFHead.VERSION(4:8),' 0.10')
    EDFHead.Dur =  fread(EDFHead.FILE.FID, 1, 'float64');   
  else
    tmp = fread(EDFHead.FILE.FID, 2, 'uint32');    
    EDFHead.Dur =  tmp(1)./tmp(2);
  end
  EDFHead.NS   = fread(EDFHead.FILE.FID, 1, 'uint32');
else
  EDFHead.T0 = [str2num(H1(168+[7 8])) ...
        str2num(H1(168+[4 5])) ...
        str2num(H1(168+[1 2])) ...
        str2num(H1(168+[9 10])) ...
        str2num(H1(168+[12 13])) ...
        str2num(H1(168+[15 16])) ]; 
  if EDFHead.T0(1) < 91
    EDFHead.T0(1) = 2000 + EDFHead.T0(1);
  elseif EDFHead.T0(1) < 100
    EDFHead.T0(1) = 1900 + EDFHead.T0(1);
  end
  H1(185:256) = setstr(fread(EDFHead.FILE.FID, 256-184, 'uchar')');
  EDFHead.HeadLen = str2num(H1(185:192));  % 8 Byte  Length of Header
  EDFHead.NRec = str2num(H1(237:244));     % 8 Byte  # of data records
  EDFHead.Dur = str2num(H1(245:252));      % 8 Byte  # duration of data record in sec
  EDFHead.NS = str2num(H1(253:256));       % 8 Byte  # of signals
end

if ~IsGDF
  idx1 = cumsum([0 H2idx]);
  idx2 = EDFHead.NS * idx1;
  h2 = zeros(EDFHead.NS, 256);
  H2 = fread(EDFHead.FILE.FID, EDFHead.NS * 256, 'uchar');
  H2(H2==0) = 32; % set zero padded strings to blanks
  for k = 1:length(H2idx)
    h2(:, (idx1(k)+1):idx1(k+1)) = reshape(H2((idx2(k)+1):idx2(k+1)), ...
        H2idx(k), EDFHead.NS)';
  end
  h2 = setstr(h2);
  EDFHead.Label      = h2(:, idx1(1)+1:idx1(2));
  EDFHead.Transducer = h2(:, idx1(2)+1:idx1(3));
  EDFHead.PhysDim    = cellstr(h2(:, idx1(3)+1:idx1(4)));
  EDFHead.PhysMin = str2num(h2(:, idx1(4)+1:idx1(5)));
  EDFHead.PhysMax = str2num(h2(:, idx1(5)+1:idx1(6)));
  EDFHead.DigMin  = str2num(h2(:, idx1(6)+1:idx1(7)));
  EDFHead.DigMax  = str2num(h2(:, idx1(7)+1:idx1(8)));
  EDFHead.PreFilt = h2(:, idx1(8)+1:idx1(9));
  EDFHead.SPR     = str2num(h2(:, idx1(9)+1:idx1(10)));
  EDFHead.GDFTYP  = 3*ones(1, EDFHead.NS);
else
  fseek(EDFHead.FILE.FID, 256, 'bof');
  EDFHead.Label      =  setstr(fread(EDFHead.FILE.FID, [16,EDFHead.NS], ...
      'char')'); 
  EDFHead.Transducer =  setstr(fread(EDFHead.FILE.FID, [80,EDFHead.NS], ...
      'char')');
  EDFHead.PhysDim    =  cellstr(setstr(fread(EDFHead.FILE.FID, [8,EDFHead.NS], ...
      'uchar')'));
  EDFHead.PhysMin    =  fread(EDFHead.FILE.FID, [EDFHead.NS,1], 'float64');
  EDFHead.PhysMax    =  fread(EDFHead.FILE.FID, [EDFHead.NS,1], 'float64');
  EDFHead.DigMin     =  fread(EDFHead.FILE.FID, [EDFHead.NS,1], 'int64');
  EDFHead.DigMax     =  fread(EDFHead.FILE.FID, [EDFHead.NS,1], 'int64');
  EDFHead.PreFilt    =  setstr(fread(EDFHead.FILE.FID, [80,EDFHead.NS], ...
      'char')');
  EDFHead.SPR        =  fread(EDFHead.FILE.FID, [1,EDFHead.NS], 'uint32')';
  EDFHead.GDFTYP     =  fread(EDFHead.FILE.FID, [1,EDFHead.NS], 'uint32');
  tmp = (EDFHead.GDFTYP == 0);
  EDFHead.PhysMax(tmp) = 1;
  EDFHead.PhysMin(tmp) = 0;
  EDFHead.DigMax(tmp)  = 1;
  EDFHead.DigMin(tmp)  = 0;
end
 
% check validity of DigMin and DigMax
if (length(EDFHead.DigMin) ~= EDFHead.NS)
  if ShowHeadErr
    waitfor(warndlg('Failing Digital Minimum', 'Open EDF file'));
  end
  EDFHead.DigMin = -(2^15) * ones(EDFHead.NS, 1);
end
if (length(EDFHead.DigMax) ~= EDFHead.NS)
  if ShowHeadErr
    waitfor(warndlg('Failing Digital Maximum', 'Open EDF file'));
  end
  EDFHead.DigMax = (2^15-1) * ones(EDFHead.NS, 1);
end
if (any(EDFHead.DigMin >= EDFHead.DigMax))
  if ShowHeadErr
    waitfor(warndlg('Digital Minimum larger than Maximum', ['Open EDF' ...
                    ' file']));
  end
end  

% check validity of PhysMin and PhysMax
if (length(EDFHead.PhysMin) ~= EDFHead.NS)
  if ShowHeadErr
    waitfor(warndlg('EDFOPEN: Failing Physical Minimum', 'Open EDF file'));
  end
  EDFHead.PhysMin = EDFHead.DigMin;
end
if (length(EDFHead.PhysMax) ~= EDFHead.NS)
  if ShowHeadErr
    waitfor(warndlg('Warning EDFOPEN: Failing Physical Maximum', ['Open EDF' ...
                    ' file']));
  end
  EDFHead.PhysMax = EDFHead.DigMax;
end
if (any(EDFHead.PhysMin >= EDFHead.PhysMax))
  if ShowHeadErr
    waitfor(warndlg('EDFOPEN: Physical Minimum larger than Maximum', ['Open' ...
                    ' EDF file']));
  end
  EDFHead.PhysMin = EDFHead.DigMin;
  EDFHead.PhysMax = EDFHead.DigMax;
end  

EDFHead.Cal = (EDFHead.PhysMax - EDFHead.PhysMin) ./ (EDFHead.DigMax - ...
    EDFHead.DigMin);
EDFHead.Off = EDFHead.PhysMin - EDFHead.Cal .* EDFHead.DigMin;
tmp = (EDFHead.Cal < 0);
EDFHead.Cal(tmp) = 1;
EDFHead.Off(tmp) = 0;

EDFHead.Calib = [EDFHead.Off'; diag(EDFHead.Cal)];
EDFHead.SampleRate = EDFHead.SPR / EDFHead.Dur;

bi = [0; cumsum(EDFHead.SPR)]; 
EDFHead.AS.spb = sum(EDFHead.SPR);
EDFHead.AS.bi = bi;
EDFHead.AS.bpb = sum(EDFHead.SPR .* GDFTYP_BYTE(EDFHead.GDFTYP+1));
EDFHead.AS.GDFbi = [0; cumsum(EDFHead.AS.bpb)];

EDFHead.Chan_Select = (EDFHead.SPR == max(EDFHead.SPR));
for k = 1:EDFHead.NS
  if EDFHead.Chan_Select(k)
    EDFHead.ChanTyp(k) = 'N';
  else
    EDFHead.ChanTyp(k) = ' ';
  end;         
  if findstr(upper(EDFHead.Label(k,:)), 'ECG')
    EDFHead.ChanTyp(k) = 'C';
  elseif findstr(upper(EDFHead.Label(k,:)), 'EKG')
    EDFHead.ChanTyp(k) = 'C';
  elseif findstr(upper(EDFHead.Label(k,:)), 'EEG')
    EDFHead.ChanTyp(k) = 'E';
  elseif findstr(upper(EDFHead.Label(k,:)), 'EOG')
    EDFHead.ChanTyp(k) = 'O';
  elseif findstr(upper(EDFHead.Label(k,:)), 'EMG')
    EDFHead.ChanTyp(k) = 'M';
  elseif findstr(upper(EDFHead.Label(k,:)), 'RESP')
    EDFHead.ChanTyp(k) = 'R';
  end
end

fseek(EDFHead.FILE.FID, 32 * EDFHead.NS, 0);
if EDFHead.NRec == -1   % unknown record size, determine correct NRec
  fseek(EDFHead.FILE.FID, 0, 'eof');
  endpos = ftell(EDFHead.FILE.FID);
  EDFHead.NRec = floor((endpos - EDFHead.HeadLen) / (sum(EDFHead.AS.bpb))); 
end
fseek(EDFHead.FILE.FID, EDFHead.HeadLen, 'bof');
else 
	EDFHead = sopen(Filename,'r',0,'OVERFLOWDETECTION:OFF FORCEALLCHANNEL');
	EDFHead.Dur = EDFHead.SPR/EDFHead.SampleRate;
	EDFHead.PhysDim = physicalunits(EDFHead.PhysDimCode);
	
	Data   = get(findobj('Tag', 'ViewEDFFigure'), 'UserData');
	numedf = length(Data.EDF) + 1;
end; 	

Data.EDF(numedf).Head = EDFHead;  
drawnow;
set(findobj('Tag', 'ViewEDFFigure'), 'UserData', Data);
LocalResetDisplay;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LocalEDFRead
% read data from EDF file
function [Record, EDFHead] = LocalEDFRead(EDFHead, recinfo)

startrec = recinfo(1);
numrec = recinfo(2);


if 0, 
%%% obsolete ??? 
 
% define GDF data types
GDF_STRING = {'uchar', 'int8', 'uint8', 'int16', 'uint16', 'int32', ...
      'uint32', 'int64', 'uint64', '', '', '', '', '', '', '', 'float32', ...
      'float64'};

GDF_STRING{256+24} = 'bit24';
GDF_STRING{512+24} = 'ubit24';

if LocalEDFSeek(EDFHead, startrec, 'bof') < 0
  if numrec ~= 0
    waitfor(warndlg('Can not seek start record', 'Read EDF data'));
  end
  Record = {};
  return;
end

LocalWatchOn;

RecLen = max(EDFHead.SPR);
for rec = 1:numrec
  for ch = 1:EDFHead.NS
    [d, count] = fread(EDFHead.FILE.FID, EDFHead.AS.SPR(ch), ...
        GDF_STRING{EDFHead.GDFTYP(ch)+1});
    if count < 1
      break;
    end
    S(EDFHead.AS.bi(ch)+1:EDFHead.AS.bi(ch+1), rec) = d;
  end
  if count < 1
    numrec = rec;
    waitfor(warndlg(sprintf('Can not read %i records', numrec), 'Read EDF data')); 
    break
  end
end
    
for ch = 1:EDFHead.NS
  S(EDFHead.AS.bi(ch)+1 : EDFHead.AS.bi(ch+1),:) = S(EDFHead.AS.bi(ch)+1 : ...
      EDFHead.AS.bi(ch+1),:) * EDFHead.Cal(ch) + EDFHead.Off(ch);
end

for ch = 1:EDFHead.NS,
  Record{ch} = reshape(S(EDFHead.AS.bi(ch)+1 : EDFHead.AS.bi(ch+1),:), ...
      EDFHead.AS.SPR(ch) * numrec, 1);
end
EDFHead.AS.numrec = numrec;
EDFHead.AS.startrec = startrec;
LocalWatchOff;

else
LocalWatchOn;

	Dur = EDFHead.SPR/EDFHead.SampleRate;
	[s,EDFHead] = sread(EDFHead,numrec*Dur,startrec*Dur);
	for ch = 1:EDFHead.NS,
		Record{ch} = s(:,ch);
	end;

LocalWatchOff;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LocalEDFSeek
% change position in EDF file
function [status] = LocalEDFSeek(EDFHead, offset, origin)
if origin == -1 | origin == 'bof'
  offset = EDFHead.HeadLen + EDFHead.AS.bpb * offset;
  status = fseek(EDFHead.FILE.FID, offset, origin);
elseif origin == 0 | origin == 'cof'
  offset = EDFHead.AS.bpb * offset;
  status = fseek(EDFHead.FILE.FID, offset, origin);
elseif origin == 1 | origin == 'eof'
  offset = EDFHead.HeadLen + EDFHead.AS.bpb * (EDFHead.NRec + offset);
  status = fseek(EDFHead.FILE.FID, offset, -1);
end;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LocalEDFClose
% close an EDF file
function LocalEDFClose()
Data = get(findobj('Tag', 'ViewEDFFigure'), 'UserData');
if length(Data.EDF)==0
  errordlg('No EDF-File is open!', 'Error');
  return;
end

[ind, cancelled] = LocalSelectEDFFile({'Close EDF file', 'Close file'}, Data.EDF);
if cancelled
  return
end
  
%fclose(Data.EDF(ind).Head.FILE.FID);
sclose(Data.EDF(ind).Head);

for i = ind:(length(Data.EDF)-1)
  Data.Display.EDF(i) = Data.Display.EDF(i+1);
  Data.EDF(i) = Data.EDF(i+1);
end
numedf = length(Data.EDF)-1;
if numedf == 0
  Data.EDF = [];
  Data.Display.EDF = [];
else
  Data.EDF = Data.EDF(1:numedf);
  Data.Display.EDF = Data.Display.EDF(1:numedf);
end
set(findobj('Tag', 'ViewEDFFigure'), 'UserData', Data);
% take care of plugins
if (length(Data.Plugin) > 0)
  waitfor(warndlg('All plugins using the closed EDF file will be removed', ...
      'Close warning'));
  for i = length(Data.Plugin):-1:1
    if (Data.Plugin(i).EDFFile == ind)
      LocalRemovePlugin(i);
    end
  end
end
LocalRepaint;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LocalSelectEDFFile
% display dialog to select a loaded EDF file
function [FileNum, Cancelled] = LocalSelectEDFFile(Title, EDFData)
dlgh = dialog(...
    'Name', Title{1}, ...
    'CloseRequestFcn', 'set(gcf,''UserData'',''Cancel'');uiresume;');
dlgpos = get(dlgh, 'Position');
set(dlgh, 'Position', [dlgpos(1),dlgpos(2),300,100]);

EDFNames = {};
% get EDF FIle names
for i = 1:length(EDFData)
  EDFNames = {EDFNames{:}, EDFData(i).Head.FileName};
end

poph = uicontrol(dlgh, ...
    'Style', 'Popup', ...
    'Units', 'Normalized', ...
    'String',EDFNames, ...
    'Value', 1, ...
    'Position', [0.05, 0.5, 0.9, 0.35]);
% buttons
uicontrol(dlgh, ...
    'Style', 'PushButton', ...
    'Units', 'Normalized', ...
    'String', Title{2}, ...
    'Position', [0.1, 0.05, 0.3, 0.3], ...
    'Callback', 'set(gco,''UserData'',''OK'');uiresume;');
uicontrol(dlgh, ...
    'Style', 'PushButton', ...
    'Units', 'Normalized', ...
    'String', 'Cancel', ...
    'Position', [0.6, 0.05, 0.3, 0.3], ...
    'Callback', 'set(gco,''UserData'',''Cancel'');uiresume;');

drawnow;
uiwait(dlgh);
if strcmp(get(gco,'UserData'),'OK')
  Cancelled = 0;
  FileNum = get(poph, 'Value');
else
  Cancelled = 1;
  FileNum = 0;  
end
delete(dlgh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LocalCloseViewEDF
% close viewer
function LocalCloseViewEDF()
Data = get(findobj('Tag', 'ViewEDFFigure'), 'UserData');
% close all files
for i = 1:length(Data.EDF);
  %fclose(Data.EDF(i).Head.FILE.FID);
  sclose(Data.EDF(i).Head);
end
delete(findobj('Tag', 'ViewEDFFigure'));
