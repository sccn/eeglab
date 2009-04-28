function [cmap] = colormapsmake(bins,fighandle)

% colorMapsMake - Interactive color map maker
%
%     Useage:  colormapsmake
%               colormapsmake(map)
%               colormapsmake(map,fighandle)
%               colormapsmake(bins)
%               colormapsmake(bins,fighandle)
%               colormapsmake(fighandle,'fig')
%
%       bins - size of colormap (default - 64)
%       map - rgb colormap (matrix)
%       fighandle - figure(s) to apply colormap (optional)
%       fighandle,'fig' - reads colormap from figure (fighandle)
%
%    Nodes are displayed as circles. To move a node use mouse
%       button #1. To add a node use mouse button #2. To delete use 
%       button #3.
%
%    Mouse buttons     button 1 - move    (3 button mouse)
%                      button 2 - add
%                      button 3 - erase
%
%    There are four graphs. The first three corresond to the relative
%       rgb values. And the fourth graph is a luminance factor
%       which the rgb values are multiplied by.
%
%     colormap = [red,green,blue] * luminance
%
%    Reset button - resets the colormap to the original values.
%    Export button - assigns the colormap to a desktop variable.
%    Close button - closes the figure window.
%
%    Runs under Matlab 5.0 through 6.0
%

% Note: the readcolormap option is not perfect. It works well if the
% colormap contains straight lines like most matlab colormaps. You can play
% with the tolerance. Also, the more nodes in an axes the slower the mouse
% handler routines will run.

% Written by Colin Humphries, Salk Institute, March 1997
%   -add lspace function (now runs under matlab 5.2) April, 1998
%   email: chumphri@uci.edu
%
% Darren.Weber_at_radiology.ucsf.edu
% Sept. 2001 - renamed from makecmap to colorMapsMake.m
%            - added help button and minor cosmetics
%            - runs on Matlab 6.0 (R12)

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:50 $


nargs = nargin;


if nargs == 0
  bins = 64;
  nargs = 1;
end


if ~isstr(bins)
  
  if size(bins,2) == 3
    map = bins;
    bins = size(map,1);
    if nargs == 1
      fighandle = [];
    end
  else
    map = [];
    if nargs == 2
      if isstr(fighandle)
        fighandle = bins;
        map = get(fighandle,'Colormap');
        bins = size(map,1);
      end
    else
      fighandle = [];
    end
  end
  
  if isempty(map)
    red = [1 0;bins 1];       % default colormap
    green = [1 0;bins 1];
    blue = [1 0;bins 1];
  else
                          % Note: this routine cycles through all
                              % the data points and assigns nodes to 
                              % the points where straight lines end.
                              % This does not work well in every 
                              % occasion. If the colormap is nonlinear
                              % in many places, a lot of nodes will be
                              % assigned.
        
    for j = 1:3           % cycle through each column
      for i = 1:bins      % cycle through each row
        if i == 1         % assign the first point as the first node
          nodes(1,:) = [1 map(1,j)];
          lastnode = nodes;   % last node assigned
        elseif i == 2
          oldtempnode = [i map(2,j)]; % assigns point 2 as the test node
        else
          tempnode = [i map(i,j)];    % look at current point
          if (abs((oldtempnode(2)-lastnode(2))/(oldtempnode(1)...
                  -lastnode(1))) <= abs((tempnode(2)-lastnode(2))...
                  /(tempnode(1)-lastnode(1)))*1.01) & ...
                  (abs((oldtempnode(2)-lastnode(2))/...
                  (oldtempnode(1)-lastnode(1))) >= ...
                  abs((tempnode(2)-lastnode(2))/...
                  (tempnode(1)-lastnode(1)))*.99) % Check whether the line
                                                  % from the last node to
                                                  % tempnode and
                                                  % oldtempnode have the
                                                  % same slope.
                                                  % Tolerance is +/- 1%
                                                    
            oldtempnode = tempnode;    % if so then move oldtempnode
          else
            nodes = [nodes;oldtempnode]; % if not then assign a new
            lastnode = oldtempnode;      % node.
            oldtempnode = tempnode;
          end
        end
      end
      nodes = [nodes;bins map(bins,j)];  % the last node is the last point
      if j == 1
        red = nodes;  % get red nodes
        nodes = []; 
      elseif j == 2
        green = nodes;% get green nodes
        nodes = [];
      elseif j == 3
        blue = nodes; % get blue nodes
        nodes = [];
      end
    end
  end
  
  mcmfighandle = figure('color','k','NumberTitle','off',...
                        'menubar','none','Tag','ColorMap Editor');
  
  if isempty(fighandle)
    fighandle = mcmfighandle;
  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up Colormap from nodes
%   Converts each set of nodes into a discrete vector


  for i = 1:size(red,1)-1
    W = lspace(red(i,2),red(i+1,2),red(i+1,1)-red(i,1)+1);
    mapred(red(i,1):red(i+1,1)-1) = W(1:red(i+1,1)-red(i,1))';
  end
  mapred(red(i+1,1)) = red(i+1,2);
  for i = 1:size(green,1)-1
    W = lspace(green(i,2),green(i+1,2),green(i+1,1)-green(i,1)+1);
    mapgreen(green(i,1):green(i+1,1)-1) = W(1:green(i+1,1)-green(i,1))';
  end
  mapgreen(green(i+1,1)) = green(i+1,2);
  for i = 1:size(blue,1)-1
    W = lspace(blue(i,2),blue(i+1,2),blue(i+1,1)-blue(i,1)+1);
    mapblue(blue(i,1):blue(i+1,1)-1) = W(1:blue(i+1,1)-blue(i,1))';
  end
  mapblue(blue(i+1,1)) = blue(i+1,2);
  map = [mapred(:),mapgreen(:),mapblue(:)]; % get new colormap
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up Graphs  


  nodenum = size(red,1);
  subplot(3,2,1)                 % Red Graph
  hold on
  for i = 1:nodenum-1
    plot(red(i,1),red(i,2),'ow')
    line([red(i,1),red(i+1)],[red(i,2),red(i+1,2)],'color','r')
  end
  plot(red(nodenum,1),red(nodenum,2),'ow') 
  axis([1 red(nodenum,1) 0 1])
  redax = gca;
  set(redax,'UserData',red,'tag','r','box','on','color','k',...
      'xcolor','w','ycolor','w','xgrid','on','ygrid','on',...
      'Ytick',(0:.2:1));
  ylabel('Red','FontSize',14)


  subplot(3,2,3)                % Green Graph
  nodenum = size(green,1);
  hold on
  for i = 1:nodenum-1
    plot(green(i,1),green(i,2),'ow')
    line([green(i,1),green(i+1)],[green(i,2),green(i+1,2)],'color','g')
  end
  plot(green(nodenum,1),green(nodenum,2),'ow') 
  axis([1 green(nodenum,1) 0 1])
  greenax = gca;
  set(greenax,'UserData',green,'tag','g','box','on','color','k',...
      'xcolor','w','ycolor','w','xgrid','on','ygrid','on',...
      'YTick',(0:.2:1));
  ylabel('Green','FontSize',14)
  
  subplot(3,2,5)                % Blue Graph
  nodenum = size(blue,1);
  hold on
  for i = 1:nodenum-1
    plot(blue(i,1),blue(i,2),'ow')
    line([blue(i,1),blue(i+1)],[blue(i,2),blue(i+1,2)],'color','b')
  end
  plot(blue(nodenum,1),blue(nodenum,2),'ow') 
  axis([1 blue(nodenum,1) 0 1])
  blueax = gca;
  set(blueax,'UserData',blue,'tag','b','box','on','color','k',...
      'xcolor','w','ycolor','w','xgrid','on','ygrid','on',...
      'YTick',(0:.2:1));
  ylabel('Blue','FontSIze',14)


  colormap(map)                % Apply Colormap
  set(fighandle,'Colormap',map)
  
  transax = axes('Position',[.578 .4056 .327 .2238]); % Transfer Graph
  hold on
  trans = [1 1;bins 1];
  nodenum = 2;
  cla
  for i = 1:nodenum-1
    plot(trans(i,1),trans(i,2),'ow')
    line([trans(i,1),trans(i+1)],[trans(i,2),trans(i+1,2)],'color','y')
  end
  plot(trans(nodenum,1),trans(nodenum,2),'ow') 
  axis([1 trans(nodenum,1) 0 1])
  trans = ones(1,bins);
  title('luminance','color','w','FontSize',14)
  set(transax,'UserData',[1 1;bins 1],'tag','y','box','on','color','k',...
      'xcolor','w','ycolor','w','xgrid','on','ygrid','on',...
      'YTick',(0:.2:1))
  set(mcmfighandle,'UserData',[map,trans(:)],'tag',int2str(fighandle))


  cbax = axes('Position',[.578 .75 .327 .1072]);  % Colorbar Graph
  colorbar(cbax)
  set(cbax,'tag','cb','xcolor','w','ycolor','w',...
      'UserData',[size(red,1),0;...
                  size(green,1),0;...
                  size(blue,1),0;...
                  red;...
                  green;...
                  blue])
  
% Title axis


  hdaxes = axes('Position',[.57 .93 .327 .05],'Visible','off',...
      'tag','hd','UserData',map);
  text(.5,0,'Colormap Editor','FontSize',16,'Color','w',...
      'HorizontalAlignment','center')


% Uicontrols
  
  uicontrol('style','pushbutton','string','Reset',...
      'Units','Normalized',...
      'Position',[.6  .2 .1 .06],...
      'Callback','ColorMapsMake(''reset'');')
  
  uicontrol('style','pushbutton','string','Export',...
      'Units','Normalized',...
      'Position',[.75 .2 .1 .06],...
      'Callback','ColorMapsMake(''export'');')
  
  uicontrol('style','pushbutton','string','Close',...
      'Units','Normalized',...
      'Position',[.6  .1 .1 .06],...
      'Callback','ColorMapsMake(''close'');')
  
  uicontrol('Style','PushButton','string','Help',...
      'Units','Normalized',...
      'Position',[.75 .1 .1 .06],...
      'Callback','ColorMapsMake(''help'');')
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up callbacks


  set(mcmfighandle,'WindowButtonDownFcn','ColorMapsMake(''down'');')
  set(mcmfighandle,'WindowButtonUpFcn','ColorMapsMake(''up'');')
  set(mcmfighandle,'WindowButtonMotionFcn','ColorMapsMake(''motion'');')
  set(mcmfighandle,'CloseRequestFcn','ColorMapsMake(''close'');')


% End of main function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


else
  if strcmp(bins,'down')           % Mouse Button Down Handler
    figh = gcbf;
    axhandle = gca;
    lcolor = get(axhandle,'tag');  % find out what axis we are in
    if strcmp(lcolor,'cb') | strcmp(lcolor,'hd')
      return           % if colorbar or title axis then return
    end
    xlimits = get(axhandle,'Xlim');
    ylimits = get(axhandle,'Ylim');
    nodes = get(axhandle,'UserData');
    nodenum = size(nodes,1);      % get nodes
    P = get(axhandle,'CurrentPoint'); % get mouse position
    X = P(1,1);
    Y = P(1,2);
    button = get(gcbf,'SelectionType'); % get button
    if X <= xlimits(1)-xlimits(2)*.05 | X >= xlimits(2)*1.05 | ...
          Y <= ylimits(1)-ylimits(2)*.07 | Y >= ylimits(2)*1.07
      return    % if outside of axis limits then return
    end


    if strcmp(button,'normal')    % left mouse button
      global M_BUTTON_DOWN M_PNTS M_LINES M_TEXT
      ii = find(X > nodes(:,1));  % find node above point
      minnode = max(ii);
      ii = find(X < nodes(:,1));  % find node below point
      maxnode = min(ii);
      if (X - nodes(minnode,1)) > (nodes(maxnode,1)-X) % if closer to maxnode
         if maxnode == nodenum;     % if maxnode is the endnode
           nodes(maxnode,2) = max(min(Y,1),0);  % only move Y-position
         else
           nodes(maxnode,1) = round(X);         % move node
           nodes(maxnode,2) = max(min(Y,1),0);
         end 
         M_BUTTON_DOWN = maxnode;               % save what node was selected
      elseif (X - nodes(minnode,1)) < (nodes(maxnode,1)-X)% if closer to minnode
         if minnode == 1     %if firstnode
           nodes(minnode,2) = max(min(Y,1),0); % only move Y-position
         else
           nodes(minnode,1) = round(X);        % move node
           nodes(minnode,2) = max(min(Y,1),0);
         end
         M_BUTTON_DOWN = minnode;              % save what node was selected
      end      
      cla
      for i = 1:nodenum-1                 % Redraw axis
        M_PNTS(i) = plot(nodes(i,1),nodes(i,2),'ow');
        M_LINES(i) = line([nodes(i,1),nodes(i+1)],...
            [nodes(i,2),nodes(i+1,2)],'color',lcolor);
      end
      M_PNTS(i+1) = plot(nodes(nodenum,1),nodes(nodenum,2),'ow');
      
      M_TEXT = text(0,1.05,...
          [' ',int2str(nodes(M_BUTTON_DOWN,1)),', ',...
            num2str(nodes(M_BUTTON_DOWN,2),2)],...
          'Color','w','clipping','off',...
          'VerticalAlignment','bottom',...
          'FontSize',10);
      
      axis([1 nodes(nodenum,1) 0 1])    


    elseif strcmp(button,'extend')     % Middle Mouse Button
      global M_BUTTON_DOWN M_PNTS M_LINES M_TEXT
      
      ii = find(X > nodes(:,1));       % find node below
      minnode = max(ii);
      ii = find(X < nodes(:,1));       % find node above
      maxnode = min(ii);
      nodes = [nodes(1:minnode,:);[round(X),max(min(Y,1),0)];...
          nodes(maxnode:nodenum,:)];   % add node shifting those above
      nodenum = size(nodes,1);
      M_BUTTON_DOWN = minnode+1;
      cla
      for i = 1:nodenum-1              % Redraw axis
        M_PNTS(i) = plot(nodes(i,1),nodes(i,2),'ow');
        M_LINES(i) = line([nodes(i,1),nodes(i+1)],...
            [nodes(i,2),nodes(i+1,2)],'color',lcolor);
      end
      M_PNTS(i+1) = plot(nodes(nodenum,1),nodes(nodenum,2),'ow'); 
      M_TEXT = text(0,1.05,...
          [' ',int2str(nodes(M_BUTTON_DOWN,1)),', ',...
            num2str(nodes(M_BUTTON_DOWN,2),2)],...
          'Color','w','clipping','off',...
          'VerticalAlignment','bottom',...
          'FontSize',10);
      
      axis([1 nodes(nodenum,1) 0 1])
      
    elseif strcmp(button,'alt')          % Right Mouse Button
      ii = find(X > nodes(:,1));         % find node below
      minnode = max(ii);
      ii = find(X < nodes(:,1));         % find node above
      maxnode = min(ii);
      if (X - nodes(minnode,1)) > (nodes(maxnode,1)-X)    % reassign nodes
        if maxnode < nodenum;
          nodes = [nodes(1:maxnode-1,:);nodes(maxnode+1:nodenum,:)];
        end
      elseif (X - nodes(minnode,1)) < (nodes(maxnode,1)-X)
        if minnode > 1
          nodes = [nodes(1:minnode-1,:);nodes(minnode+1:nodenum,:)];
        end
      end
      nodenum = size(nodes,1);
      cla
      for i = 1:nodenum-1       %redraw axis
         plot(nodes(i,1),nodes(i,2),'ow')
         line([nodes(i,1),nodes(i+1)],[nodes(i,2),nodes(i+1,2)],'color',lcolor)
      end
      plot(nodes(nodenum,1),nodes(nodenum,2),'ow') 
      axis([1 nodes(nodenum,1) 0 1])
    end
    set(axhandle,'UserData',nodes)  % update nodes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


elseif strcmp(bins,'up')            % Mouse Button Up Handler
    figh = gcbf;
    global M_TEXT
    delete(M_TEXT)
    clear global M_BUTTON_DOWN M_LINES M_PNTS M_TEXT % clear variables
    UserData = get(figh,'UserData');
    fighandle = str2num(get(figh,'tag')); % get figure to apply colormap
    axhandle = gca;
    nodes = get(axhandle,'UserData');    % get nodes
    tagval = get(axhandle,'tag');
    if strcmp(tagval,'cb') | strcmp(tagval,'hd')
      return   % if colorbar or title axis then return
    end
    for i = 1:size(nodes,1)-1
      W = lspace(nodes(i,2),nodes(i+1,2),nodes(i+1,1)-nodes(i,1)+1);
      map(nodes(i,1):nodes(i+1,1)-1) = W(1:nodes(i+1,1)-nodes(i,1))';
    end
    map(nodes(i+1,1)) = nodes(i+1,2);
    if strcmp(tagval,'r')
      UserData(:,1) = map';
    elseif strcmp(tagval,'g')
      UserData(:,2) = map';
    elseif strcmp(tagval,'b')
      UserData(:,3) = map';
    elseif strcmp(tagval,'y')
      UserData(:,4) = map';
    end
    set(figh,'UserData',UserData,'Colormap',...      % assign new colormap
        UserData(:,1:3).*(UserData(:,4)*ones(1,3)))
    if fighandle ~= figh
      set(fighandle,'Colormap',UserData(:,1:3).*(UserData(:,4)*ones(1,3)))
    end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  elseif strcmp(bins,'motion')      % Mouse Motion Handler
    global M_BUTTON_DOWN M_LINES M_PNTS M_TEXT
    if isempty(M_BUTTON_DOWN)
      return     % If mouse button is not down then return
    end   
    axhandle = gca;
    nodes = get(gca,'UserData');
    tagval = get(gca,'tag');
    P = get(axhandle,'CurrentPoint');
    X = P(1,1);
    Y = P(1,2);
    nodenum = size(nodes,1);
    if M_BUTTON_DOWN == 1 | M_BUTTON_DOWN == nodenum; % move selected node to 
      nodes(M_BUTTON_DOWN,2) = max(min(Y,1),0);       % new mouse position
    else
      nodes(M_BUTTON_DOWN,1) = max(nodes(M_BUTTON_DOWN-1,1), ...
          min(nodes(M_BUTTON_DOWN+1,1),round(X)));
      nodes(M_BUTTON_DOWN,2) = max(min(Y,1),0);
    end
    
    delete(M_PNTS(M_BUTTON_DOWN))
    M_PNTS(M_BUTTON_DOWN) = plot(nodes(M_BUTTON_DOWN,1),...
        nodes(M_BUTTON_DOWN,2),'ow');
    if M_BUTTON_DOWN > 1
      delete(M_LINES(M_BUTTON_DOWN-1))
      M_LINES(M_BUTTON_DOWN-1) = line([nodes(M_BUTTON_DOWN-1,1),...
            nodes(M_BUTTON_DOWN,1)],...
          [nodes(M_BUTTON_DOWN-1,2),nodes(M_BUTTON_DOWN,2)],'color',tagval);
    end
    if M_BUTTON_DOWN < nodenum
      delete(M_LINES(M_BUTTON_DOWN))
      M_LINES(M_BUTTON_DOWN) = line([nodes(M_BUTTON_DOWN,1),...
            nodes(M_BUTTON_DOWN+1,1)],...
          [nodes(M_BUTTON_DOWN,2),nodes(M_BUTTON_DOWN+1,2)],'color',tagval);
    end  
    
    delete(M_TEXT)
    M_TEXT = text(0,1.05,...
          [' ',int2str(nodes(M_BUTTON_DOWN,1)),', ',...
            num2str(nodes(M_BUTTON_DOWN,2),2)],...
          'Color','w','clipping','off',...
          'VerticalAlignment','bottom',...
          'FontSize',10);
    
    axis([1 nodes(nodenum,1) 0 1])
    set(axhandle,'UserData',nodes)
    
  elseif strcmp(bins,'export')   % export colormap to desktop
    fig = gcbf;
    pos = get(fig,'Position');
    figx = 400;                  % set up dialog figure
    figy = 200;
    figh = figure('Units','pixels',...
        'Position',...
        [pos(1)+pos(3)/2-figx/2 pos(2)+pos(4)/2-figy/2 figx figy],...
        'Resize','off','CloseRequestFcn','','menubar','none',...
        'numbertitle','off');
    uicolor = get(figh,'Color');
    
    % text uicontrol
    uicontrol('Style','Text','Units','Pixels',...
        'String','Output Colormap to Desktop Variable:',...
        'Position',[20 figy-40 300 22],'HorizontalAlignment','left',...
        'FontSize',14,'BackGroundColor',uicolor)
    
    % edit uicontrol
    ui1 = uicontrol('Style','Edit','Units','Pixels',...
        'String','cmap','FontSize',12,...
        'Position',[120 figy-100 150 30]);
    
    TIMESTRING = ['[OBJ1,FIGH1] = gcbo;',...
                  'OBJHAN = get(OBJ1,''UserData'');',...
                  'LAB1 = get(OBJHAN(1),''string'');',...
                  'DATA1 = get(OBJHAN(2),''Colormap'');',...
                  'eval([LAB1,'' = DATA1;'']);',...
                  'delete(FIGH1);',...
                  'clear OBJ1 FIGH1 OBJHAN DATA1 LAB1'];
    
    % OK Button
    uicontrol('Style','PushButton','Units','Pixels',...
        'String','OK','FontSize',14,...
        'Position',[figx/4-20 10 65 30],...
        'UserData',[ui1 fig],'Callback',TIMESTRING)
    
    TIMESTRING = ['[OBJ1,FIGH1] = gcbo;',...
                'delete(FIGH1);',...
                'clear OBJ1 FIGH1;'];
    
    % Cancel Button
    uicontrol('Style','PushButton','Units','Pixels',...
        'String','Cancel','FontSize',14,...
        'Position',[3*figx/4-20 10 65 30],...
        'Callback',TIMESTRING)
    
  elseif strcmp(bins,'help')   % Help function
      
      helpstr = sprintf(['Left Button - move RGB node\n',...
                         'Middle Button - add node to RGB\n',...
                         'Right Button - remove node (>2)\n\n',...
                         'Export button - assign colormap to variable']);
      
      help = helpdlg(helpstr,'ColorMap HELP');
      
      
  elseif strcmp(bins,'close')   % Close Request function
    figh = gcbf;
    clear global M_BUTTON_DOWN M_LINES M_PNTS M_TEXT % get rid of global vars
    delete(figh)
  
  elseif strcmp(bins,'reset')   % reset Colormap to original values
    fighandle = gcbf;
    obj = findobj('tag','hd','parent',fighandle);
    map = get(obj,'UserData');  % get original colormap stored in hd axes
    set(fighandle,'Colormap',map,...
        'UserData',[map , ones(size(map,1),1)])
    
    obj = findobj('tag','cb','parent',fighandle);
    nodes = get(obj,'UserData');  % get node values stored in cb axes
    red = nodes(4:nodes(1,1)+3,:);
    green = nodes(nodes(1,1)+4:nodes(1,1)+nodes(2,1)+3,:);
    blue = nodes(nodes(1,1)+nodes(2,1)+4:nodes(1,1)+nodes(2,1)+...
        nodes(3,1)+3,:);
    
    obj = findobj('tag','r','parent',fighandle);
    axes(obj)                         % reset colormap graphs
    cla
    nodenum = nodes(1,1);
    for i = 1:nodenum-1
      plot(red(i,1),red(i,2),'ow')
      line([red(i,1),red(i+1)],[red(i,2),red(i+1,2)],'color','r')
    end
    plot(red(nodenum,1),red(nodenum,2),'ow')
    set(obj,'UserData',red)
    
    obj = findobj('tag','g','parent',fighandle);
    axes(obj)
    cla
    nodenum = nodes(2,1);
    for i = 1:nodenum-1
      plot(green(i,1),green(i,2),'ow')
      line([green(i,1),green(i+1)],[green(i,2),green(i+1,2)],'color','g')
    end
    plot(green(nodenum,1),green(nodenum,2),'ow')
    set(obj,'UserData',green)
    
    obj = findobj('tag','b','parent',fighandle);
    axes(obj)
    cla
    nodenum = nodes(3,1);
    for i = 1:nodenum-1
      plot(blue(i,1),blue(i,2),'ow')
      line([blue(i,1),blue(i+1)],[blue(i,2),blue(i+1,2)],'color','b')
    end
    plot(blue(nodenum,1),blue(nodenum,2),'ow')
    set(obj,'UserData',blue)
    
    obj = findobj('tag','y','parent',fighandle);
    axes(obj)
    cla
    trans = [1 1;size(map,1) 1];
    nodenum = 2;
    for i = 1:nodenum-1
      plot(trans(i,1),trans(i,2),'ow')
      line([trans(i,1),trans(i+1)],[trans(i,2),trans(i+1,2)],'color','y')
    end
    plot(trans(nodenum,1),trans(nodenum,2),'ow') 
    set(obj,'UserData',trans)
    figh = str2num(get(fighandle,'tag'));
    if figh~=fighandle
      set(figh,'Colormap',map)
    end
  end
end


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = lspace(x,y,N)
% LSPACE
%   Unfortunately, linspace.m changed in Matlab 5.2, and it doesn't work
%   the same. The new function cannot handle the command linspace(0,0,24). 
%   It just returns an empty matrix.


out = [x+(0:N-2)*(y-x)/(N-1) y];

return
