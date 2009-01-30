function [Xpoint,Ypoint,XYhist] = crosshair(action);
%  CROSSHAIR:  A gui interface for reading (x,y) values from a plot.
%  
%  [Xpoint,Ypoint,XYhist] = crosshair(action);
%  
%  A set of mouse driven crosshairs is placed on the current axes,
%  and displays the current (x,y) values of the line plot.  There is an
%  option to provide data specific values or interpolated values.  The
%  resolution of the data specific values depends on both the data
%  resolution and the GUI interface (mainly mouse movement resolution).
%  The interpolated values appear to provide a more continuous function,
%  however they too depend on the GUI interface resolution.  There are 
%  no options for extrapolation.
%  
%  For multiple traces, plots with the same length(xdata) are
%  tracked. Each mouse click returns Xpoint,Ypoint values and selecting 
%  'done' will remove the GUI and restore the mouse buttons to previous 
%  values.  Selecting 'exit' will remove the GUI and close the figure.
%  
%  In this version (Dec 2002), there are new options to define functions
%  of X/Y and monitor their results.  There is also a new STORE button 
%  that creates and updates an XYhist struct in the base workspace, which
%  contains value labels and values.  This version has better controls 
%  of X/Y movements, including better interpolation movement options.  
%  This version attempts to respond correctly to keyboard entries also.
%  
%  Some further help is given in the tool tips of the GUI.
%
%  Note: crosshair always updates the Xpoint,Ypoint in the base 
%  workspace. Here is an example of how to get return values within
%  a script/function after pressing the exit button of crosshair:
%       function [X,Y] = crosshair_returnXY
%		x = [1:10]; y(1,:) = sin(x); y(2,:) = cos(x);
%		figure; plot(x,y); crosshair;
%		uiwait
%		X = evalin('base','Xpoint');
%		Y = evalin('base','Ypoint');
%		return
%  Copy this text to a function .m file and then call it from the
%  base workspace with [X,Y] = crosshair_returnXY
%  
%  Useage:  x = [1:10]; y(1,:) = sin(x); y(2,:) = cos(x); x2 = x.^2;
%           figure; plot(x2,y); crosshair
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2002  Darren L. Weber
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

%  History: 03/96, Richard G. Cobb <cobbr@plk.af.mil>
%           08/01, Darren.Weber@flinders.edu.au
%                  replaced obsolete 'table1' with 'interp1'; fixed bug
%                  with number of 'traces'; rationalized calculations into
%                  a common subfunction for x,y point calc in 'down','up', 
%                  & 'move' button functions; added option to turn on/off
%                  interpolation and the exit button; simplified updates 
%                  to graphics using global GUI handle structure.
%           11/01, Darren.Weber@flinders.edu.au
%                  added tooltips for several GUI handles
%                  added multiple interpolation methods
%                  added GUI for data matrix indices (given no interpolation)
%                  added option to select trace nearest to mouse click point
%                  reversed order of lines in data matrix to be consistent
%                    with the value returned from the nearest trace subfunction
%                  create crosshair lines after finding all plot lines to
%                    avoid confusing them with the plot lines
%           01/02, Darren.Weber@flinders.edu.au
%                  should now work across multiple plot figures, given
%                    that all gui handles and data are now stored in the
%                    plot figure 'userdata' handle.
%                  added functionality to move smoothly from interpolation
%                    back to the delimited data via the "next/previous" 
%                    buttons.
%           06/02, Darren.Weber@flinders.edu.au
%                  learned how to get values back to a script/function with
%                    evalin command and updated help above.
%           12/02, Darren.Weber@flinders.edu.au
%                  added Store uicontrol and associated XYhist variable 
%                   updates in base workspace to provide storage of 
%                   consecutive XY values (thanks to C.A.Swenne@lumc.nl &
%                   H.van_de_Vooren@lumc.nl for their suggestion/assistance).
%                  added keyboard controls to left/right arrow and return.
%                  added prev/next X interpolation interval.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if ~exist('action','var') 
    action = 'init';
elseif isempty(action)
    action = 'init';
end

XHR_HANDLES = get(gcbf,'userdata');


% Check for specific keys and assign reasonable actions
if strcmp(action, 'keypress'),
    
    CC = get(XHR_HANDLES.gui,'CurrentCharacter');
    cc = double(CC);
    if cc,
        switch cc,
        case 27, action = 'done';  % ESC
        case 28, action = 'prevx'; % left
        case 29, action = 'nextx'; % right
        case 30, action = 'ygt';   % up
        case 31, action = 'ylt';   % down
        case 13, action = 'store'; % return/enter
        otherwise, action = 'up';  % all other keys
        end
    end
end

action = lower(action);

switch action,

case 'init',
    
    % Paint GUI
    XHR_HANDLES = INIT;
    
    % Update and return values
    XHR_HANDLES = updateDATA(XHR_HANDLES);
    Xpoint = get(XHR_HANDLES.xvalue,'Value');
    Ypoint = get(XHR_HANDLES.yvalue,'Value');
    %updateXYhistory(Xpoint,Ypoint);
    
% Mouse Click Down
case 'down',
    
    set(XHR_HANDLES.gui,'WindowButtonMotionFcn','crosshair(''move'');');
    set(XHR_HANDLES.gui,'WindowButtonUpFcn','[Xpoint,Ypoint] = crosshair(''up'');');
    
    XHR_HANDLES = updateDATA(XHR_HANDLES);
    Xpoint = get(XHR_HANDLES.xvalue,'Value');
    Ypoint = get(XHR_HANDLES.yvalue,'Value');
    %updateXYhistory(Xpoint,Ypoint);
    
% Mouse Drag Motion
case 'move',
    
    XHR_HANDLES = updateDATA(XHR_HANDLES);
    Xpoint = get(XHR_HANDLES.xvalue,'Value');
    Ypoint = get(XHR_HANDLES.yvalue,'Value');
    %updateXYhistory(Xpoint,Ypoint);
    
% Mouse Click Up
case 'up',
    
    set(XHR_HANDLES.gui,'WindowButtonMotionFcn',' ');
    set(XHR_HANDLES.gui,'WindowButtonUpFcn',' ');
    
    XHR_HANDLES = updateDATA(XHR_HANDLES);
    Xpoint = get(XHR_HANDLES.xvalue,'Value');
    Ypoint = get(XHR_HANDLES.yvalue,'Value');
    %updateXYhistory(Xpoint,Ypoint);
    
% Next or Previous X point
case {'nextx','prevx','changex','nexty','prevy','changey','ylt','ygt'}, % Change X/Y
    
    XHR_HANDLES = moveXY(XHR_HANDLES,action);
    Xpoint = get(XHR_HANDLES.xvalue,'Value');
    Ypoint = get(XHR_HANDLES.yvalue,'Value');
    %updateXYhistory(Xpoint,Ypoint);
    
% Store XY values into a history array
case 'store',
    
    Xpoint = get(XHR_HANDLES.xvalue,'Value');
    Ypoint = get(XHR_HANDLES.yvalue,'Value');
    updateXYhistory(XHR_HANDLES);
    
% Exit crosshairs GUI
case {'done','exit'},
    
    XHR_HANDLES = get(gcbf,'userdata');
    
    Xpoint = get(XHR_HANDLES.xvalue,'Value');
    Ypoint = get(XHR_HANDLES.yvalue,'Value');
    %updateXYhistory(Xpoint,Ypoint);
    
    handles = fieldnames(XHR_HANDLES);
    for i=1:length(handles),
        switch handles{i},
        case {'axis','datalines','gui'},
        otherwise,
            h = getfield(XHR_HANDLES,handles{i});
            if ishandle(h), delete(h); end
        end
    end
    
    if strcmp(action,'exit');
        if ishandle(XHR_HANDLES.gui),
            close(XHR_HANDLES.gui);
        end
    else
        
        set(XHR_HANDLES.gui,'WindowButtonUpFcn','');
        set(XHR_HANDLES.gui,'WindowButtonMotionFcn','');
        set(XHR_HANDLES.gui,'WindowButtonDownFcn',XHR_HANDLES.button);
        
        refresh(XHR_HANDLES.gui);
    end
    
    clear XHR_HANDLES;
    return;
    
end

set(gcf,'userdata',XHR_HANDLES);
if ishandle(XHR_HANDLES.gui),
    figure(XHR_HANDLES.gui);
end

return;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateXYhistory(H),
    
    Ch = get(H.yindex,'Value');
    X  = get(H.xvalue,'Value');
    Y  = get(H.yvalue,'Value');
    fX = get(H.fxvalue,'Value');
    fY = get(H.fyvalue,'Value');
    
    fXeq = get(H.fxeq,'String');
    fYeq = get(H.fyeq,'String');
    
    XY.labels = {'Channel','X','Y',['f(x) = ',fXeq],['f(y) = ',fYeq]};
    XY.data   = [Ch,X,Y,fX,fY];
    
    if evalin('base','exist(''XYhist'',''var'');'),
        XYhist = evalin('base','XYhist');
        % get the current history set
        set = getfield(XYhist,['set',num2str(XYhist.set)]);
        if isequal(set.labels,XY.labels),
            set.data(end+1,:) = XY.data;
            XYhist = setfield(XYhist,['set',num2str(XYhist.set)],set);
            assignin('base','XYhist',XYhist);
        else
            fprintf('\nWarning: creating new set of XYhist in base workspace.\n\n');
            XYhist.set = XYhist.set + 1;
            XYhist = setfield(XYhist,['set',num2str(XYhist.set)],XY);
            assignin('base','XYhist',XYhist);
        end
    else
        XYhist.set = 1;
        XYhist.set1 = XY;
        assignin('base','XYhist',XYhist);
    end
    
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H] = moveXY(H,move)
    
    interp  = get(H.interp,'Value');
    xinterp = get(H.xinterp,'Value');
    
    if (xinterp > 0) & (interp > 1),
        
        % Use incremental interpolation of x values
        switch move,
        case 'nexty', H.data.yindex = H.data.yindex + 1;
        case 'prevy', H.data.yindex = H.data.yindex - 1;
        case 'ygt',
            ydata = interpYall(H);
            [ysort,yi] = sort(ydata);
            currYI = find(ysort > H.data.ypoint);
            if min(currYI),
                H.data.yindex = yi(min(currYI));
            end
        case 'ylt',
            ydata = interpYall(H);
            [ysort,yi] = sort(ydata);
            currYI = find(ysort < H.data.ypoint);
            if max(currYI),
                H.data.yindex = yi(max(currYI));
            end
        case 'nextx',
            H.data.xpoint = H.data.xpoint + xinterp;
        case 'prevx',
            H.data.xpoint = H.data.xpoint - xinterp;
        end
        H = checkdatarange(H);
        H = interpY(H);
        updateGUI(H);
        return
    end
    
    
    % No interpolation of x values...
    
    if (interp > 1)
        xdata = H.data.xdata(:,H.data.yindex);
        [H.data.xindex] = NearestXYArrayPoint( xdata, H.data.xpoint, move );
    end
    
    switch move,
    case 'nextx',
        % Increase current xindex by one
        if(interp == 1), H.data.xindex = H.data.xindex + 1; end
    case 'prevx',
        % Decrease current xindex by one
        if(interp == 1), H.data.xindex = H.data.xindex - 1; end
    case 'nexty', H.data.yindex = H.data.yindex + 1;
    case 'prevy', H.data.yindex = H.data.yindex - 1;
    case 'ygt',
        ydata = H.data.ydata(H.data.xindex,:);
        [ysort,yi] = sort(ydata);
        currYI = find(ysort == H.data.ypoint);
        if currYI < length(yi),
            H.data.yindex = yi(currYI+1);
        end
    case 'ylt',
        ydata = H.data.ydata(H.data.xindex,:);
        [ysort,yi] = sort(ydata);
        currYI = find(ysort == H.data.ypoint);
        if currYI > 1,
            H.data.yindex = yi(currYI-1);
        end
    otherwise
    end
    
    H = checkdatarange(H);
    
    % Get x/y value at new x/y index
    H.data.xpoint = H.data.xdata(H.data.xindex,H.data.yindex);
    H.data.ypoint = H.data.ydata(H.data.xindex,H.data.yindex);
    
    set(H.interp,'Value',1);
    H = updateGUI(H);
    
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ H ] = checkdatarange(H),

    % Ensure that x/y index is within data range
    s = size(H.data.xdata,1);
    if( H.data.xindex < 1 ),
        H.data.xindex = 1;
    elseif( H.data.xindex >= s ),
        H.data.xindex = s;
    end
    s = size(H.data.ydata,2);
    if( H.data.yindex < 1 ),
        H.data.yindex = 1;
    elseif( H.data.yindex >= s ),
        H.data.yindex = s;
    end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ H ] = updateDATA( H )

    % Only update if mouse pointer is in the
    % axis limits
    set(H.gui,'units','normalized');
    axpos = get(H.axis,'position');
    figcp = get(H.gui,'Currentpoint');
    axlim = axpos(1) + axpos(3);
    aylim = axpos(2) + axpos(4);
    if or(figcp(1) > (axlim+.01), figcp(1) < (axpos(1)-.01)),
        return;
    elseif or(figcp(2) > (aylim+.01), figcp(2) < (axpos(2)-.01)),
        return;
    end
    
    CurrentPoint  = get(H.axis,'Currentpoint');
    H.data.xpoint = CurrentPoint(1,1);
    H.data.ypoint = CurrentPoint(1,2);
    
    doNearTrace   = get(H.traceNearest,'Value');
    
    if (doNearTrace > 0)
        
        % Get new yindex for nearest trace
        [ H.data.xpoint, ...
          H.data.xindex, ...
          H.data.ypoint, ...
          H.data.yindex ] = NearestXYMatrixPoint( H.data.xdata,...
                                                  H.data.ydata,...
                                                  H.data.xpoint,...
                                                  H.data.ypoint);
    else
        H.data.yindex = get(H.trace,'Value');
    end
    
    CurrentPoint  = get(H.axis,'Currentpoint');
    H.data.xpoint = CurrentPoint(1,1);
    H.data.ypoint = CurrentPoint(1,2);
    
    H = interpY(H);
    H = updateGUI(H);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ H ] = updateGUI( H )

    InterpMethod = get(H.interp,'Value');
    if (InterpMethod > 1)
        % There is no specific matrix x-index for 
        % an interpolated point, but the nearest xindex
        % is always returned from the interp function below
        % so that next/prev move function works correctly
        set(H.xindex,'String','interp');
    else
        set(H.xindex,'String',num2str(H.data.xindex));
    end
    set(H.xindex,'Value',H.data.xindex);
    
    tracestr = sprintf('%d',H.data.yindex);
    set(H.yindex,'String',tracestr,'Value',uint16(H.data.yindex));
    set(H.trace,'Value',uint16(H.data.yindex));
    
    % Create the crosshair lines on the figure, crossing at the x,y point
    x_rng  = get(H.axis,'Xlim');
    y_rng  = get(H.axis,'Ylim');
    set(H.xline,'Xdata',[H.data.xpoint H.data.xpoint],'Ydata',y_rng);
    set(H.yline,'Ydata',[H.data.ypoint H.data.ypoint],'Xdata',x_rng);
    
    % Update the x,y values displayed for the x,y point
    xstring = sprintf('%g',H.data.xpoint);
    ystring = sprintf('%g',H.data.ypoint);
    set(H.xvalue,'String',xstring,'Value',H.data.xpoint);
    set(H.yvalue,'String',ystring,'Value',H.data.ypoint);
    
    % Calculate the f(x) function
    fyeq = get(H.fyeq,'string');
    fyval = eval(strrep(fyeq,'y',ystring));
    fystr = sprintf('%g',fyval);
    set(H.fyvalue,'String',fystr,'Value',fyval);
    
    % Calculate the f(y) function
    fxeq = get(H.fxeq,'string');
    fxval = eval(strrep(fxeq,'x',xstring));
    fxstr = sprintf('%g',fxval);
    set(H.fxvalue,'String',fxstr,'Value',fxval);
    
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ H ] = interpY( H )
    
    % In this function, xdata & ydata are arrays, not matrices
    xdata = H.data.xdata(:,H.data.yindex);
    ydata = H.data.ydata(:,H.data.yindex);
    
    if      H.data.xpoint >= max(xdata)                
            H.data.xpoint  = max(xdata);
            H.data.xindex  = find(xdata == max(xdata));
            H.data.ypoint  = ydata(H.data.xindex);
            return;
    elseif  H.data.xpoint <= min(xdata)
            H.data.xpoint  = min(xdata);
            H.data.xindex  = find(xdata == min(xdata));
            H.data.ypoint  = ydata(H.data.xindex);
            return;
    end
    
    % 'none|nearest|linear|spline|cubic'
    interp = get(H.interp,'Value');
    
    switch interp
    case 1
        % Given that xdata & ydata are the same length arrays,
        % we can find the ypoint given the nearest xpoint.
        [H.data.xindex, H.data.xpoint] = NearestXYArrayPoint( xdata, H.data.xpoint );
        H.data.ypoint = ydata(H.data.xindex);
    case 2
        H.data.ypoint = interp1( xdata, ydata, H.data.xpoint, 'nearest' );
    case 3
        H.data.ypoint = interp1( xdata, ydata, H.data.xpoint, 'linear' );
    case 4
        H.data.ypoint = interp1( xdata, ydata, H.data.xpoint, 'spline' );
    case 5
        H.data.ypoint = interp1( xdata, ydata, H.data.xpoint, 'cubic' );
    otherwise
        %use default (linear in matlabR12)
        H.data.ypoint = interp1( xdata, ydata, H.data.xpoint );
    end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Yall ] = interpYall( H )
    
    xdata = H.data.xdata(:,H.data.yindex);
    Yall  = H.data.ydata;
    
    if      H.data.xpoint >= max(xdata),
            H.data.xpoint  = max(xdata);
            H.data.xindex  = find(xdata == max(xdata));
            Yall = ydata(:,H.data.xindex);
            return;
    elseif  H.data.xpoint <= min(xdata),
            H.data.xpoint  = min(xdata);
            H.data.xindex  = find(xdata == min(xdata));
            Yall = ydata(:,H.data.xindex);
            return;
    end
    
    % 'none|nearest|linear|spline|cubic'
    interp = get(H.interp,'Value');
    
    switch interp,
    case 1
        % do nothing in this case
    case 2
        Yall = interp1( xdata, Yall, H.data.xpoint, 'nearest' );
    case 3
        Yall = interp1( xdata, Yall, H.data.xpoint, 'linear' );
    case 4
        Yall = interp1( xdata, Yall, H.data.xpoint, 'spline' );
    case 5
        Yall = interp1( xdata, Yall, H.data.xpoint, 'cubic' );
    otherwise
        %use default (linear in matlabR12)
        Yall = interp1( xdata, Yall, H.data.xpoint );
    end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ index, point ] = NearestXYArrayPoint( data_array, point, type )
    
    if ~exist('type','var') type = ''; end

    % In this function, input data_array is an array, not a matrix.
    % This function returns the data point in the array
    % that has the closest value to the value given (point).  In
    % the context of 'crosshair' the point is a mouse position.
    
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
    
    if strcmp(type,'nextx'),
        index = greater_index;
    elseif strcmp(type,'prevx'),
        index = lesser_index;
    else
        if (greater_dif < lesser_dif)
            index = find(data_array == data_sorted(greater_index));
        else
            index = find(data_array == data_sorted(lesser_index));
        end
    end
    point = data_array(index);
    
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ xpoint, xindex, ypoint, yindex ] = NearestXYMatrixPoint( Xdata, Ydata, xpoint, ypoint )

    % In this function, Xdata & Ydata are matrices of the same dimensions.
    % This function attempts to find the nearest values in Xdata & Ydata
    % to the mouse position (xpoint, ypoint).
    
    % It is assumed that Xdata has identical columns, so we only really
    % need the first column to find the nearest value to xpoint.
    
    [ xindex, xpoint ] = NearestXYArrayPoint( Xdata(:,1), xpoint );
    
    % Now, given the xpoint, we can select just that row of the
    % Ydata matrix corresponding to the xpoint.
    ydata = Ydata(xindex,:);
    
    % The ydata array is searched in same manner as the xdata
    % array for the nearest value.
    [ yindex, ypoint ] = NearestXYArrayPoint( ydata, ypoint );
    
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H] = INIT

    H.gui  = gcf; % Get current figure handles
    H.axis = gca; % Get current axis handles
    
    H.axis = gca; % Get current axis handles
    H.ylim = get(H.axis,'ylim');
    
    % store current button fcn
    H.button = get(H.gui,'WindowButtonDownFcn');
    % set XHR button down fcn
    set(H.gui,'WindowButtonDownFcn','crosshair(''down'');');
    set(H.gui,'KeyPressFcn','crosshair(''keypress'');');
    
    % Match background figure colour
    bgcolor = get(H.gui,'Color');
    % Try to adapt the foreground colour a little
    black = find(bgcolor <= .6);
    fgcolor = [0 0 0]; %black text
    if length(black)>2, fgcolor = [1 1 1]; end
    
    
    Font.FontName   = 'Helvetica';
    Font.FontUnits  = 'Pixels';
    Font.FontSize   = 10;
    Font.FontWeight = 'normal';
    Font.FontAngle  = 'normal';
    
    H.yflip       =   uicontrol(H.gui,'Style','pushbutton','Units','Normalized',Font,...
                              'Position',[.00 .70 .08 .05],...
                              'Tag','YFLIP',...
                              'TooltipString','Flip Y Axis', ...
                              'String','Flip',...
                              'Callback',strcat('XHR = get(gcbf,''userdata'');',...
                                                'ydir = get(XHR.axis,''YDir'');',...
                                                'if isequal(ydir,''normal''),',...
                                                '   set(XHR.axis,''YDir'',''reverse'');',...
                                                'else,',...
                                                '   set(XHR.axis,''YDir'',''normal'');',...
                                                'end;',...
                                                'set(gcbf,''userdata'',XHR); figure(XHR.gui); clear XHR;'));
    
    H.yreset      =   uicontrol(H.gui,'Style','pushbutton','Units','Normalized',Font,...
                              'Position',[.00 .65 .08 .05],...
                              'Tag','YRESET',...
                              'TooltipString','Reset Axis Limits', ...
                              'String','Reset',...
                              'Callback',strcat('XHR = get(gcbf,''userdata'');',...
                                                'XHR.ylim = get(XHR.axis,''ylim'');',...
                                                'set(XHR.ymin,''string'',sprintf(''%g'',XHR.ylim(1)));',...
                                                'set(XHR.ymax,''string'',sprintf(''%g'',XHR.ylim(2)));',...
                                                'set(gcbf,''userdata'',XHR); figure(XHR.gui); clear XHR;'));
    
    H.ymin        =   uicontrol(H.gui,'Style','edit','Units','Normalized',Font,...
                              'Position',[.00 .60 .08 .05],...
                              'HorizontalAlign','left',...
                              'Tag','YMIN',...
                              'TooltipString','Set Y min', ...
                              'String',sprintf('%g',H.ylim(1)),...
                              'Callback',strcat('XHR = get(gcbf,''userdata'');',...
                                                'ymin = str2num(get(XHR.ymin,''string''));',...
                                                'XHR.ylim(1) = ymin;',...
                                                'set(XHR.axis,''ylim'',XHR.ylim);',...
                                                'set(XHR.ymin,''string'',sprintf(''%g'',XHR.ylim(1)));',...
                                                'set(gcbf,''userdata'',XHR); figure(XHR.gui); clear XHR ymin;'));
    
    H.ymax        =   uicontrol(H.gui,'Style','edit','Units','Normalized',Font,...
                              'Position',[.00 .55 .08 .05],...
                              'HorizontalAlign','left',...
                              'Tag','YMAX',...
                              'TooltipString','Set Y max', ...
                              'String',sprintf('%g',H.ylim(2)),...
                              'Callback',strcat('XHR = get(gcbf,''userdata'');',...
                                                'ymax = str2num(get(XHR.ymax,''string''));',...
                                                'XHR.ylim(2) = ymax;',...
                                                'set(XHR.axis,''ylim'',XHR.ylim);',...
                                                'set(XHR.ymax,''string'',sprintf(''%g'',XHR.ylim(2)));',...
                                                'set(gcbf,''userdata'',XHR); figure(XHR.gui); clear XHR ymax;'));
    
    H.grid        =   uicontrol(H.gui,'Style','checkbox','Units','Normalized',Font,...
                              'Position',[.00 .50 .08 .05],...
                              'BackgroundColor',bgcolor,'ForegroundColor',fgcolor,...
                              'Tag','GRID',...
                              'TooltipString','Toggle plot grid on/off.', ...
                              'String','grid',...
                              'Callback',strcat('XHR = get(gcbf,''userdata'');',...
                                                'grid(XHR.axis);',...
                                                'set(gcbf,''userdata'',XHR); figure(XHR.gui); clear XHR;'));
    
    H.xvalue      = uicontrol(H.gui,'Style','edit','Units','Normalized',Font,...
                              'Position',[.13 .95 .15 .05],...
                              'BackGroundColor',[ 0 .9 0],'ForeGroundColor',[ 0 0 0],...
                              'Tag','XVALUE',...
                              'TooltipString','X value (Read Only)',...
                              'String',' ');
    H.yvalue      = uicontrol(H.gui,'Style','edit','Units','Normalized',Font,...
                              'Position',[.28 .95 .15 .05],...
                              'BackGroundColor',[ 0 0 .9],'ForeGroundColor',[ 1 1 1],...
                              'Tag','YVALUE',...
                              'TooltipString','Y value (Read Only)',...
                              'String',' ');
    
    
    H.fxeq        = uicontrol('Style','edit','Units','Normalized',...
                              'Position',[.45 .95 .10 .05],...
                              'Tag','FXDEFINE',...
                              'TooltipString','Define f(x) equation here',...
                              'BackGroundColor',[ 0 .9 0],...
                              'ForeGroundColor',[ 0 0 0],'String','x');
    H.fxvalue     = uicontrol('Style','edit','Units','Normalized',...
                              'Position',[.55 .95 .15 .05],...
                              'Tag','FYVALUE',...
                              'TooltipString','f(x) result',...
                              'BackGroundColor',[ 0 .9 0],...
                              'ForeGroundColor',[ 0 0 0],'String',' ');
    
    H.fyeq        = uicontrol('Style','edit','Units','Normalized',...
                              'Position',[.70 .95 .10 .05],...
                              'Tag','FYDEFINE',...
                              'TooltipString','Define f(y) equation here',...
                              'BackGroundColor',[ 0 0 .9],...
                              'ForeGroundColor',[ 1 1 1],'String','y');
    H.fyvalue     = uicontrol('Style','edit','Units','Normalized',...
                              'Position',[.80 .95 .15 .05],...
                              'Tag','FYVALUE',...
                              'TooltipString','f(y) result',...
                              'BackGroundColor',[ 0 0 .9],...
                              'ForeGroundColor',[ 1 1 1],'String',' ');
                              
    H.yindex      = uicontrol('Style','edit','Units','Normalized',Font,...
                              'Position',[.92 .87 .08 .05],...
                              'BackGroundColor',[ 0 0 .9],'ForeGroundColor',[ 1 1 1],...
                              'Tag','YINDEX',...
                              'TooltipString','Enter Y index into plot data matrix.  Same as trace number.',...
                              'String','1',...
                              'Value',1,...
                              'Callback',strcat('H = get(gcbf,''userdata'');',...
                                                'yi = str2num(get(H.yindex,''String''));',...
                                                'H.data.yindex = yi;',...
                                                'set(H.gui,''userdata'',H); clear H yi; ',...
                                                '[Xpoint,Ypoint] = crosshair(''changey'');'));
    H.yprev       = uicontrol('Style','Push','Units','Normalized',Font,...
                              'Position',[.92 .82 .04 .05],...
                              'String','<',...
                              'Tag','YPREV',...
                              'TooltipString','Goto Previous Y Index (channel).',...
                              'CallBack','[Xpoint,Ypoint] = crosshair(''prevy'');');
    H.ynext       = uicontrol('Style','Push','Units','Normalized',Font,...
                              'Position',[.96 .82 .04 .05],...
                              'String','>',...
                              'Tag','YNEXT',...
                              'TooltipString','Goto Next Y Index (channel).',...
                              'CallBack','[Xpoint,Ypoint] = crosshair(''nexty'');');
    H.yLT         = uicontrol('Style','Push','Units','Normalized',Font,...
                              'Position',[.92 .77 .04 .05],...
                              'String','LT',...
                              'Tag','YLT',...
                              'TooltipString','Goto next Y Less Than current Y.',...
                              'CallBack','[Xpoint,Ypoint] = crosshair(''yLT'');');
    H.yGT         = uicontrol('Style','Push','Units','Normalized',Font,...
                              'Position',[.96 .77 .04 .05],...
                              'String','GT',...
                              'Tag','YGT',...
                              'TooltipString','Goto next Y Greater Than current Y.',...
                              'CallBack','[Xpoint,Ypoint] = crosshair(''yGT'');');
    
    H.xindex      = uicontrol('Style','edit','Units','Normalized',Font,...
                              'Position',[.92 .70 .08 .05],...
                              'BackGroundColor',[ 0 .9 0],'ForeGroundColor',[ 0 0 0],...
                              'Tag','XINDEX',...
                              'TooltipString','Enter X index into plot data matrix.  Only available for interpolation = ''none''.',...
                              'String','1',...
                              'Value',1,...
                              'Callback',strcat('XHR = get(gcbf,''userdata'');',...
                                                'xi = str2num(get(XHR.xindex,''String''));',...
                                                'XHR.data.xindex = xi;',...
                                                'set(XHR.xinterp,''value'',0); ',...
                                                'set(XHR.xinterp,''string'',''0''); ',...
                                                'set(XHR.interp, ''value'',1); ',...
                                                'set(XHR.gui,''userdata'',XHR); clear XHR; ',...
                                                '[Xpoint,Ypoint] = crosshair(''changex'');'));
    H.xprev       = uicontrol('Style','Push','Units','Normalized',Font,...
                              'Position',[.92 .65 .04 .05],...
                              'String','<',...
                              'Tag','XPREV',...
                              'TooltipString','Goto Previous X Index (no interpolation).',...
                              'CallBack','[Xpoint,Ypoint] = crosshair(''prevx'');');
    H.xnext       = uicontrol('Style','Push','Units','Normalized',Font,...
                              'Position',[.96 .65 .04 .05],...
                              'String','>',...
                              'Tag','XNEXT',...
                              'TooltipString','Goto Next X Index (no interpolation).',...
                              'CallBack','[Xpoint,Ypoint] = crosshair(''nextx'');');
    
    H.xinterp     = uicontrol('Style','Edit','Units','Normalized',...
                              'Position',[.92 .60 .08 .05],...
                              'String','0',...
                              'Value',0,....
                              'Tag','XINTERP',...
                              'TooltipString','Interpolation X increment (zero = nearest X).',...
                              'Callback',strcat('XHR = get(gcf,''userdata''); ',...
                                                'xint = str2num(get(XHR.xinterp,''string'')); ',...
                                                'set(XHR.xinterp,''value'',xint); figure(XHR.gui); clear XHR; '));
                                            
    interpstr     =   'none|nearest|linear|spline|cubic';
    H.interp      =   uicontrol('Style','popup','Units','Normalized',Font,...
                              'Position',[.92 .55 .08 .05],...
                              'Tag','INTERP',...
                              'TooltipString','INTERP1 methods (none = raw values).', ...
                              'String',interpstr,...
                              'Callback',strcat('XHR = get(gcbf,''userdata'');',...
                                                'xint = get(XHR.xinterp,''Value''); ',...
                                                'if xint == 0, ',...
                                                '   xint = 0.5; ',...
                                                '   set(XHR.xinterp,''value'',xint); ',...
                                                '   set(XHR.xinterp,''string'',num2str(xint)); ',...
                                                '   set(XHR.xindex, ''string'',''interp''); ',...
                                                'end; ',...
                                                'set(XHR.gui,''userdata'',XHR); figure(XHR.gui); clear XHR; '));
    
    H.store       = uicontrol('Style','Push','Units','Normalized',...
                              'Position',[.92 .40 .08 .05],...
                              'String','Store',...
                              'Tag','STORE',...
                              'TooltipString','Store current Channel & XY values into base XYhist array.', ...
                              'CallBack','crosshair(''store'');');
                          
    H.done        = uicontrol('Style','Push','Units','Normalized',...
                              'Position',[.80 .00 .10 .05],...
                              'BackgroundColor',[.8 .5 0],...
                              'ForegroundColor',[1 1 1],...
                              'FontWeight','bold',...
                              'String','Done',...
                              'Tag','DONE',...
                              'TooltipString','Close crosshair', ...
                              'CallBack','crosshair(''done'');');
    H.exit        = uicontrol('Style','Push','Units','Normalized',...
                              'Position',[.90 .00 .10 .05],...
                              'BackgroundColor',[.8 0 0],...
                              'ForegroundColor',[1 1 1],...
                              'FontWeight','bold',...
                              'String','Exit',...
                              'Tag','EXIT',...
                              'TooltipString','Close crosshair and Figure', ...
                              'CallBack','crosshair(''exit'');');
                          
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Get Line Data from Plot
    
    % Lines are referenced as axis children, among other
    % axis children; so first get all axis children
    sibs = get(H.axis,'Children');
    
    % Now search axis children for any line types.
    % Because the columns of the y data matrix in a plot
    % command seem to be reversed in the axis children, 
    % count down from max sibs to the first sib.
    lines = 0;
    H.data.xdata = [];
    H.data.ydata = [];
    H.data.xpoint = [];
    H.data.ypoint = [];
    H.data.xindex = 1;
    H.data.yindex = 1;
    i = max(size(sibs));
    while i >= 1
        if strcmp(get(sibs(i),'Type'),'line')
            
            % OK, found a line among the axis children.
            lines = lines + 1;
            
            % put line data into a column of DATA.xdata|DATA.ydata
            H.data.xdata(:,lines) = get(sibs(i),'XData').';
            H.data.ydata(:,lines) = get(sibs(i),'YData').';
        end
        i = i - 1;
    end
    
    % 'traces' string variable must be in ascending order
    traces  = '';
    i = 1;
    while i <= lines;
        if i < lines
            tracelabel = sprintf('Column %4d|',i);            
        else
            tracelabel = sprintf('Column %4d',i);
        end
        traces = strcat(traces,tracelabel);
        i = i + 1;
    end
    
    % If more than one line, provide GUI for line selection
    
    % Switch off||on Trace Selection GUI
    Vis = 'Off';
    if lines > 1,
        Vis = 'On';
    elseif lines == 0
        error('No lines found in the current plot window\n');
    end
    
    % Create Trace Index GUI
    H.traceLabel   = uicontrol('Style','Edit', 'Units','Normalized',...
                               'Position',[.00 .00 .15 .05],...
                               'Tag','TRACELABEL',...
                               'Visible',Vis,'String','Select Trace :',...
                               'TooltipString','Select trace to follow with crosshairs.');
    H.trace        = uicontrol('Style','Popup','Units','Normalized',...
                               'Position',[.15 .00 .20 .05],...
                               'Tag','TRACESWITCH',...
                               'BackGroundColor','w','String',traces,...
                               'Visible',Vis,...
                               'CallBack','crosshair(''up'');');
    H.traceNearest = uicontrol('Style','checkbox', 'Units','Normalized',...
                               'Position',[.36 .00 .19 .05],...
                               'BackgroundColor',bgcolor,'ForegroundColor',fgcolor,...
                               'Tag','TRACENEAREST',...
                               'Visible',Vis,'String','Nearest Trace','Value',1,...
                               'TooltipString','Trace nearest to mouse click; switch off to keep trace constant.');
    
    
    % Set X,Y cross hair lines
    % Do this after finding all the line axis children
    % to avoid confusing these lines with those of the
    % plot itself (counted above).
    x_rng = get(H.axis,'Xlim');
    y_rng = get(H.axis,'Ylim');
    H.xline = line(x_rng,[y_rng(1) y_rng(1)]);
    H.yline = line(x_rng,[y_rng(1) y_rng(1)]);
    set(H.xline,'Color','r','EraseMode','xor','Tag','XLINE');
    set(H.yline,'Color','r','EraseMode','xor','Tag','YLINE');
    
return
