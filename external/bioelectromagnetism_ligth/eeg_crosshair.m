function [Xpoint,Ypoint] = eeg_crosshair(action,p,parent);

% eeg_crosshair - A gui interface for reading (x,y) values from a plot.
%
%  A set of mouse driven eeg_crosshairs is placed on the current axes,
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
%  Some further help is given in the tool tips of the GUI.
%
%  Usage: x = [1:10]; y(1,:) = sin(x); y(2,:) = cos(x); x2 = x.^2;
%          figure; plot(x2,y); eeg_crosshair
%


% $Revision: 1.1 $ $Date: 2009-04-28 22:13:51 $

%
%  Licence:        GNU GPL, no express or implied warranties
%  History: 03/96, Richard G. Cobb <cobbr@plk.af.mil>
%           08/01, Darren.Weber_at_radiology.ucsf.edu
%                  replaced obsolete 'table1' with 'interp1'; fixed bug
%                  with number of 'traces'; rationalized calculations into
%                  a common subfunction for x,y point calc in 'down','up', 
%                  & 'move' button functions; added option to turn on/off
%                  interpolation and the exit button; simplified updates 
%                  to graphics using global GUI handle structure.
%           11/01, Darren.Weber_at_radiology.ucsf.edu
%                  added tooltips for several GUI handles
%                  added multiple interpolation methods
%                  added GUI for data matrix indices (given no interpolation)
%                  added option to select trace nearest to mouse click point
%                  reversed order of lines in data matrix to be consistent
%                    with the value returned from the nearest trace subfunction
%                  create eeg_crosshair lines after finding all plot lines to
%                    avoid confusing them with the plot lines
%           01/02, Darren.Weber_at_radiology.ucsf.edu
%                  should now work across multiple plot figures, given
%                    that all gui handles and data are now stored in the
%                    plot figure 'userdata' handle.
%                  added functionality to move smoothly from interpolation
%                    back to the delimited data via the "next/previous" 
%                    buttons.
%           02/02, Darren.Weber_at_radiology.ucsf.edu
%                  added functionality to accept general eeg parameter
%                    structure and call eeg_contours for topo mapping.
%                  added functions for eeg peak detection and visualization.
%           04/02, Darren.Weber_at_radiology.ucsf.edu
%                  only update crosshair if cursor is inside axis limits
%                  optimised peak updating logic and next/prev peak method
%                  added option to 'add channels' when using 'view channels'
%                  added LT/GT for yindex movement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('action','var') 
    action = 'init';
elseif isempty(action)
    action = 'init';
end

XHR = get(gcbf,'userdata');

action = lower(action);

% Check for specific keys and assign reasonable actions
if strcmp(action, 'keypress'),
    
    CC = get(XHR.gui,'CurrentCharacter');
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


switch action,
case 'init',
    
    if ~isempty(XHR),
        if isfield(XHR,'data'),
            fprintf('\nWarning...Crosshair already initialised in current figure.\n');
            return;
        end
    end
    
    % Paint GUI
    XHR = INIT;
    
    % Pass input parameters to the figure handles
    if exist('p','var'),
        XHR.p = p;
        %XHR.p = eeg_peaks(XHR.p); % slows init too much
    end
    if exist('parent','var'), XHR.parent.gui = parent; end
    
    % Update and return values
    XHR = updateDATA(XHR);
    
case 'down',    % Mouse Click Down
    
    set(XHR.gui,'WindowButtonMotionFcn','[Xpoint,Ypoint] = eeg_crosshair(''move'');');
    set(XHR.gui,'WindowButtonUpFcn','[Xpoint,Ypoint] = eeg_crosshair(''up'');');
    
    XHR = updateDATA(XHR);
    
case 'move',    % Mouse Drag Motion
    
    XHR = updateDATA(XHR);
    
case 'up',      % Mouse Click Up
    
    set(XHR.gui,'WindowButtonMotionFcn',' ');
    set(XHR.gui,'WindowButtonUpFcn',' ');
    
    XHR = updateDATA(XHR);
    
case {'nextx','prevx','changex','nexty','prevy','changey','ylt','ygt'}, % Change X/Y
    
    XHR = moveXY(XHR,action);
    
case {'nextpeak','prevpeak','nearpeak'}, % Next or Previous Peak
    
    XHR = movepeak(XHR,action);

case {'viewpeaks'}, % View Peaks
    
    XHR = viewpeaks(XHR);

% Store XY values into a history array
case 'store',
    
    Xpoint = get(XHR.handles.xvalue,'Value');
    Ypoint = get(XHR.handles.yvalue,'Value');
    updateXYhistory(XHR.handles);
    
case {'done','exit'},   % Exit eeg_crosshairs GUI
    
    handles = fieldnames(XHR.handles);
    for i=1:length(handles),
        switch handles{i},
        case {'axis','datalines','gui'},
        otherwise,
            h = getfield(XHR.handles,handles{i});
            if ishandle(h), delete(h); end
        end
    end
    
    if strcmp(action,'exit');
        
        close(XHR.gui);
    else
        
        set(XHR.gui,'WindowButtonUpFcn','');
        set(XHR.gui,'WindowButtonMotionFcn','');
        set(XHR.gui,'WindowButtonDownFcn',XHR.handles.button);
        set(XHR.gui,'HandleVisibility','on');
        set(XHR.gui,'MenuBar','figure');
        set(XHR.gui,'userdata',[]);
        refresh(XHR.gui);
    end
    
    if isfield(XHR,'parent'),
        if isfield(XHR.parent,'gui'),
            if ishandle(XHR.parent.gui),
                % Get the userdata from the parent
                parent = get(XHR.parent.gui,'userdata');
                if isfield(parent,'p') & isfield(XHR,'p'),
                    % Update the parent p structure
                    XHR.p.volt.peaks = [];
                    parent.p = XHR.p;
                    set(XHR.parent.gui,'userdata',parent);
                    if isfield(parent,'handles'),
                        if isfield(parent.handles,'EvoltFile'),
                            set(parent.handles.EvoltFile,'String',parent.p.volt.file);
                        end
                        if isfield(parent.handles,'EvoltPath'),
                            set(parent.handles.EvoltPath,'String',parent.p.volt.path);
                        end
                    end
                end
            end
        end
    end
    
    Xpoint = XHR.data.xpoint;
    Ypoint = XHR.data.ypoint;
    clear XHR;
    return;
    
otherwise,
    
end

set(XHR.gui,'userdata',XHR);
Xpoint = XHR.data.xpoint;
Ypoint = XHR.data.ypoint;
if ishandle(XHR.gui),
    figure(XHR.gui);
end

return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateXYhistory(H),
    
    Ch = get(H.yindex,'Value');
    X  = get(H.xvalue,'Value');
    Y  = get(H.yvalue,'Value');
    
    XY.labels = {'Channel','X','Y'};
    XY.data   = [Ch,X,Y];
    
    if evalin('base','exist(''XYhist'',''var'');'),
        XYhist = evalin('base','XYhist');
        if size(XYhist,2) == 3,
            XYhist.data(end+1,:) = XY.data;
            assignin('base','XYhist',XYhist);
        else
            fprintf('\nWarning: creating new XYhist in base workspace.\n\n');
            assignin('base','XYhist',XY);
        end
    else
        assignin('base','XYhist',XY);
    end
    
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ H ] = updateGUI( H )

    InterpMethod = get(H.handles.interp,'Value');
    if (InterpMethod > 1)
        % There is no specific matrix x-index for 
        % an interpolated point, but the nearest xindex 'value'
        % is always returned from the interp function
        % so that next/prev move function works correctly
        set(H.handles.xindex,'String','interp');
    else
        set(H.handles.xindex,'String',num2str(H.data.xindex));
    end
    set(H.handles.xindex,'Value',H.data.xindex);
    
    % Switch on/off views of traces
    if (get(H.handles.traceView,'Value')),
        prevY = get(H.handles.yindex,'Value');
        currY = uint16(H.data.yindex);
        if ~isequal(prevY,currY),
            addchan = get(H.handles.traceAdd,'Value');
            if ~addchan,
                set(H.handles.datalines(prevY),'Visible','off');
            end
            set(H.handles.datalines(currY),'Visible','on');
        else
            set(H.handles.datalines(currY),'Visible','on');
        end
    else
        set(H.handles.datalines,'Visible','on');
    end
    
    tracestr = sprintf('%d',H.data.yindex);
    set(H.handles.yindex,'String',tracestr,'Value',uint16(H.data.yindex));
    %set(H.handles.trace,'Value',uint16(H.data.yindex));
    
    % Create the crosshair lines on the figure, crossing at the x,y point
    x_rng  = get(H.handles.axis,'Xlim');
    y_rng  = get(H.handles.axis,'Ylim');
    set(H.handles.xline,'Xdata',[H.data.xpoint H.data.xpoint],'Ydata',y_rng);
    set(H.handles.yline,'Ydata',[H.data.ypoint H.data.ypoint],'Xdata',x_rng);
    
    % Update the x,y values displayed for the x,y point
    xstring = sprintf('%14.4f',H.data.xpoint);
    ystring = sprintf('%14.4f',H.data.ypoint);
    set(H.handles.xvalue,'String',xstring,'Value',H.data.xpoint);
    set(H.handles.yvalue,'String',ystring,'Value',H.data.ypoint);
    
    % Update the view of the peak scatterplot
    H = viewpeaks(H);
    
    set(H.gui,'userdata',H);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H] = viewpeaks(H)
    
    view = get(H.handles.peaksView,'Value');
    if ~view,
        % switch off any current peak handles from plot
        if isfield(H.handles,'peakhandles'),
            h = H.handles.peakhandles;
            if ishandle(h(1)), delete(h); end
        end
        return;
    end
    
    % Make sure data required for peaks is available
    H = checkpeaks(H);
    
    % delete any current peak handles from plot
    if isfield(H.handles,'peakhandles'),
        h = H.handles.peakhandles;
        if ishandle(h(1)), delete(h); end
    end
    
    elec = H.data.yindex;
    peaksindex = find(H.p.volt.peaks.data(:,elec) ~= 0);
    
    time = H.p.volt.timeArray(peaksindex,elec);
    volt = H.p.volt.data(peaksindex,elec);
    
    figure(H.gui);
    hold on;
    H.handles.peakhandles = scatter(time,volt);
    hold off;
    
    %if view,
    %    set(H.handles.peakhandles,'Visible','on');
    %else
    %    set(H.handles.peakhandles,'Visible','off');
    %end
    
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H] = movepeak(H,move)
    
    % Make sure data required for peaks is available
    H = checkpeaks(H);
    
    elec = H.data.yindex;
    peaksindex = find(H.p.volt.peaks.data(:,elec) ~= 0);
    
    time = H.p.volt.timeArray(peaksindex,elec);
    volt = H.p.volt.data(peaksindex,elec);
    
    % Find nearest peakindex to currentpoint
    x = H.data.xpoint;
    peaks = [time volt];
    nindex = 1; pindex = 1;
    for i=1:length(peaks),
        peak = peaks(i,:);
        if (x > peak(1)),
            pindex = i;
            if (i+1 <= length(peaks)), nindex = i+1; 
            else                       nindex = i;
            end
        elseif (x == peak(1)),
            pindex = i-1;
            if (i-1 > 0), pindex = i-1;
            else          pindex = 1;
            end
            if (i+1 <= length(peaks)), nindex = i+1;
            else                       nindex = i;
            end
        end
    end
    % Identify prev/next peak
    if strcmp(move,'nextpeak'), index = nindex;
    else                        index = pindex;
    end
    
    H.data.xpoint = time(index);
    H.data.ypoint = volt(index);
    
    xdata = H.data.xdata(:,H.data.yindex);
    H.data.xindex = NearestXYArrayPoint( xdata, H.data.xpoint, 'exact' );
    
    set(H.handles.interp,'Value',1);
    H = updateGUI(H);
    
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H] = checkpeaks(H)
    % Make sure peaks data required is available
    if isfield(H,'p'),
        if isfield(H.p.volt,'peaks'),
            if isempty(H.p.volt.peaks),
                H.p.volt.timeArray = H.data.xdata;
                H.p.volt.data = H.data.ydata;
                H.p = eeg_peaks(H.p);
            end
        else
            H.p.volt.timeArray = H.data.xdata;
            H.p.volt.data = H.data.ydata;
            H.p = eeg_peaks(H.p);
        end
    else
        H.p.volt.timeArray = H.data.xdata;
        H.p.volt.data = H.data.ydata;
        H.p = eeg_peaks(H.p);
    end
    
    % Check the volt.timeArray
    if isempty(H.p.volt.timeArray),
        H.p.volt.timeArray = H.data.xdata;
    end
    if ~isequal(size(H.p.volt.timeArray),size(H.p.volt.data)),
        H.p.volt.timeArray = repmat(H.p.volt.timeArray,1,size(H.p.volt.data,2));
    end
    
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H] = moveXY(H,move)
    
    interp  = get(H.handles.interp,'Value');
    xinterp = get(H.handles.xinterp,'Value');
    
    
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
    
    % Ensure that x/y index is within data range
    H = checkdatarange(H);
    
    % Get x/y value at new x/y index
    H.data.xpoint = H.data.xdata(H.data.xindex,H.data.yindex);
    H.data.ypoint = H.data.ydata(H.data.xindex,H.data.yindex);
    
    set(H.handles.interp,'Value',1);
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
    axpos = get(H.handles.axis,'position');
    figcp = get(H.gui,'Currentpoint');
    axlim = axpos(1) + axpos(3);
    aylim = axpos(2) + axpos(4);
    if or(figcp(1) > (axlim+.01),figcp(1) < (axpos(1)-.01)),
        return;
    elseif or(figcp(2) > (aylim+.01),figcp(2) < (axpos(2)-.01)),
        return;
    end
    
    % OK, update...
    CurrentPoint  = get(H.handles.axis,'Currentpoint');
    H.data.xpoint = CurrentPoint(1,1);
    H.data.ypoint = CurrentPoint(1,2);
    
    doNearTrace   = get(H.handles.traceNearest,'Value');
    
    if (doNearTrace > 0)
        
        % Get new yindex for nearest trace
        [ H.data.xpoint, ...
          H.data.xindex, ...
          H.data.ypoint, ...
          H.data.yindex ] = NearestXYMatrixPoint( H.data.xdata,...
                                                  H.data.ydata,...
                                                  H.data.xpoint,...
                                                  H.data.ypoint);
%     else
%         H.data.yindex = get(H.handles.trace,'Value');
    end
    
    CurrentPoint  = get(H.handles.axis,'Currentpoint');
    H.data.xpoint = CurrentPoint(1,1);
    H.data.ypoint = CurrentPoint(1,2);
    
    H = interpY(H);
    H = updateGUI(H);

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
    interp = get(H.handles.interp,'Value');
    
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
    interp = get(H.handles.interp,'Value');
    
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
    % the context of 'eeg_crosshair' the point is a mouse position.
    
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
    lesser_index = lesser(end);
    
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
    
    H.gui = gcf;          % Get current figure handle
    
    % Match background figure colour
    bgcolor = get(H.gui,'Color');
    % Try to adapt the foreground colour a little
    black = find(bgcolor <= .6);
    fgcolor = [0 0 0]; %black text
    if length(black)>2, fgcolor = [1 1 1]; end
    
    
    H.handles.axis = gca; % Get current axis handles
    axis tight;
    set(H.handles.axis,'YDir','reverse');
    H.handles.ylim = get(H.handles.axis,'ylim');
    
    % enable right click access to eeg_crosshair
    if isempty(get(H.handles.axis,'uicontextmenu')),
        amenu=uicontextmenu;
        a=uimenu(amenu,'Label','Crosshair','Callback','eeg_crosshair; ');
        set(H.handles.axis,'uicontextmenu',amenu);
    end
    
    % store current button fcn
    H.handles.button = get(H.gui,'WindowButtonDownFcn');
    % set XHR button down fcn
    set(H.gui,'WindowButtonDownFcn','[Xpoint,Ypoint] = eeg_crosshair(''down'');');
    set(H.gui,'KeyPressFcn','[Xpoint,Ypoint] = eeg_crosshair(''keypress'');');
    
    Font.FontName   = 'Helvetica';
    Font.FontUnits  = 'Pixels';
    Font.FontSize   = 10;
    Font.FontWeight = 'normal';
    Font.FontAngle  = 'normal';
    
    
    H.handles.yflip       =   uicontrol(H.gui,'Style','pushbutton','Units','Normalized',Font,...
                              'Position',[.00 .70 .08 .05],...
                              'Tag','YFLIP',...
                              'TooltipString','Flip Y Axis', ...
                              'String','Flip',...
                              'Callback',strcat('XHR = get(gcbf,''userdata'');',...
                                                'ydir = get(XHR.handles.axis,''YDir'');',...
                                                'if isequal(ydir,''normal''),',...
                                                '   set(XHR.handles.axis,''YDir'',''reverse'');',...
                                                'else,',...
                                                '   set(XHR.handles.axis,''YDir'',''normal'');',...
                                                'end;',...
                                                'set(gcbf,''userdata'',XHR); figure(XHR.gui); clear XHR;'));
    
    H.handles.yreset      =   uicontrol(H.gui,'Style','pushbutton','Units','Normalized',Font,...
                              'Position',[.00 .65 .08 .05],...
                              'Tag','YRESET',...
                              'TooltipString','Reset Axis Limits', ...
                              'String','Reset',...
                              'Callback',strcat('XHR = get(gcbf,''userdata'');',...
                                                'axis tight;',...
                                                'XHR.handles.ylim = get(XHR.handles.axis,''ylim'');',...
                                                'set(XHR.handles.ymin,''string'',sprintf(''%7.1f'',XHR.handles.ylim(1)));',...
                                                'set(XHR.handles.ymax,''string'',sprintf(''%7.1f'',XHR.handles.ylim(2)));',...
                                                'set(gcbf,''userdata'',XHR); figure(XHR.gui); clear XHR;'));
    
    H.handles.ymin        =   uicontrol(H.gui,'Style','edit','Units','Normalized',Font,...
                              'Position',[.00 .60 .08 .05],...
                              'HorizontalAlign','left',...
                              'Tag','YMIN',...
                              'TooltipString','Set Y min', ...
                              'String',sprintf('%7.1f',H.handles.ylim(1)),...
                              'Callback',strcat('XHR = get(gcbf,''userdata'');',...
                                                'ymin = str2num(get(XHR.handles.ymin,''string''));',...
                                                'XHR.handles.ylim(1) = ymin;',...
                                                'set(XHR.handles.axis,''ylim'',XHR.handles.ylim);',...
                                                'set(XHR.handles.ymin,''string'',sprintf(''%7.1f'',XHR.handles.ylim(1)));',...
                                                'set(gcbf,''userdata'',XHR); figure(XHR.gui); clear XHR ymin;'));
    
    H.handles.ymax        =   uicontrol(H.gui,'Style','edit','Units','Normalized',Font,...
                              'Position',[.00 .55 .08 .05],...
                              'HorizontalAlign','left',...
                              'Tag','YMAX',...
                              'TooltipString','Set Y max', ...
                              'String',sprintf('%7.1f',H.handles.ylim(2)),...
                              'Callback',strcat('XHR = get(gcbf,''userdata'');',...
                                                'ymax = str2num(get(XHR.handles.ymax,''string''));',...
                                                'XHR.handles.ylim(2) = ymax;',...
                                                'set(XHR.handles.axis,''ylim'',XHR.handles.ylim);',...
                                                'set(XHR.handles.ymax,''string'',sprintf(''%7.1f'',XHR.handles.ylim(2)));',...
                                                'set(gcbf,''userdata'',XHR); figure(XHR.gui); clear XHR ymax;'));
    
    H.handles.grid        =   uicontrol(H.gui,'Style','checkbox','Units','Normalized',Font,...
                              'Position',[.00 .50 .08 .05],...
                              'BackgroundColor',bgcolor,'ForegroundColor',fgcolor,...
                              'Tag','GRID',...
                              'TooltipString','Toggle plot grid on/off.', ...
                              'String','grid',...
                              'Callback',strcat('XHR = get(gcbf,''userdata'');',...
                                                'grid(XHR.handles.axis);',...
                                                'set(gcbf,''userdata'',XHR); figure(XHR.gui); clear XHR;'));
    
    H.handles.xvalue      = uicontrol(H.gui,'Style','edit','Units','Normalized',Font,...
                              'Position',[.13 .95 .15 .05],...
                              'Tag','XVALUE',...
                              'TooltipString','X value (Read Only)',...
                              'BackGroundColor',[ 0 .7 .7],'String',' ');
    H.handles.yvalue      = uicontrol(H.gui,'Style','edit','Units','Normalized',Font,...
                              'Position',[.28 .95 .15 .05],...
                              'Tag','YVALUE',...
                              'TooltipString','Y value (Read Only)',...
                              'BackGroundColor',[ 0 .7 .7],'String',' ');
    
    H.handles.yindex      = uicontrol(H.gui,'Style','edit','Units','Normalized',Font,...
                              'Position',[.92 .87 .08 .05],...
                              'BackGroundColor',[ 0 .7 .7],...
                              'Tag','YINDEX',...
                              'TooltipString','Enter Y index into plot data matrix.  Same as trace number.',...
                              'String','1',...
                              'Value',1,...
                              'Callback',strcat('XHR = get(gcbf,''userdata'');',...
                                                'yi = str2num(get(XHR.handles.yindex,''String''));',...
                                                'XHR.data.yindex = yi;',...
                                                'set(XHR.gui,''userdata'',XHR); figure(XHR.gui); clear XHR yi; ',...
                                                '[Xpoint,Ypoint] = eeg_crosshair(''changey'');'));
    H.handles.yprev       = uicontrol(H.gui,'Style','Push','Units','Normalized',Font,...
                              'Position',[.92 .82 .04 .05],...
                              'String','<',...
                              'Tag','YPREV',...
                              'TooltipString','Goto Previous Y Index (channel).',...
                              'CallBack','[Xpoint,Ypoint] = eeg_crosshair(''prevy'');');
    H.handles.ynext       = uicontrol(H.gui,'Style','Push','Units','Normalized',Font,...
                              'Position',[.96 .82 .04 .05],...
                              'String','>',...
                              'Tag','YNEXT',...
                              'TooltipString','Goto Next Y Index (channel).',...
                              'CallBack','[Xpoint,Ypoint] = eeg_crosshair(''nexty'');');
    H.handles.yLT         = uicontrol(H.gui,'Style','Push','Units','Normalized',Font,...
                              'Position',[.92 .77 .04 .05],...
                              'String','LT',...
                              'Tag','YLT',...
                              'TooltipString','Goto next Y Less Than current Y.',...
                              'CallBack','[Xpoint,Ypoint] = eeg_crosshair(''yLT'');');
    H.handles.yGT         = uicontrol(H.gui,'Style','Push','Units','Normalized',Font,...
                              'Position',[.96 .77 .04 .05],...
                              'String','GT',...
                              'Tag','YGT',...
                              'TooltipString','Goto next Y Greater Than current Y.',...
                              'CallBack','[Xpoint,Ypoint] = eeg_crosshair(''yGT'');');
    
    H.handles.xindex      = uicontrol(H.gui,'Style','edit','Units','Normalized',Font,...
                              'Position',[.92 .70 .08 .05],...
                              'BackGroundColor',[ 0 .7 .7],...
                              'Tag','XINDEX',...
                              'TooltipString','Enter X index into plot data matrix.  Only available for interpolation = ''none''.',...
                              'String','1',...
                              'Value',1,...
                              'Callback',strcat('XHR = get(gcbf,''userdata'');',...
                                                'xi = str2num(get(XHR.handles.xindex,''String''));',...
                                                'XHR.data.xindex = xi;',...
                                                'set(XHR.handles.xinterp,''value'',0); ',...
                                                'set(XHR.handles.xinterp,''string'',''0''); ',...
                                                'set(XHR.handles.interp, ''value'',1); ',...
                                                'set(XHR.gui,''userdata'',XHR); figure(XHR.gui); clear XHR; ',...
                                                '[Xpoint,Ypoint] = eeg_crosshair(''changex'');'));
    H.handles.xprev       = uicontrol(H.gui,'Style','Push','Units','Normalized',Font,...
                              'Position',[.92 .65 .04 .05],...
                              'String','<',...
                              'Tag','XPREV',...
                              'TooltipString','Goto Previous X Index (no interpolation).',...
                              'CallBack','[Xpoint,Ypoint] = eeg_crosshair(''prevx'');');
    H.handles.xnext       = uicontrol(H.gui,'Style','Push','Units','Normalized',Font,...
                              'Position',[.96 .65 .04 .05],...
                              'String','>',...
                              'Tag','XNEXT',...
                              'TooltipString','Goto Next X Index (no interpolation).',...
                              'CallBack','[Xpoint,Ypoint] = eeg_crosshair(''nextx'');');
    H.handles.xinterp     = uicontrol('Style','Edit','Units','Normalized',...
                              'Position',[.92 .60 .08 .05],...
                              'String','0',...
                              'Value',0,....
                              'Tag','XINTERP',...
                              'TooltipString','Interpolation X increment (zero = nearest X).',...
                              'Callback',strcat('XHR = get(gcbf,''userdata''); ',...
                                                'xint = str2num(get(XHR.handles.xinterp,''string'')); ',...
                                                'set(XHR.handles.xinterp,''value'',xint); figure(XHR.gui); clear XHR;'));
    
    interpstr     =   'none|nearest|linear|spline|cubic';
    H.handles.interp      =   uicontrol(H.gui,'Style','popup','Units','Normalized',Font,...
                              'Position',[.92 .55 .08 .05],...
                              'Tag','INTERP',...
                              'TooltipString','INTERP1 methods (none = raw values).', ...
                              'String',interpstr,...
                              'Callback',strcat('XHR = get(gcbf,''userdata'');',...
                                                'xint = get(XHR.handles.xinterp,''Value''); ',...
                                                'if xint == 0, ',...
                                                '   xint = 0.5; ',...
                                                '   set(XHR.handles.xinterp,''value'',xint); ',...
                                                '   set(XHR.handles.xinterp,''string'',num2str(xint)); ',...
                                                '   set(XHR.handles.xindex, ''string'',''interp''); ',...
                                                'end; ',...
                                                'set(XHR.gui,''userdata'',XHR); figure(XHR.gui); clear XHR; '));
    
    Font.FontWeight = 'bold';
    H.handles.store       = uicontrol('Style','Push','Units','Normalized',Font,...
                              'Position',[.92 .45 .08 .05],...
                              'BackgroundColor',[0 0.7 0],'ForegroundColor', [1 1 1],...
                              'String','Store',...
                              'Tag','STORE',...
                              'TooltipString','Store current values into base XYhist array.', ...
                              'CallBack','eeg_crosshair(''store'');');
    Font.FontWeight = 'normal';
    
    H.handles.topo        = uicontrol(H.gui,'Style','Push','Units','Normalized',Font,...
                              'Position',[.92 .35 .08 .05],...
                              'String','Topo',...
                              'Tag','TOPO',...
                              'TooltipString','Show Topographic Map (or GUI)', ...
                              'CallBack',strcat('XHR = get(gcbf,''userdata'');',...
                                                'if isfield(XHR,''p''), ',...
                                                    'XHR.p.clickTimePoint = 0;',...
                                                    'XHR.p.volt.samplePoint = XHR.data.xindex;',...
                                                    'XHR.p.volt.sampleTime = XHR.data.xpoint;',...
                                                    'if get(XHR.handles.topo_opt,''Value''),',...
                                                        '[G,XHR.p] = gui_eeg_contours(XHR.p,XHR.gui);',...
                                                        'clear G; ',...
                                                    'else, ',...
                                                        'XHR.p = eeg_contours_engine(XHR.p);',...
                                                    'end;',...
                                                'else ',...
                                                    '[G,XHR.p] = gui_eeg_contours('''',XHR.gui);',...
                                                    'clear G; ',...
                                                'end;[p] = XHR.p;',...
                                                'set(XHR.gui,''userdata'',XHR); clear XHR;'));                                            
    H.handles.topo_opt    = uicontrol(H.gui,'Style','checkbox', 'Units','Normalized',Font,...
                              'Position',[.92 .30 .08 .05],...
                              'BackgroundColor',bgcolor,'ForegroundColor',fgcolor,...
                              'Tag','TOPO_OPT',...
                              'String','GUI','Value',0,...
                              'TooltipString','View/Modify Topography Parameters');
                          
    H.handles.peaks       = uicontrol(H.gui,'Style','Push','Units','Normalized',Font,...
                              'Position',[.92 .25 .08 .05],...
                              'String','Peaks',...
                              'Tag','PEAKS',...
                              'TooltipString','Calculate Peaks for All Channels (RunOnce)', ...
                              'CallBack',strcat('XHR = get(gcbf,''userdata'');',...
                                                'if isfield(XHR,''p''), ',...
                                                    'XHR.p = eeg_peaks(XHR.p);',...
                                                    'set(XHR.handles.peaksView,''Value'',1); ',...
                                                'else, ',...
                                                    'XHR.p.volt.data = XHR.data.ydata;',...
                                                    'XHR.p.volt.timeArray = XHR.data.xdata;',...
                                                    'XHR.p = eeg_peaks(XHR.p);',...
                                                    'set(XHR.handles.peaksView,''Value'',1); ',...
                                                'end; ',...
                                                'set(XHR.gui,''userdata'',XHR); clear XHR;',...
                                                '[Xpoint,Ypoint] = eeg_crosshair(''viewpeaks''); '));
    H.handles.pprev       = uicontrol(H.gui,'Style','Push','Units','Normalized',Font,...
                              'Position',[.92 .20 .04 .05],...
                              'String','<',...
                              'Tag','PPREV',...
                              'TooltipString','Goto Previous Peak of Current Channel',...
                              'CallBack','[Xpoint,Ypoint] = eeg_crosshair(''prevpeak'');');                          
    H.handles.pnext       = uicontrol(H.gui,'Style','Push','Units','Normalized',Font,...
                              'Position',[.96 .20 .04 .05],...
                              'String','>',...
                              'Tag','PNEXT',...
                              'TooltipString','Goto Next Peak of Current Channel',...
                              'CallBack','[Xpoint,Ypoint] = eeg_crosshair(''nextpeak'');');
    H.handles.peaksView   = uicontrol(H.gui,'Style','checkbox','Units','Normalized',Font,...
                               'Position',[.92 .15 .08 .05],...
                               'BackgroundColor',bgcolor,'ForegroundColor',fgcolor,...
                               'String','View','Value',0,...
                               'Tag','VIEWPEAK',...
                               'TooltipString','View peaks for current channel.',...
                               'CallBack','[Xpoint,Ypoint] = eeg_crosshair(''viewpeaks'');');
    H.handles.peaksMove   = uicontrol(H.gui,'Style','checkbox','Units','Normalized',Font,...
                               'Position',[.92 .10 .08 .05],...
                               'BackgroundColor',bgcolor,'ForegroundColor',fgcolor,...
                               'String','Move','Value',0,...
                               'Tag','MOVEPEAK','Visible','off',...
                               'TooltipString','Auto move to nearest peak for current channel.',...
                               'CallBack',strcat('[Xpoint,Ypoint] = eeg_crosshair(''nearpeak'');'));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Get Line Data from Plot
    
    % Lines are referenced as axis children, among other
    % axis children; so first get all axis children
    sibs = get(H.handles.axis,'Children');
    
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
            % Found a line child, but check its size
            getline = 1;
            if ~isempty(H.data.xdata),
                if isequal(size(get(sibs(i),'XData').',1),size(H.data.xdata,1)),
                    getline = 1;
                else
                    getline = 0;
                end
            end
            if getline,
                % OK, found a line among the axis children.
                lines = lines + 1;
                H.handles.datalines(lines) = sibs(i);
                
                % put line data into a column of data.xdata|data.ydata
                H.data.xdata(:,lines) = get(sibs(i),'XData').';
                H.data.ydata(:,lines) = get(sibs(i),'YData').';
            end
        end
        i = i - 1;
    end
    
    % Switch off||on Trace GUI controls
    Vis = 'Off';
    if lines > 1,
        Vis = 'On';
    elseif lines == 0,
        error('No lines found in the current plot window\n');
    end
    
    % REMOVING THE FOLLOWING, AS NOW EFFECTIVELY AVAILABLE FROM YINDEX edit box
    
%     % 'traces' string variable must be in ascending order
%     traces  = '';
%     i = 1;
%     while i <= lines;
%         if i < lines
%             tracelabel = sprintf('Channel %4d|',i);            
%         else
%             tracelabel = sprintf('Channel %4d',i);
%         end
%         traces = strcat(traces,tracelabel);
%         i = i + 1;
%     end
%     
%     % If more than one line, provide GUI for line selection
%     
%     
%     % Create Trace Index GUI
%     H.handles.traceLabel   = uicontrol(H.gui,'Style','Edit', 'Units','Normalized',Font,...
%                                'Position',[.00 .00 .09 .05],...
%                                'Tag','TRACELABEL',...
%                                'Visible',Vis,'String','Select :',...
%                                'TooltipString','Select channel to follow with crosshairs.');
%     H.handles.trace        = uicontrol(H.gui,'Style','Popup','Units','Normalized',Font,...
%                                'Position',[.10 .00 .19 .05],...
%                                'Tag','TRACESWITCH',...
%                                'BackGroundColor','w','String',traces,...
%                                'Visible',Vis,...
%                                'CallBack',strcat('XHR = get(gcbf,''userdata'');',...
%                                                  'XHR.data.yindex = get(XHR.handles.trace,''Value'');',...
%                                                  'set(XHR.gui,''userdata'',XHR); clear XHR; ',...
%                                                  '[Xpoint,Ypoint] = eeg_crosshair(''changey'');'));
    
    
    H.handles.traceNearest = uicontrol(H.gui,'Style','checkbox', 'Units','Normalized',Font,...
                               'Position',[.30 .00 .14 .05],...
                               'BackgroundColor',bgcolor,'ForegroundColor',fgcolor,...
                               'Tag','TRACENEAREST',...
                               'Visible',Vis,'String','Nearest','Value',1,...
                               'TooltipString','Channel nearest to mouse click selection. Switch off for fixed channel.');
    H.handles.traceView    = uicontrol(H.gui,'Style','checkbox','Units','Normalized',Font,...
                               'Position',[.45 .00 .19 .05],...
                               'BackgroundColor',bgcolor,'ForegroundColor',fgcolor,...
                               'String','View Channel','Value',0,...
                               'Tag','VIEW','Visible',Vis,...
                               'TooltipString','View only selected channel (on) or all channels (off).', ...
                               'CallBack',strcat('XHR = get(gcbf,''userdata'');',...
                                                 'if (get(XHR.handles.traceView,''Value'')),',...
                                                     'set(XHR.handles.datalines,''Visible'',''off''); ',...
                                                     'set(XHR.handles.datalines(XHR.data.yindex),''Visible'',''on''); ',...
                                                 'else, ',...
                                                     'set(XHR.handles.datalines,''Visible'',''on''); ',...
                                                 'end; ',...
                                                 'figure(XHR.gui); clear XHR;'));
    H.handles.traceAdd    = uicontrol(H.gui,'Style','checkbox', 'Units','Normalized',Font,...
                               'Position',[.65 .00 .19 .05],...
                               'BackgroundColor',bgcolor,'ForegroundColor',fgcolor,...
                               'Tag','TRACEADD',...
                               'Visible',Vis,'String','Add Channel','Value',0,...
                               'TooltipString','Add channels when using ''view channel'' rather than replace them.');
    
    Font.FontWeight = 'bold';
    
    % ASCII Parameters
    H.handles.Bparam      = uicontrol(H.gui,'Style','pushbutton','Units','Normalized',Font,...
                              'Position',[.69 .95 .1 .05],'BusyAction','queue',...
                              'String','PARAM',...
                              'TooltipString','Define ERP sample rate, epoch, etc.',...
                              'BackgroundColor',[0.0 0.0 0.75],...
                              'ForegroundColor',[1 1 1], 'HorizontalAlignment', 'center',...
                              'Callback',strcat('XHR = get(gcbf,''userdata'');',...
                                                'tempgui = gui_eeg_ascii_parameters(XHR.gui); clear tempgui; ',...
                                                'XHR = get(gcbf,''userdata'');',...
                                                'figure(XHR.gui); [Xpoint,Ypoint] = eeg_crosshair(''done''); clf;',...
                                                'plot(XHR.p.volt.timeArray,XHR.p.volt.data);',...
                                                '[Xpoint,Ypoint] = eeg_crosshair(''init'',XHR.p,XHR.parent.gui); clear XHR; '));
    % Interface to electrodes
    H.handles.Belec       = uicontrol(H.gui,'Style','pushbutton','Units','Normalized',Font,...
                              'Position',[.79 .95 .1 .05],...
                              'String','ELEC','HorizontalAlignment','center',...
                              'TooltipString','Load associated electrode coordinates.',...
                              'BackgroundColor',[0.0 0.0 0.75],'ForegroundColor', [1 1 1],...
                              'Callback',strcat('XHR = get(gcbf,''userdata'');',...
                                                'p = gui_elec_open(XHR.p,''init'',XHR.gui);',...
                                                'clear XHR;'));
    % Interface to meshes
    H.handles.Bmesh       = uicontrol(H.gui,'Style','pushbutton','Units','Normalized',Font,...
                              'Position',[.89 .95 .1 .05],...
                              'String','MESH','HorizontalAlignment','center',...
                              'visible','on',...
                              'TooltipString','Load associated MRI tesselations.',...
                              'BackgroundColor',[0.0 0.0 0.75],'ForegroundColor', [1 1 1],...
                              'Callback',strcat('XHR = get(gcbf,''userdata'');',...
                                                'XHR.p = gui_mesh_open(XHR.p,''init'',XHR.gui);',...
                                                'clear XHR;'));
    
    
    H.handles.done        = uicontrol(H.gui,'Style','Push','Units','Normalized',Font,...
                              'Position',[.84 .00 .08 .05],...
                              'BackgroundColor',[0.7 0 0],'ForegroundColor', [1 1 1],...
                              'String','Done',...
                              'Tag','DONE',...
                              'TooltipString','Close eeg_crosshair', ...
                              'CallBack','[Xpoint,Ypoint] = eeg_crosshair(''done'');');
    H.handles.exit        = uicontrol(H.gui,'Style','Push','Units','Normalized',Font,...
                              'Position',[.92 .00 .08 .05],...
                              'BackgroundColor',[0.8 0 0],'ForegroundColor', [1 1 1],...
                              'String','Exit',...
                              'Tag','EXIT',...
                              'TooltipString','Close eeg_crosshair and Figure', ...
                              'CallBack','[Xpoint,Ypoint] = eeg_crosshair(''exit'');');
                          

    % Set X,Y cross hair lines
    % Do this after finding all the line axis children
    % to avoid confusing these lines with those of the
    % plot itself (counted above).
    x_rng = get(H.handles.axis,'Xlim');
    y_rng = get(H.handles.axis,'Ylim');
    axes(H.handles.axis);
    H.handles.xline = line(x_rng,[y_rng(1) y_rng(1)]);
    H.handles.yline = line(x_rng,[y_rng(1) y_rng(1)]);
    set(H.handles.xline,'Color','r','EraseMode','xor','Tag','XLINE');
    set(H.handles.yline,'Color','r','EraseMode','xor','Tag','YLINE');
    
    %set(H.gui,'HandleVisibility','callback'); % hide fig from other commands
    
return
