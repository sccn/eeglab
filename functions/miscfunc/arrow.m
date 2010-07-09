
% arrow() - Draw a line with an arrowhead.
%
% Usage:
%  >> arrow('Property1',PropVal1,'Property2',PropVal2,...)
%  >> arrow(H,'Prop1',PropVal1,...)
%  >> arrow(Start,Stop)
%  >> arrow(Start,Stop,Length,BaseAngle,TipAngle,Width,Page,CrossDir)
%  >> arrow demo  %2-D demos of the capabilities of arrow()
%  >> arrow demo2 %3-D demos of the capabilities of arrow()
%
% Inputs:
%    H               - vector of handles to previously created arrows and/or 
%                      line objects, will update the previously-created
%                      arrows according to the current view and any specified 
%                      properties, and will convert two-point line objects to 
%                      corresponding arrows. Note that arrow(H) will update the 
%                      arrows if the current view has changed. 
%    Start           - vectors of length 2 or 3, or matrices with 2 or 3 columns
%    Stop            - same as Start
%    'Start'         - The starting points.                        B
%    'Stop'          - The end points.                            /|\           ^
%    'Length'        - Length of the arrowhead in pixels.        /|||\          |
%    'BaseAngle'     - Base angle in degrees (ADE).             //|||\\        L|
%    'TipAngle'      - Tip angle in degrees (ABC).             ///|||\\\       e|
%    'Width'         - Width of the base in pixels.           ////|||\\\\      n|
%    'Page'          - Use hardcopy proportions.             /////|D|\\\\\     g|
%    'CrossDir'      - Vector || to arrowhead plane.        ////  |||  \\\\    t|
%    'NormalDir'     - Vector out of arrowhead plane.      ///    |||    \\\   h|
%    'Ends'          - Which end has an arrowhead.        //<----->||      \\   |
%    'ObjectHandles' - Vector of handles to update.      /   base |||        \  V
%    'LineStyle'     - The linestyle of the arrow.      E    angle||<-------->C
%    'LineWidth'     - Line thicknesses.                          |||tipangle
%    'FaceColor'     - FaceColor of patch arrows.                 |||
%    'EdgeColor'     - EdgeColor/Color of patch/line arrows.      |||
%    'Color'         - Set FaceColor & EdgeColor properties.   -->|A|<-- width
%
% Notes:
%  A property list can follow any specified normal argument list, e.g.,
%  ARROW([1 2 3],[0 0 0],36,'BaseAngle',60) creates an arrow from (1,2,3) to
%  the origin, with an arrowhead of length 36 pixels and 60-degree base angle.
%
%  The basic arguments or properties can generally be vectorized to create
%  multiple arrows with the same call.  This is done by passing a property
%  with one row per arrow, or, if all arrows are to have the same property
%  value, just one row may be specified.
%
%  You may want to execute AXIS(AXIS) before calling ARROW so it doesn't change
%  the axes on you; ARROW determines the sizes of arrow components BEFORE the
%  arrow is plotted, so if ARROW changes axis limits, arrows may be malformed.
%
%  ARROW uses features of Matlab 4.2c and later, so earlier versions may be
%  incompatible; call ARROW VERSION for more details.
%
% Author: Erik A. Johnson <johnsone@uxh.cso.uiuc.edu>, 1995 

% Copyright (c)1995, Erik A. Johnson <johnsone@uxh.cso.uiuc.edu>, 10/6/95

% Many thanks to Keith Rogers <kerog@ai.mit.edu> for his many excellent
% suggestions and beta testing.  Check out his shareware package MATDRAW.
% He has permission to distribute ARROW with MATDRAW.


% $Log: arrow.m,v $
% Revision 1.5  2009/10/20 22:19:56  dev
%
% added log to file
%


function [h,yy,zz] = arrow(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8, ...
                           arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16, ...
                           arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24)

% Are we doing the demo?
c = sprintf('\n');
if (nargin==1),
	if (ischar(arg1)),
		arg1 = lower([arg1 '                ']);
		if (strcmp(arg1(1:4),'prop')),
			disp([c ...
			      'ARROW Properties:  Default values are given in [square brackets], and other' c ...
			      '                   acceptable equivalent property names are in (parenthesis).' c c ...
			      '  Start           The starting points. For N arrows,            B' c ...
			      '                  this should be a Nx2 or Nx3 matrix.          /|\           ^' c ...
			      '  Stop            The end points. For N arrows, this          /|||\          |' c ...
			      '                  should be a Nx2 or Nx3 matrix.             //|||\\        L|' c ...
			      '  Length          Length of the arrowhead (in pixels on     ///|||\\\       e|' c ...
			      '                  screen, points on a page). [16] (Len)    ////|||\\\\      n|' c ...
			      '  BaseAngle       Angle (degrees) of the base angle       /////|D|\\\\\     g|' c ...
			      '                  ADE.  For a simple stick arrow, use    ////  |||  \\\\    t|' c ...
			      '                  BaseAngle=TipAngle. [90] (Base)       ///    |||    \\\   h|' c ...
			      '  TipAngle        Angle (degrees) of tip angle ABC.    //<----->||      \\   |' c ...
			      '                  [16] (Tip)                          /   base |||        \  V' c ...
			      '  Width           Width of the base in pixels.  Not  E   angle ||<-------->C' c ...
			      '                  the ''LineWidth'' prop. [0] (Wid)              |||tipangle' c ...
			      '  Page            If provided, non-empty, and not NaN,         |||' c ...
			      '                  this causes ARROW to use hardcopy            |||' c ...
			      '                  rather than onscreen proportions.             A' c ...
			      '                  This is important if screen aspect        -->   <-- width' c ...
			      '                  ratio and hardcopy aspect ratio are    ----CrossDir---->' c ...
				  '                  vastly different. []' c...
			      '  CrossDir        A vector giving the direction towards which the fletches' c ...
			      '                  on the arrow should go.  [computed such that it is perpen-' c ...
			      '                  dicular to both the arrow direction and the view direction' c ...
			      '                  (i.e., as if it was pasted on a normal 2-D graph)]  (Note' c ...
			      '                  that CrossDir is a vector; if an axis is plotted on a log' c ...
			      '                  scale, then the corresponding component of CrossDir must' c ...
			      '                  also be set appropriately, i.e., to 1 for no change in' c ...
			      '                  that direction, >1 for a positive change, >0 and <1 for' c ...
			      '                  negative change.)' c ...
			      '  NormalDir       A vector normal to the fletch direction (CrossDir is then' c ...
			      '                  computed by the vector cross product {Line}x{NormalDir}). []' c ...
				  '  Ends            Set which end has an arrowhead.  Valid values are ''none'',' c ...
				  '                  ''stop'', ''start'', and ''both''. [''stop''] (End)' c...
			      '  ObjectHandles   Vector of handles to previously-created arrows to be' c ...
			      '                  updated or line objects to be converted to arrows.' c ...
			      '                  [] (Object,Handle)' c ...
			      '  LineStyle       The linestyle of the arrow.  If anything other than ''-'',' c ...
			      '                  the arrow will be drawn with a line object, otherwise it' c ...
			      '                  will be drawn with a patch.  [''-''] (LineS)' c ...
			      '  LineWidth       Same as used in SET commands, but may be passed a vector' c ...
			      '                  for multiple arrows. [default of the line or patch]' c ...
			      '  FaceColor       Set the FaceColor of patch arrows.  May be one of' c ...
			      '                  ''ymcrgbwkfn'' (f=flat,n=none) or a column-vector colorspec' c ...
			      '                  (or Nx3 matrix for N arrows). [1 1 1] (FaceC)' c ...
			      '  EdgeColor       Set the EdgeColor of patch arrows and Color of line' c ...
			      '                  arrows. [1 1 1] (EdgeC)' c ...
			      '  Color           Sets both FaceColor and EdgeColor properties.' c ...
			      '                  Included for compatibility with line objects.' c c ...
			      '  ARROW(''Start'',P,''Stop'',[]) or ARROW(''Start'',[],''Stop'',P), with size(P)=[N 3]' c ...
			      '  will create N-1 arrows, just like ARROW(''Start'',P1,''Stop'',P2), where' c c ...
			      '          / x1 y1 z1 \         /  x1   y1   z1  \         / x2 y2 z2 \' c ...
			      '          | .  .  .  |         |  .    .    .   |         | .  .  .  |' c ...
			      '      P = | .  .  .  |    P1 = |  .    .    .   |    P2 = | .  .  .  |' c ...
			      '          | .  .  .  |         |  .    .    .   |         | .  .  .  |' c ...
			      '          \ xN yN zN /         \ xN-1 yN-1 zN-1 /         \ xN yN zN /' c]);
		elseif (strcmp(arg1(1:4),'vers')),
			disp([c ...
			      'ARROW Version:  There are two compatibility problems for ARROW in Matlab' c ...
			      '                versions earlier than 4.2c:' c c c ...
			      '   1) In Matlab <4.2, the ''Tag'' property does not exist.  ARROW uses this' c ...
				  '      ''Tag'' to identify arrow objects for subsequent calls to ARROW.' c c ...
			      '      Solution: (a) delete all instances of the string' c c ...
			      '                         ,''Tag'',ArrowTag' c c ...
			      '                (b) and replace all instances of the string' c c ...
			      '                         strcmp(get(oh,''Tag''),ArrowTag)' c c ...
			      '                    with the string' c c ...
			      '                         all(size(ud)==[1 15])' c c c ...
			      '   2) In Matlab <4.2c, FINDOBJ is buggy (it can cause Matlab to crash if' c ...
				  '      more than 100 objects are returned).  ARROW uses FINDOBJ to make sure' c ...
				  '      that any handles it receives are truly handles to existing objects.' c c ...
			      '      Solution: replace the line' c c ...
			      '                         objs = findobj;' c c ...
			      '                with the lines' c c ...
			      '                         objs=0; k=1;' c ...
			      '                         while (k<=length(objs)),' c ...
			      '                             objs = [objs; get(objs(k),''Children'')];' c ...
			      '                             if (strcmp(get(objs(k),''Type''),''axes'')),' c ...
			      '                                 objs=[objs;get(objs(k),''XLabel''); ...' c ...
			      '                                            get(objs(k),''YLabel''); ...' c ...
			      '                                            get(objs(k),''ZLabel''); ...' c ...
			      '                                            get(objs(k),''Title'')];' c ...
			      '                             end;' c ...
			      '                             k=k+1;' c ...
			      '                         end;' c c]);
		elseif (strcmp(arg1(1:4),'demo')),
		 	% demo
			% create the data
			[x,y,z] = peaks;
			[ddd,iii]=max(z(:));
			axlim = [min(x(:)) max(x(:)) min(y(:)) max(y(:)) min(z(:)) max(z(:))];
			
			% modify it by inserting some NaN's
			[m,n] = size(z);
			m = floor(m/2);
			n = floor(n/2);
			z(1:m,1:n) = NaN*ones(m,n);
			
			% graph it
			clf('reset');
			hs=surf(x,y,z);
			xlabel('x'); ylabel('y');
			
			if (~strcmp(arg1(1:5),'demo2')),
				% set the view
				axis(axlim);
				zlabel('z');
				%shading('interp'); set(hs,'EdgeColor','k');
				view(viewmtx(-37.5,30,20));
				title('Demo of the capabilities of the ARROW function in 3-D');
				
				% Normal yellow arrow
				h1 = arrow([axlim(1) axlim(4) 4],[-.8 1.2 4], ...
				           'EdgeColor','y','FaceColor','y');
				
				% Normal white arrow, clipped by the surface
				h2 = arrow(axlim([1 4 6]),[0 2 4]);
				t=text(-2.4,2.7,7.7,'arrow clipped by surf');
				
				% Baseangle<90
				h3 = arrow([3 .125 3.5],[1.375 0.125 3.5],30,50);
				t2=text(3.1,.125,3.5,'local maximum');
				
				% Baseangle<90, fill and edge colors different
				h4 = arrow(axlim(1:2:5)*.5,[0 0 0],36,60,25, ...
				           'EdgeColor','b','FaceColor','c');
				t3=text(axlim(1)*.5,axlim(3)*.5,axlim(5)*.5-.75,'origin');
				set(t3,'HorizontalAlignment','center');
				
				% Baseangle>90, black fill
				h5 = arrow([-2.9 2.9 3],[-1.3 .4 3.2],30,120,[],6, ...
				           'EdgeColor','r','FaceColor','k','LineWidth',2);
				
				% Baseangle>90, no fill
				h6 = arrow([-2.9 2.9 1.3],[-1.3 .4 1.5],30,120,[],6, ...
				           'EdgeColor','r','FaceColor','none','LineWidth',2);
				
				% Stick arrow
				h7 = arrow([-1.6 -1.65 -6.5],[0 -1.65 -6.5],[],16,16);
				t4=text(-1.5,-1.65,-7.25,'global mininum');
				set(t4,'HorizontalAlignment','center');
				
				% Normal, black fill
				h8 = arrow([-1.4 0 -7.2],[-1.4 0 -3],'FaceColor','k');
				t5=text(-1.5,0,-7.75,'local minimum');
				set(t5,'HorizontalAlignment','center');
				
				% Gray fill, crossdir specified
				h9 = arrow([-3 2.2 -6],[-3 2.2 -.05],36,[],27,6,[],[0 -1 0], ...
				           'EdgeColor','w','FaceColor',.2*[1 1 1]);
				
				% a series of normal arrows, linearly spaced, crossdir specified
				h10y=(0:4)'/3;
				h10 = arrow([-3*ones(size(h10y)) h10y -6.5*ones(size(h10y))], ...
				            [-3*ones(size(h10y)) h10y -.05*ones(size(h10y))], ...
				            12,[],[],[],[],[0 -1 0]);
				
				% a series of normal arrows, linearly spaced
				h11x=(1:.33:2.8)';
				h11 = arrow([h11x -3*ones(size(h11x)) 6.5*ones(size(h11x))], ...
				            [h11x -3*ones(size(h11x)) -.05*ones(size(h11x))]);
				
				% series of black-filled arrows, radially oriented, crossdir specified
				h12x=2; h12y=-3; h12z=axlim(5)/2; h12xr=1; h12zr=h12z; ir=.15;or=.81;
				h12t=(0:11)'/6*pi;
				h12 = arrow([h12x+h12xr*cos(h12t)*ir h12y*ones(size(h12t))       ...
				            h12z+h12zr*sin(h12t)*ir],[h12x+h12xr*cos(h12t)*or    ...
				            h12y*ones(size(h12t)) h12z+h12zr*sin(h12t)*or],      ...
				            10,[],[],[],[],                                      ...
				            [-h12xr*sin(h12t) zeros(size(h12t)) h12zr*cos(h12t)],...
				            'FaceColor','none','EdgeColor','m');
				
				% series of normal arrows, tangentially oriented, crossdir specified
				or13=.91; h13t=(0:.5:12)'/6*pi;
				h13 = arrow([h12x+h12xr*cos(h13t)*or13 ...
				             h12y*ones(size(h13t))     ...
				             h12z+h12zr*sin(h13t)*or13],[],6);
				
				% arrow with no line ==> oriented upwards
				h14 = arrow([3 3 3],[3 3 3],30);
				t6=text(3,3,3.6,'no line'); set(t6,'HorizontalAlignment','center');
				
				% arrow with -- linestyle
				h15 = arrow([-.5 -3 -3],[1 -3 -3],'LineStyle','--','EdgeColor','g');
				
				if (nargout>=1), h=[h1;h2;h3;h4;h5;h6;h7;h8;h9;h10;h11;h12;h13;h14;h15]; end;
			else,	
				set(hs,'YData',10.^get(hs,'YData'));
				shading('interp');
				view(2);
				title('Demo of the capabilities of the ARROW function in 2-D');
				hold on; [C,H]=contour(x,y,z,20); hold off;
				for k=H', set(k,'ZData',(axlim(6)+1)*ones(size(get(k,'XData'))),...
				                'YData',10.^get(k,'YData'),'Color','k'); end;
				set(gca,'YScale','log');
				axis([axlim(1:2) 10.^axlim(3:4)]);
				
				% Normal yellow arrow
				h1 = arrow([axlim(1) 10^axlim(4) axlim(6)+2],[x(iii) 10^y(iii) axlim(6)+2], ...
				           'EdgeColor','y','FaceColor','y');
				
				% three arrows with varying fill, width, and baseangle
				h2 = arrow([-3 10^(-3) 10; -3 10^(-1.5) 10; -1.5 10^(-3) 10], ...
				           [-.03 10^(-.03) 10; -.03 10^(-1.5) 10; -1.5 10^(-.03) 10], ...
				           24,[90;60;120],[],[0;0;4]);
				set(h2(2),'EdgeColor','g','FaceColor','c');
				set(h2(3),'EdgeColor','m','FaceColor','r');
				if (nargout>=1), h=[h1;h2]; end;
			end;
		else,
			error(['ARROW got an unknown single-argument string ''' deblank(arg1) '''.']);
		end;
		return;
	end;
end;

% Check # of arguments
if (nargin==0), help arrow ; return;
elseif (nargout>3), error('ARROW produces at most 3 output arguments.');
end;

% find first property number
firstprop = nargin+1;
if (nargin<=3), % to speed things up a bit
	if (nargin==1),
	elseif (ischar(arg1)), firstprop=1;
	elseif (ischar(arg2)), firstprop=2;
	elseif (nargin==3),
		if (ischar(arg3)), firstprop=3; end;
	end;
else,
	for k=1:nargin,
		curarg = eval(['arg' num2str(k)]);
		if (ischar(curarg)),
			firstprop = k;
			break;
		end;
	end;
end;

% check property list
if (firstprop<=nargin),
	for k=firstprop:2:nargin,
		curarg = eval(['arg' num2str(k)]);
		if ((~ischar(curarg))|(min(size(curarg))~=1)),
			error('ARROW requires that a property name be a single string.');
		end;
	end;
	if (rem(nargin-firstprop,2)~=1),
		error(['ARROW requires that the property ''' eval(['arg' num2str(nargin)]) ...
		       ''' be paired with a property value.']);
	end;
end;

% default output
if (nargout>0), h=[]; end;
if (nargout>1), yy=[]; end;
if (nargout>2), zz=[]; end;

% set values to empty matrices
start      = [];
stop       = [];
len        = [];
baseangle  = [];
tipangle   = [];
wid        = [];
page       = [];
crossdir   = [];
ends       = [];
linewidth  = [];
linestyle  = [];
edgecolor  = [];
facecolor  = [];
ax         = [];
oldh       = [];
defstart      = [NaN NaN NaN];
defstop       = [NaN NaN NaN];
deflen        = 16;
defbaseangle  = 90;
deftipangle   = 16;
defwid        = 0;
defpage       = 0;
defcrossdir   = [NaN NaN NaN];
defends       = 1;
deflinewidth  = NaN;
deflinestyle  = '- ';
defedgecolor  = [1 1 1];
deffacecolor  = [1 1 1];
defoldh       = [];

% The 'Tag' we'll put on our arrows
ArrowTag = 'Arrow';

% check for oldstyle arguments
if (firstprop>2),
	% if old style arguments, get them
	if     (firstprop==3), start=arg1; stop=arg2;
	elseif (firstprop==4), start=arg1; stop=arg2; len=arg3(:);
	elseif (firstprop==5), start=arg1; stop=arg2; len=arg3(:); baseangle=arg4(:);
	elseif (firstprop==6), start=arg1; stop=arg2; len=arg3(:); baseangle=arg4(:); tipangle=arg5(:);
	elseif (firstprop==7), start=arg1; stop=arg2; len=arg3(:); baseangle=arg4(:); tipangle=arg5(:); wid=arg6(:);
	elseif (firstprop==8), start=arg1; stop=arg2; len=arg3(:); baseangle=arg4(:); tipangle=arg5(:); wid=arg6(:); page=arg7(:);
	elseif (firstprop==9), start=arg1; stop=arg2; len=arg3(:); baseangle=arg4(:); tipangle=arg5(:); wid=arg6(:); page=arg7(:); crossdir=arg8;
	else, error('ARROW takes at most 8 non-property arguments.');
	end;
elseif (firstprop==2),
	% assume arg1 is a set of handles
	oldh = arg1(:);
end;

% parse property pairs
extraprops=char([]);
for k=firstprop:2:nargin,
	prop = eval(['arg' num2str(k)]);
	val  = eval(['arg' num2str(k+1)]);
	prop = [lower(prop(:)') '      '];
	if     (all(prop(1:5)=='start')),      start      = val;
	elseif (all(prop(1:4)=='stop')),       stop       = val;
	elseif (all(prop(1:3)=='len')),        len        = val(:);
	elseif (all(prop(1:4)=='base')),       baseangle  = val(:);
	elseif (all(prop(1:3)=='tip')),        tipangle   = val(:);
	elseif (all(prop(1:3)=='wid')),        wid        = val(:);
	elseif (all(prop(1:4)=='page')),       page       = val;
	elseif (all(prop(1:5)=='cross')),      crossdir   = val;
	elseif (all(prop(1:4)=='norm')),       if (ischar(val)), crossdir=val; else, crossdir=val*sqrt(-1); end;
	elseif (all(prop(1:3)=='end')),        ends       = val;
	elseif (all(prop(1:5)=='linew')),      linewidth  = val(:);
	elseif (all(prop(1:5)=='lines')),      linestyle  = val;
	elseif (all(prop(1:5)=='color')),      edgecolor  = val; facecolor=val;
	elseif (all(prop(1:5)=='edgec')),      edgecolor  = val;
	elseif (all(prop(1:5)=='facec')),      facecolor  = val;
	elseif (all(prop(1:6)=='object')),     oldh       = val(:);
	elseif (all(prop(1:6)=='handle')),     oldh       = val(:);
	elseif (all(prop(1:5)=='userd')),      %ignore it
	else, extraprops=[extraprops ',arg' num2str(k) ',arg' num2str(k+1)];
	end;
end;

% Check if we got 'default' values
if (ischar(start     )),  s=lower([start(:)'      '   ']); if (all(s(1:3)=='def')), start      = defstart;      else, error(['ARROW does not recognize ''' start(:)'     ''' as a valid ''Start'' string.']); end; end;
if (ischar(stop      )),  s=lower([stop(:)'       '   ']); if (all(s(1:3)=='def')), stop       = defstop;       else, error(['ARROW does not recognize ''' stop(:)'      ''' as a valid ''Stop'' string.']); end; end;
if (ischar(len       )),  s=lower([len(:)'        '   ']); if (all(s(1:3)=='def')), len        = deflen;        else, error(['ARROW does not recognize ''' len(:)'       ''' as a valid ''Length'' string.']); end; end;
if (ischar(baseangle )),  s=lower([baseangle(:)'  '   ']); if (all(s(1:3)=='def')), baseangle  = defbaseangle;  else, error(['ARROW does not recognize ''' baseangle(:)' ''' as a valid ''BaseAngle'' string.']); end; end;
if (ischar(tipangle  )),  s=lower([tipangle(:)'   '   ']); if (all(s(1:3)=='def')), tipangle   = deftipangle;   else, error(['ARROW does not recognize ''' tipangle(:)'  ''' as a valid ''TipAngle'' string.']); end; end;
if (ischar(wid       )),  s=lower([wid(:)'        '   ']); if (all(s(1:3)=='def')), wid        = defwid;        else, error(['ARROW does not recognize ''' wid(:)'       ''' as a valid ''Width'' string.']); end; end;
if (ischar(crossdir  )),  s=lower([crossdir(:)'   '   ']); if (all(s(1:3)=='def')), crossdir   = defcrossdir;   else, error(['ARROW does not recognize ''' crossdir(:)'  ''' as a valid ''CrossDir'' or ''NormalDir'' string.']); end; end;
if (ischar(page      )),  s=lower([page(:)'       '   ']); if (all(s(1:3)=='def')), page       = defpage;       else, error(['ARROW does not recognize ''' page(:)'      ''' as a valid ''Page'' string.']); end; end;
if (ischar(ends      )),  s=lower([ends(:)'       '   ']); if (all(s(1:3)=='def')), ends       = defends;       end; end;
if (ischar(linewidth )),  s=lower([linewidth(:)'  '   ']); if (all(s(1:3)=='def')), linewidth  = deflinewidth;  else, error(['ARROW does not recognize ''' linewidth(:)' ''' as a valid ''LineWidth'' string.']); end; end;
if (ischar(linestyle )),  s=lower([linestyle(:)'  '   ']); if (all(s(1:3)=='def')), linestyle  = deflinestyle;  end; end;
if (ischar(edgecolor )),  s=lower([edgecolor(:)'  '   ']); if (all(s(1:3)=='def')), edgecolor  = defedgecolor;  end; end;
if (ischar(facecolor )),  s=lower([facecolor(:)'  '   ']); if (all(s(1:3)=='def')), facecolor  = deffacecolor;  end; end;
if (ischar(oldh      )),  s=lower([oldh(:)'       '   ']); if (all(s(1:3)=='def')), oldh       = [];            else, error(['ARROW does not recognize ''' oldh(:)'      ''' as a valid ''ObjectHandles'' string.']); end; end;

% check transpose on arguments; convert strings to numbers
if (((size(start   ,1)==2)|(size(start   ,1)==3))&((size(start   ,2)==1)|(size(start   ,2)>3))), start    = start';    end;
if (((size(stop    ,1)==2)|(size(stop    ,1)==3))&((size(stop    ,2)==1)|(size(stop    ,2)>3))), stop     = stop';     end;
if (((size(crossdir,1)==2)|(size(crossdir,1)==3))&((size(crossdir,2)==1)|(size(crossdir,2)>3))), crossdir = crossdir'; end;
if ((size(linestyle,2)>2)&(size(linestyle,1)<=2)), linestyle=linestyle'; end;
if (all(size(edgecolor))),
	if (ischar(edgecolor)),
		col = lower(edgecolor(:,1));
		edgecolor = zeros(length(col),3);
		ii=find(col=='y'); if (all(size(ii))), edgecolor(ii,:)=ones(length(ii),1)*[1 1 0]; end;
		ii=find(col=='m'); if (all(size(ii))), edgecolor(ii,:)=ones(length(ii),1)*[1 0 1]; end;
		ii=find(col=='c'); if (all(size(ii))), edgecolor(ii,:)=ones(length(ii),1)*[0 1 1]; end;
		ii=find(col=='r'); if (all(size(ii))), edgecolor(ii,:)=ones(length(ii),1)*[1 0 0]; end;
		ii=find(col=='g'); if (all(size(ii))), edgecolor(ii,:)=ones(length(ii),1)*[0 1 0]; end;
		ii=find(col=='b'); if (all(size(ii))), edgecolor(ii,:)=ones(length(ii),1)*[0 0 1]; end;
		ii=find(col=='w'); if (all(size(ii))), edgecolor(ii,:)=ones(length(ii),1)*[1 1 1]; end;
		ii=find(col=='k'); if (all(size(ii))), edgecolor(ii,:)=ones(length(ii),1)*[0 0 0]; end;
		ii=find(col=='f'); if (all(size(ii))), edgecolor(ii,:)=ones(length(ii),1)*[1 1 1]*inf; end;
		ii=find(col=='n'); if (all(size(ii))), edgecolor(ii,:)=ones(length(ii),1)*[1 1 1]*(-inf); end;
	elseif ((size(edgecolor,2)~=3)&(size(edgecolor,1)==3)),
		edgecolor=edgecolor';
	elseif (size(edgecolor,2)~=3),
		error('ARROW requires that color specifications must be a ?x3 RGB matrix.');
	end;
end;
if (all(size(facecolor))),
	if (ischar(facecolor)),
		col = lower(facecolor(:,1));
		facecolor = zeros(length(col),3);
		ii=find(col=='y'); if (all(size(ii))), facecolor(ii,:)=ones(length(ii),1)*[1 1 0]; end;
		ii=find(col=='m'); if (all(size(ii))), facecolor(ii,:)=ones(length(ii),1)*[1 0 1]; end;
		ii=find(col=='c'); if (all(size(ii))), facecolor(ii,:)=ones(length(ii),1)*[0 1 1]; end;
		ii=find(col=='r'); if (all(size(ii))), facecolor(ii,:)=ones(length(ii),1)*[1 0 0]; end;
		ii=find(col=='g'); if (all(size(ii))), facecolor(ii,:)=ones(length(ii),1)*[0 1 0]; end;
		ii=find(col=='b'); if (all(size(ii))), facecolor(ii,:)=ones(length(ii),1)*[0 0 1]; end;
		ii=find(col=='w'); if (all(size(ii))), facecolor(ii,:)=ones(length(ii),1)*[1 1 1]; end;
		ii=find(col=='k'); if (all(size(ii))), facecolor(ii,:)=ones(length(ii),1)*[0 0 0]; end;
		ii=find(col=='f'); if (all(size(ii))), facecolor(ii,:)=ones(length(ii),1)*[1 1 1]*inf; end;
		ii=find(col=='n'); if (all(size(ii))), facecolor(ii,:)=ones(length(ii),1)*[1 1 1]*(-inf); end;
	elseif ((size(facecolor,2)~=3)&(size(facecolor,1)==3)),
		facecolor=facecolor';
	elseif (size(facecolor,2)~=3),
		error('ARROW requires that color specifications must be a ?x3 RGB matrix.');
	end;
end;
if (all(size(ends))),
	if (ischar(ends)),
		endsorig = ends;
		col = lower([ends(:,1:min(3,size(ends,2))) ones(size(ends,1),max(0,3-size(ends,2)))*' ']);
		ends = NaN*ones(size(ends,1),1);
		oo = ones(1,size(ends,1));
		ii=find(all(col'==['non']'*oo)'); if (all(size(ii))), ends(ii)=ones(length(ii),1)*0; end;
		ii=find(all(col'==['sto']'*oo)'); if (all(size(ii))), ends(ii)=ones(length(ii),1)*1; end;
		ii=find(all(col'==['sta']'*oo)'); if (all(size(ii))), ends(ii)=ones(length(ii),1)*2; end;
		ii=find(all(col'==['bot']'*oo)'); if (all(size(ii))), ends(ii)=ones(length(ii),1)*3; end;
		if (any(isnan(ends))),
			ii = min(find(isnan(ends)));
			error(['ARROW does not recognize ''' deblank(endsorig(ii,:)) ''' as a valid ''Ends'' value.']);
		end;
	else,
		ends = ends(:);
	end;
end;
oldh = oldh(:);

% check object handles
if (all(size(oldh))),
	oldh = oldh.';
	objs = findobj;
	if (length(objs)==0), error('ARROW found no graphics handles.');
	elseif (length(objs)==1), objs=[objs;objs]; end;
	if (~all(any(objs(:,ones(1,length(oldh)))==oldh(ones(length(objs),1),:)))),
		error('ARROW got invalid object handles.');
	end;
	oldh = oldh.';
end;

% Check for an empty Start or Stop (but not both) with no object handles
if ((~all(size(oldh)))&(all(size(start))~=all(size(stop)))),
	if (~all(size(start))), start=stop; end;
	ii = find(all(diff(start)'==0)');
	if (size(start,1)==1),
		stop = start;
	elseif (length(ii)==size(start,1)-1)
		stop = start(1,:);
		start = stop;
	else,
		if (all(size(ii))),
			jj = (1:size(start,1))';
			jj(ii) = zeros(length(ii),1);
			jj = jj(find(jj>0));
			start = start(jj,:);
		end;
		stop = start(2:size(start,1),:);
		start = start(1:size(start,1)-1,:);
	end;
end;

% largest argument length
argsizes = [length(oldh) size(start,1) size(stop,1)                ...
            length(len) length(baseangle) length(tipangle)         ...
			length(wid) length(page) size(crossdir,1) length(ends) ...
            length(linewidth) size(edgecolor,1) size(facecolor,1)];
args=['length(ObjectHandle)  '; ...
      '#rows(Start)          '; ...
      '#rows(Stop)           '; ...
      'length(Length)        '; ...
      'length(BaseAngle)     '; ...
      'length(TipAngle)      '; ...
      'length(Width)         '; ...
      'length(Page)          '; ...
      '#rows(CrossDir)       '; ...
	  '#rows(Ends)           '; ...
      'length(LineWidth)     '; ...
      '#colors in EdgeColor  '; ...
      '#colors in FaceColor  '];
if (any(imag(crossdir(:))~=0)),
	args(9,:) = '#rows(NormalDir)      ';
end;
if (~all(size(oldh))),
	narrows = max(argsizes);
else,
	narrows = length(oldh);
end;

% Check size of arguments
ii = find((argsizes~=0)&(argsizes~=1)&(argsizes~=narrows));
if (all(size(ii))),
	s = args(ii',:);
	while ((size(s,2)>1)&((abs(s(:,size(s,2)))==0)|(abs(s(:,size(s,2)))==abs(' ')))),
		s = s(:,1:size(s,2)-1);
	end;
	s = [ones(length(ii),1)*'ARROW requires that  ' s ...
	     ones(length(ii),1)*['  equal the # of arrows (' num2str(narrows) ').' c]];
	s = s';
	s = s(:)';
	s = s(1:length(s)-1);
	error(char(s));
end;

% check element length in Start, Stop, and CrossDir
if (all(size(start))),
	if (size(start,2)==2),
		start = [start NaN*ones(size(start,1),1)];
	elseif (size(start,2)~=3),
		error('ARROW requires 2- or 3-element Start points.');
	end;
end;
if (all(size(stop))),
	if (size(stop,2)==2),
		stop = [stop NaN*ones(size(stop,1),1)];
	elseif (size(stop,2)~=3),
		error('ARROW requires 2- or 3-element Stop points.');
	end;
end;
if (all(size(crossdir))),
	if (size(crossdir,2)<3),
		crossdir = [crossdir NaN*ones(size(crossdir,1),3-size(crossdir,2))];
	elseif (size(crossdir,2)~=3),
		if (all(imag(crossdir(:))==0)),
			error('ARROW requires 2- or 3-element CrossDir vectors.');
		else,
			error('ARROW requires 2- or 3-element NormalDir vectors.');
		end;
	end;
end;

% fill empty arguments
if (~all(size(start     ))),   start      = [NaN NaN NaN];      end;
if (~all(size(stop      ))),   stop       = [NaN NaN NaN];      end;
if (~all(size(len       ))),   len        = NaN;                end;
if (~all(size(baseangle ))),   baseangle  = NaN;                end;
if (~all(size(tipangle  ))),   tipangle   = NaN;                end;
if (~all(size(wid       ))),   wid        = NaN;                end;
if (~all(size(page      ))),   page       = NaN;                end;
if (~all(size(crossdir  ))),   crossdir   = [NaN NaN NaN];      end;
if (~all(size(ends      ))),   ends       = NaN;                end;
if (~all(size(linewidth ))),   linewidth  = NaN;                end;
if (~all(size(linestyle ))),   linestyle  = char(['-']);      end; % was NaN
if (~all(size(edgecolor ))),   edgecolor  = [NaN NaN NaN];      end;
if (~all(size(facecolor ))),   facecolor  = [NaN NaN NaN];      end;

% expand single-column arguments
o = ones(narrows,1);
if (size(start     ,1)==1),   start      = o * start     ;   end;
if (size(stop      ,1)==1),   stop       = o * stop      ;   end;
if (length(len       )==1),   len        = o * len       ;   end;
if (length(baseangle )==1),   baseangle  = o * baseangle ;   end;
if (length(tipangle  )==1),   tipangle   = o * tipangle  ;   end;
if (length(wid       )==1),   wid        = o * wid       ;   end;
if (length(page      )==1),   page       = o * page      ;   end;
if (size(crossdir  ,1)==1),   crossdir   = o * crossdir  ;   end;
if (length(ends      )==1),   ends       = o * ends      ;   end;
if (length(linewidth )==1),   linewidth  = o * linewidth ;   end;
if (size(linestyle ,1)==1),   linestyle  = o * linestyle ;   end;
if (size(edgecolor ,1)==1),   edgecolor  = o * edgecolor ;   end;
if (size(facecolor ,1)==1),   facecolor  = o * facecolor ;   end;
ax = o * gca;
if (size(linestyle ,2)==1),   linestyle  = char([linestyle o*' ']); end;
linestyle = char(linestyle);

% if we've got handles, get the defaults from the handles
oldlinewidth = NaN*ones(narrows,1);
oldlinestyle = char(zeros(narrows,2));
oldedgecolor = NaN*ones(narrows,3);
oldfacecolor = NaN*ones(narrows,3);
if (all(size(oldh))),
	fromline = zeros(narrows,1);
	for k=1:narrows,
		oh = oldh(k);
		ud = get(oh,'UserData');
		isarrow = strcmp(get(oh,'Tag'),ArrowTag);
		ohtype = get(oh,'Type');
		ispatch = strcmp(ohtype,'patch');
		isline  = strcmp(ohtype,'line');
		if (isarrow|isline),
			% arrow UserData format: [start' stop' len base tip wid page crossdir' ends]
			if (isarrow),
				start0 = ud(1:3);
				stop0  = ud(4:6);
				if (isnan(len(k))),           len(k)        = ud( 7);   end;
				if (isnan(baseangle(k))),     baseangle(k)  = ud( 8);   end;
				if (isnan(tipangle(k))),      tipangle(k)   = ud( 9);   end;
				if (isnan(wid(k))),           wid(k)        = ud(10);   end;
				if (isnan(page(k))),          page(k)       = ud(11);   end;
				if (isnan(crossdir(k,1))),    crossdir(k,1) = ud(12);   end;
				if (isnan(crossdir(k,2))),    crossdir(k,2) = ud(13);   end;
				if (isnan(crossdir(k,3))),    crossdir(k,3) = ud(14);   end;
				if (isnan(ends(k))),          ends(k)       = ud(15);   end;
			end;
			if (isline),
				fromline(k) = 1;
				if (isarrow),
					fc = -inf*[1 1 1];
				else,
					fc = get(oh,'Color');
					x  = get(oh,'XData');
					y  = get(oh,'YData');
					z  = get(oh,'ZData');
					if (any(size(x)~=[1 2])|any(size(y)~=[1 2])),
						error('ARROW only converts two-point lines.');
					end;
					if (~all(size(z))), z=NaN*ones(size(x)); end;
					start0 = [x(1) y(1) z(1)];
					stop0  = [x(2) y(2) z(2)];
				end;
				ec = get(oh,'Color');
				ls = [get(oh,'LineStyle') '  ']; ls=char(ls(1:2));
				lw = get(oh,'LineWidth');
			else, % an arrow patch
				fc = get(oh,'FaceColor');if (ischar(fc)),
					if (strcmp(fc,'none')), fc=-inf*[1 1 1];
					elseif (strcmp(fc,'flat')), fc=inf*[1 1 1];
					else, fc=[1 1 1]; end;
					end;
				ec = get(oh,'EdgeColor');if (ischar(ec)),
					if (strcmp(ec,'none')), ec=-inf*[1 1 1];
					elseif (strcmp(ec,'flat')), ec=inf*[1 1 1];
					else, ec=[1 1 1]; end;
					end;
				ls = char('- ');
				lw = get(oh,'LineWidth');
			end;
			ax(k) = get(oh,'Parent');
		else,
			error(['ARROW cannot convert ' ohtype ' objects.']);
		end;
		oldlinewidth(k)   = lw;
		oldlinestyle(k,:) = ls;
		oldedgecolor(k,:) = ec;
		oldfacecolor(k,:) = fc;
		ii=find(isnan(start(k,:)));      if (all(size(ii))),  start(k,ii)=start0(ii);  end;
		ii=find(isnan(stop( k,:)));      if (all(size(ii))),  stop( k,ii)=stop0( ii);  end;
		if (isnan(linewidth(k))),        linewidth(k)     = lw;  end;
		if (isnan(linestyle(k,1))),      linestyle(k,1:2) = ls;  end;
		if (any(isnan(facecolor(k,:)))), facecolor(k,:)   = fc;  end;
		if (any(isnan(edgecolor(k,:)))), edgecolor(k,:)   = ec;  end;
	end;
else
	fromline = [];
end;

% set up the UserData data
% (do it here so it is not corrupted by log10's and such)
ud = [start stop len baseangle tipangle wid page crossdir ends];

% Set Page defaults
if isnan(page)
    page = ~isnan(page);
end;

% Get axes limits, range, min; correct for aspect ratio and log scale
axm      = zeros(3,narrows);
axr      = zeros(3,narrows);
ap       = zeros(2,narrows);
xyzlog   = zeros(3,narrows);
limmin   = zeros(2,narrows);
limrange = zeros(2,narrows);
oneax = all(ax==ax(1));
if (oneax),
	T = zeros(4,4);
	invT = zeros(4,4);
else,
	T        = zeros(16,narrows);
	invT     = zeros(16,narrows);
end;
axnotdone = ones(size(ax));
while (any(axnotdone)),
	ii = min(find(axnotdone));
	curax = ax(ii);
	curpage = page(ii);
	% get axes limits and aspect ratio
	axl = [get(curax,'XLim'); get(curax,'YLim'); get(curax,'ZLim')];
	ar = get(curax,'DataAspectRatio');
	% get axes size in pixels (points)
	u = get(curax,'Units');
	axposoldunits = get(curax,'Position');
	if (curpage),
		curfig = get(curax,'Parent');
		pu = get(curfig,'PaperUnits');
		set(curfig,'PaperUnits','points');
		pp = get(curfig,'PaperPosition');
		set(curfig,'PaperUnits',pu);
		set(curax,'Units','normalized');
		curap = get(curax,'Position');
		curap = pp.*curap;
	else,
		set(curax,'Units','pixels');
		curap = get(curax,'Position');
	end;
	set(curax,'Units',u);
	set(curax,'Position',axposoldunits);
	% adjust limits for log scale on axes
	curxyzlog = [strcmp(get(curax,'XScale'),'log'); ...
	             strcmp(get(curax,'YScale'),'log'); ...
	             strcmp(get(curax,'ZScale'),'log')];
	if (any(curxyzlog)),
		ii = find([curxyzlog;curxyzlog]);
		if (any(axl(ii)<=0)),
			error('ARROW does not support non-positive limits on log-scaled axes.');
		else,
			axl(ii) = log10(axl(ii));
		end;
	end;
	% correct for aspect ratio
	if (~isnan(ar(1))),
		if (curap(3) < ar(1)*curap(4)),
			curap(2) = curap(2) + (curap(4)-curap(3)/ar(1))/2;
			curap(4) = curap(3)/ar(1);
		else,
			curap(1) = curap(1) + (curap(3)-curap(4)*ar(1))/2;
			curap(3) = curap(4)*ar(1);
		end;
	end;
	% correct for 'equal'
	% may only want to do this for 2-D views, but seems right for 3-D also
	if (~isnan(ar(2))),
		if ((curap(3)/(axl(1,2)-axl(1,1)))/(curap(4)/(axl(2,2)-axl(2,1)))>ar(2)),
			incr = curap(3)*(axl(2,2)-axl(2,1))/(curap(4)*ar(2)) - (axl(1,2)-axl(1,1));
			axl(1,:) = axl(1,:) + incr/2*[-1 1];
		else,
			incr = ar(2)*(axl(1,2)-axl(1,1))*curap(4)/curap(3) - (axl(2,2)-axl(2,1));
			axl(2,:) = axl(2,:) + incr/2*[-1 1];
		end;
	end;
	% compute the range of 2-D values
	curT = get(curax,'Xform');
	lim = curT*[0 1 0 1 0 1 0 1;0 0 1 1 0 0 1 1;0 0 0 0 1 1 1 1;1 1 1 1 1 1 1 1];
	lim = lim(1:2,:)./([1;1]*lim(4,:));
	curlimmin = min(lim')';
	curlimrange = max(lim')' - curlimmin;
	curinvT = inv(curT);
	if (~oneax),
		curT = curT.';
		curinvT = curinvT.';
		curT = curT(:);
		curinvT = curinvT(:);
	end;
	% check which arrows to which cur corresponds
	ii = find((ax==curax)&(page==curpage));
	oo = ones(1,length(ii));
	axr(:,ii)      = diff(axl')' * oo;
	axm(:,ii)      = axl(:,1)    * oo;
	ap(:,ii)       = curap(3:4)' * oo;
	xyzlog(:,ii)   = curxyzlog   * oo;
	limmin(:,ii)   = curlimmin   * oo;
	limrange(:,ii) = curlimrange * oo;
	if (oneax),
		T    = curT;
		invT = curinvT;
	else,
		T(:,ii)    = curT    * oo;
		invT(:,ii) = curinvT * oo;
	end;
	axnotdone(ii) = zeros(1,length(ii));
end;

% correct for log scales
curxyzlog = xyzlog.';
ii = find(curxyzlog(:));
if (all(size(ii))),
	start(   ii) = real(log10(start(   ii)));
	stop(    ii) = real(log10(stop(    ii)));
	if (all(imag(crossdir(ii))==0)),
		crossdir(ii) = real(log10(crossdir(ii)));
	else,
		jj = find(imag(crossdir(ii))==0);
		if (all(size(jj))), crossdir(jj) = real(log10(crossdir(jj))); end;
		jj = find(imag(crossdir(ii))~=0);
		if (all(size(jj))), crossdir(jj) = real(log10(imag(crossdir(jj))))*sqrt(-1); end;
	end;
end;

% take care of defaults, page was done above
ii=find(isnan(start(:)       ));  if (all(size(ii))),  start(ii)       = axm(ii)+axr(ii)/2;                end;
ii=find(isnan(stop(:)        ));  if (all(size(ii))),  stop(ii)        = axm(ii)+axr(ii)/2;                end;
ii=find(isnan(crossdir(:)    ));  if (all(size(ii))),  crossdir(ii)    = zeros(length(ii),1);              end;
ii=find(isnan(len            ));  if (all(size(ii))),  len(ii)         = ones(length(ii),1)*deflen;        end;
ii=find(isnan(baseangle      ));  if (all(size(ii))),  baseangle(ii)   = ones(length(ii),1)*defbaseangle;  end;
ii=find(isnan(tipangle       ));  if (all(size(ii))),  tipangle(ii)    = ones(length(ii),1)*deftipangle;   end;
ii=find(isnan(wid            ));  if (all(size(ii))),  wid(ii)         = ones(length(ii),1)*defwid;        end;
ii=find(isnan(ends           ));  if (all(size(ii))),  ends(ii)        = ones(length(ii),1)*defends;       end;
ii=find(isnan(linewidth      ));  if (all(size(ii))),  linewidth(ii)   = ones(length(ii),1)*deflinewidth;  end;
ii=find(any(isnan(edgecolor')));  if (all(size(ii))),  edgecolor(ii,:) = ones(length(ii),1)*defedgecolor;  end;
ii=find(any(isnan(facecolor')));  if (all(size(ii))),  facecolor(ii,:) = ones(length(ii),1)*deffacecolor;  end;
ii=find(isnan(linestyle(:,1) ));  if (all(size(ii))),  linestyle(ii,:) = ones(length(ii),1)*deflinestyle;  end;
ii=find(isnan(linestyle(:,2) ));  if (all(size(ii))),  linestyle(ii,2) = ones(length(ii),1)*' ';           end; %just in case

% transpose all values
start     = start.';
stop      = stop.';
len       = len.';
baseangle = baseangle.';
tipangle  = tipangle.';
wid       = wid.';
page      = page.';
crossdir  = crossdir.';
ends      = ends.';
linewidth = linewidth.';
linestyle = linestyle.';
facecolor = facecolor.';
edgecolor = edgecolor.';
fromline  = fromline.';
ax        = ax.';
oldlinewidth = oldlinewidth.';
oldlinestyle = oldlinestyle.';
oldedgecolor = oldedgecolor.';
oldfacecolor = oldfacecolor.';

% given x, a 3xN matrix of points in 3-space;
% want to convert to X, the corresponding 4xN 2-space matrix
%
%   tmp1=[(x-axm)./axr; ones(1,size(x,1))];
%   if (oneax), X=T*tmp1;
%   else, tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=T.*tmp1;
%         tmp2=zeros(4,4*N); tmp2(:)=tmp1(:);
%         X=zeros(4,N); X(:)=sum(tmp2)'; end;
%   X = X ./ (X(:,4)*ones(1,4));

% for all points with start==stop, start=stop-(verysmallvalue)*(up-direction);
ii = find(all(start==stop));
if (all(size(ii))),
	% find an arrowdir vertical on screen and perpendicular to viewer
	%	transform to 2-D
		tmp1 = [(stop(:,ii)-axm(:,ii))./axr(:,ii);ones(1,length(ii))];
		if (oneax), twoD=T*tmp1;
		else, tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=T(:,ii).*tmp1;
		      tmp2=zeros(4,4*length(ii)); tmp2(:)=tmp1(:);
		      twoD=zeros(4,length(ii)); twoD(:)=sum(tmp2)'; end;
		twoD=twoD./(ones(4,1)*twoD(4,:));
	%	move the start point down just slightly
		tmp1 = twoD + [0;-1/1000;0;0]*(limrange(2,ii)./ap(2,ii));
	%	transform back to 3-D
		if (oneax), threeD=invT*tmp1;
		else, tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=invT(:,ii).*tmp1;
		      tmp2=zeros(4,4*length(ii)); tmp2(:)=tmp1(:);
		      threeD=zeros(4,length(ii)); threeD(:)=sum(tmp2)'; end;
		start(:,ii) = (threeD(1:3,:)./(ones(3,1)*threeD(4,:))).*axr(:,ii)+axm(:,ii);
end;

% compute along-arrow points
%	transform Start points
	tmp1=[(start-axm)./axr;ones(1,narrows)];
	if (oneax), X0=T*tmp1;
	else, tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=T.*tmp1;
	      tmp2=zeros(4,4*narrows); tmp2(:)=tmp1(:);
	      X0=zeros(4,narrows); X0(:)=sum(tmp2)'; end;
	X0=X0./(ones(4,1)*X0(4,:));
%	transform Stop points
	tmp1=[(stop-axm)./axr;ones(1,narrows)];
	if (oneax), Xf=T*tmp1;
	else, tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=T.*tmp1;
	      tmp2=zeros(4,4*narrows); tmp2(:)=tmp1(:);
	      Xf=zeros(4,narrows); Xf(:)=sum(tmp2)'; end;
	Xf=Xf./(ones(4,1)*Xf(4,:));
%	compute pixel distance between points
	D = sqrt(sum(((Xf(1:2,:)-X0(1:2,:)).*(ap./limrange)).^2));
%	compute and modify along-arrow distances
	len1 = len;
	len2 = len - (len.*tan(tipangle/180*pi)-wid/2).*tan((90-baseangle)/180*pi);
	slen0 = zeros(1,narrows);
	slen1 = len1 .* ((ends==2)|(ends==3));
	slen2 = len2 .* ((ends==2)|(ends==3));
	len0 = zeros(1,narrows);
	len1  = len1 .* ((ends==1)|(ends==3));
	len2  = len2 .* ((ends==1)|(ends==3));
	%	for no start arrowhead
		ii=find((ends==1)&(D<len2));
		if (all(size(ii))),
			slen0(ii) = D(ii)-len2(ii);
		end;
	%	for no end arrowhead
		ii=find((ends==2)&(D<slen2));
		if (all(size(ii))),
			len0(ii) = D(ii)-slen2(ii);
		end;
	len1  = len1  + len0;
	len2  = len2  + len0;
	slen1 = slen1 + slen0;
	slen2 = slen2 + slen0;
%	compute stoppoints
	tmp1=X0.*(ones(4,1)*(len0./D))+Xf.*(ones(4,1)*(1-len0./D));
	if (oneax), tmp3=invT*tmp1;
	else, tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=invT.*tmp1;
	      tmp2=zeros(4,4*narrows); tmp2(:)=tmp1(:);
	      tmp3=zeros(4,narrows); tmp3(:)=sum(tmp2)'; end;
	stoppoint = tmp3(1:3,:)./(ones(3,1)*tmp3(4,:)).*axr+axm;
%	compute tippoints
	tmp1=X0.*(ones(4,1)*(len1./D))+Xf.*(ones(4,1)*(1-len1./D));
	if (oneax), tmp3=invT*tmp1;
	else, tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=invT.*tmp1;
	      tmp2=zeros(4,4*narrows); tmp2(:)=tmp1(:);
	      tmp3=zeros(4,narrows); tmp3(:)=sum(tmp2)'; end;
	tippoint = tmp3(1:3,:)./(ones(3,1)*tmp3(4,:)).*axr+axm;
%	compute basepoints
	tmp1=X0.*(ones(4,1)*(len2./D))+Xf.*(ones(4,1)*(1-len2./D));
	if (oneax), tmp3=invT*tmp1;
	else, tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=invT.*tmp1;
	      tmp2=zeros(4,4*narrows); tmp2(:)=tmp1(:);
	      tmp3=zeros(4,narrows); tmp3(:)=sum(tmp2)'; end;
	basepoint = tmp3(1:3,:)./(ones(3,1)*tmp3(4,:)).*axr+axm;
%	compute startpoints
	tmp1=X0.*(ones(4,1)*(1-slen0./D))+Xf.*(ones(4,1)*(slen0./D));
	if (oneax), tmp3=invT*tmp1;
	else, tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=invT.*tmp1;
	      tmp2=zeros(4,4*narrows); tmp2(:)=tmp1(:);
	      tmp3=zeros(4,narrows); tmp3(:)=sum(tmp2)'; end;
	startpoint = tmp3(1:3,:)./(ones(3,1)*tmp3(4,:)).*axr+axm;
%	compute stippoints
	tmp1=X0.*(ones(4,1)*(1-slen1./D))+Xf.*(ones(4,1)*(slen1./D));
	if (oneax), tmp3=invT*tmp1;
	else, tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=invT.*tmp1;
	      tmp2=zeros(4,4*narrows); tmp2(:)=tmp1(:);
	      tmp3=zeros(4,narrows); tmp3(:)=sum(tmp2)'; end;
	stippoint = tmp3(1:3,:)./(ones(3,1)*tmp3(4,:)).*axr+axm;
%	compute sbasepoints
	tmp1=X0.*(ones(4,1)*(1-slen2./D))+Xf.*(ones(4,1)*(slen2./D));
	if (oneax), tmp3=invT*tmp1;
	else, tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=invT.*tmp1;
	      tmp2=zeros(4,4*narrows); tmp2(:)=tmp1(:);
	      tmp3=zeros(4,narrows); tmp3(:)=sum(tmp2)'; end;
	sbasepoint = tmp3(1:3,:)./(ones(3,1)*tmp3(4,:)).*axr+axm;

% compute cross-arrow directions for arrows with NormalDir specified
if (any(imag(crossdir(:))~=0)),
	ii = find(any(imag(crossdir)~=0));
	crossdir(:,ii) = cross((stop(:,ii)-start(:,ii))./axr(:,ii), ...
	                       imag(crossdir(:,ii))).*axr(:,ii);
end;

% compute cross-arrow directions
basecross  = crossdir + basepoint;
tipcross   = crossdir + tippoint;
sbasecross = crossdir + sbasepoint;
stipcross  = crossdir + stippoint;
ii = find(all(crossdir==0)|any(isnan(crossdir)));
if (all(size(ii))),
	numii = length(ii);
	%	transform start points
		tmp1 = [basepoint(:,ii) tippoint(:,ii) sbasepoint(:,ii) stippoint(:,ii)];
		tmp1 = (tmp1-axm(:,[ii ii ii ii])) ./ axr(:,[ii ii ii ii]);
		tmp1 = [tmp1; ones(1,4*numii)];
		if (oneax), X0=T*tmp1;
		else, tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=T(:,[ii ii ii ii]).*tmp1;
		      tmp2=zeros(4,16*numii); tmp2(:)=tmp1(:);
		      X0=zeros(4,4*numii); X0(:)=sum(tmp2)'; end;
		X0=X0./(ones(4,1)*X0(4,:));
	%	transform stop points
		tmp1 = [(2*stop(:,ii)-start(:,ii)-axm(:,ii))./axr(:,ii);ones(1,numii)];
		tmp1 = [tmp1 tmp1 tmp1 tmp1];
		if (oneax), Xf=T*tmp1;
		else, tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=T(:,[ii ii ii ii]).*tmp1;
		      tmp2=zeros(4,16*numii); tmp2(:)=tmp1(:);
		      Xf=zeros(4,4*numii); Xf(:)=sum(tmp2)'; end;
		Xf=Xf./(ones(4,1)*Xf(4,:));
	%	compute perpendicular directions
		pixfact = ((limrange(1,ii)./limrange(2,ii)).*(ap(2,ii)./ap(1,ii))).^2;
		pixfact = [pixfact pixfact pixfact pixfact];
		pixfact = [pixfact;1./pixfact];
		[dummyval,jj] = max(abs(Xf(1:2,:)-X0(1:2,:)));
		jj1 = ((1:4)'*ones(1,length(jj))==ones(4,1)*jj);
		jj2 = ((1:4)'*ones(1,length(jj))==ones(4,1)*(3-jj));
		jj3 = jj1(1:2,:);
		Xp = X0;
		Xp(jj2) = X0(jj2) + ones(sum(jj2(:)),1);
		Xp(jj1) = X0(jj1) - (Xf(jj2)-X0(jj2))./(Xf(jj1)-X0(jj1)) .* pixfact(jj3);
	%	inverse transform the cross points
		if (oneax), Xp=invT*Xp;
		else, tmp1=[Xp;Xp;Xp;Xp]; tmp1=invT(:,[ii ii ii ii]).*tmp1;
		      tmp2=zeros(4,16*numii); tmp2(:)=tmp1(:);
		      Xp=zeros(4,4*numii); Xp(:)=sum(tmp2)'; end;
		Xp=(Xp(1:3,:)./(ones(3,1)*Xp(4,:))).*axr(:,[ii ii ii ii])+axm(:,[ii ii ii ii]);
		basecross(:,ii)  = Xp(:,0*numii+(1:numii));
		tipcross(:,ii)   = Xp(:,1*numii+(1:numii));
		sbasecross(:,ii) = Xp(:,2*numii+(1:numii));
		stipcross(:,ii)  = Xp(:,3*numii+(1:numii));
end;

% compute all points
%	compute start points
	axm11 = [axm axm axm axm axm axm axm axm axm axm axm];
	axr11 = [axr axr axr axr axr axr axr axr axr axr axr];
	st = [stoppoint tippoint basepoint sbasepoint stippoint startpoint stippoint sbasepoint basepoint tippoint stoppoint];
	tmp1 = (st - axm11) ./ axr11;
	tmp1 = [tmp1; ones(1,size(tmp1,2))];
	if (oneax), X0=T*tmp1;
	else, tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=[T T T T T T T T T T T].*tmp1;
	      tmp2=zeros(4,44*narrows); tmp2(:)=tmp1(:);
	      X0=zeros(4,11*narrows); X0(:)=sum(tmp2)'; end;
	X0=X0./(ones(4,1)*X0(4,:));
%	compute stop points
	tmp1 = ([start tipcross basecross sbasecross stipcross stop stipcross sbasecross basecross tipcross start] ...
	     - axm11) ./ axr11;
	tmp1 = [tmp1; ones(1,size(tmp1,2))];
	if (oneax), Xf=T*tmp1;
	else, tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=[T T T T T T T T T T T].*tmp1;
	      tmp2=zeros(4,44*narrows); tmp2(:)=tmp1(:);
	      Xf=zeros(4,11*narrows); Xf(:)=sum(tmp2)'; end;
	Xf=Xf./(ones(4,1)*Xf(4,:));
%	compute lengths
	len0  = len.*((ends==1)|(ends==3)).*tan(tipangle/180*pi);
	slen0 = len.*((ends==2)|(ends==3)).*tan(tipangle/180*pi);
	le = [zeros(1,narrows) len0 wid/2 wid/2 slen0 zeros(1,narrows) -slen0 -wid/2 -wid/2 -len0 zeros(1,narrows)];
	aprange = ap./limrange;
	aprange = [aprange aprange aprange aprange aprange aprange aprange aprange aprange aprange aprange];
	D = sqrt(sum(((Xf(1:2,:)-X0(1:2,:)).*aprange).^2));
	Dii=find(D==0); if (all(size(Dii))), D=D+(D==0); le(Dii)=zeros(1,length(Dii)); end; %should fix DivideByZero warnings
	tmp1 = X0.*(ones(4,1)*(1-le./D)) + Xf.*(ones(4,1)*(le./D));
%	inverse transform
	if (oneax), tmp3=invT*tmp1;
	else, tmp1=[tmp1;tmp1;tmp1;tmp1]; tmp1=[invT invT invT invT invT invT invT invT invT invT invT].*tmp1;
	      tmp2=zeros(4,44*narrows); tmp2(:)=tmp1(:);
	      tmp3=zeros(4,11*narrows); tmp3(:)=sum(tmp2)'; end;
	pts = tmp3(1:3,:)./(ones(3,1)*tmp3(4,:)) .* axr11 + axm11;

% correct for ones where the crossdir was specified
ii = find(~(all(crossdir==0)|any(isnan(crossdir))));
if (all(size(ii))),
	D1 = [pts(:,1*narrows+ii)-pts(:,9*narrows+ii) ...
	      pts(:,2*narrows+ii)-pts(:,8*narrows+ii) ...
	      pts(:,3*narrows+ii)-pts(:,7*narrows+ii) ...
	      pts(:,4*narrows+ii)-pts(:,6*narrows+ii) ...
	      pts(:,6*narrows+ii)-pts(:,4*narrows+ii) ...
	      pts(:,7*narrows+ii)-pts(:,3*narrows+ii) ...
	      pts(:,8*narrows+ii)-pts(:,2*narrows+ii) ...
	      pts(:,9*narrows+ii)-pts(:,1*narrows+ii)]/2;
	ii = ii'*ones(1,8) + ones(length(ii),1)*[1:4 6:9]*narrows;
	ii = ii(:)';
	pts(:,ii) = st(:,ii) + D1;
end;


% readjust for log scale on axes
tmp1=[xyzlog xyzlog xyzlog xyzlog xyzlog xyzlog xyzlog xyzlog xyzlog xyzlog xyzlog];
ii = find(tmp1(:)); if (all(size(ii))), pts(ii)=10.^pts(ii); end;

% compute the x,y,z coordinates of the patches;
ii = narrows*(0:10)'*ones(1,narrows) + ones(11,1)*(1:narrows);
ii = ii(:)';
x = zeros(11,narrows);
y = zeros(11,narrows);
z = zeros(11,narrows);
x(:) = pts(1,ii)';
y(:) = pts(2,ii)';
z(:) = pts(3,ii)';

% do the output
if (nargout<=1),
%	% create or modify the patches
	p = (linestyle(1,:)=='-')&(linestyle(2,:)==' ');
	if (all(size(oldh))),
		H = oldh;
	else,
		fromline = p;
		H = zeros(narrows,1);
	end;
%	% new patches
	ii = find(p&fromline);
	if (all(size(ii))),
		if (all(size(oldh))), delete(oldh(ii)); end;
		H(ii) = patch(x(:,ii),y(:,ii),z(:,ii),zeros(11,length(ii)),'Tag',ArrowTag);
	end;
%	% new lines
	ii = find((~p)&(~fromline));
	if (all(size(ii))),
		if (all(size(oldh))), delete(oldh(ii)); end;
		H(ii) = line(x(:,ii),y(:,ii),z(:,ii),'Tag',ArrowTag);
	end;
%	% additional properties
	for k=1:narrows,
		s = 'set(H(k),''UserData'',ud(k,:)';
		if (p(k)~=fromline(k)), % modifying the data
			s = [s ',''XData'',x(:,k)'',''YData'',y(:,k)'',''ZData'',z(:,k)'''];
		end;
		if ((~isnan(linewidth(k)))&(linewidth(k)~=oldlinewidth(k))),
			s=[s ',''LineWidth'',linewidth(k)'];
		end;
		if (p(k)), % a patch
			if (any(edgecolor(:,k)~=oldedgecolor(:,k))|(fromline(k))),
				if (edgecolor(1,k)==inf),
					s=[s ',''EdgeColor'',''flat'''];
				elseif (edgecolor(1,k)==-inf),
					s=[s ',''EdgeColor'',''none'''];
				else,
					s=[s ',''EdgeColor'',edgecolor(:,k)'''];
				end;
			end;
			if (any(facecolor(:,k)~=oldfacecolor(:,k))|(fromline(k))),
				if (facecolor(1,k)==inf),
					s=[s ',''FaceColor'',''flat'''];
				elseif (facecolor(1,k)==-inf),
					s=[s ',''FaceColor'',''none'''];
				else,
					s=[s ',''FaceColor'',facecolor(:,k)'''];
				end;
			end;
		else, % a line
			s = [s ',''LineStyle'',''' deblank(char(linestyle(:,k)')) ''''];
			if (any(facecolor(:,k)~=oldfacecolor(:,k))|any(edgecolor(:,k)~=oldedgecolor(:,k))), % a line
				if (isinf(edgecolor(1,k))&(facecolor(1,k)~=-inf)),
					s = [s ',''Color'',facecolor(:,k)'''];
				else,
					s = [s ',''Color'',edgecolor(:,k)'''];
				end;
			end;
		end;
		eval([s extraprops ');']);
	end;
	% set the output
	if (nargout>0), h=H; end;
else,
	% don't create the patch, just return the data
	h=x;
	yy=y;
	zz=z;
end;
