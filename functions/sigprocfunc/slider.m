% slider() - add slider to a figure
%
% Usage:
%   >>  slider( handler, horiz, vert, horizmag, vertmag);
%
% Inputs:
%   handler  - figure handler (for the current figure, use gcf)
%   horiz    - [0|1] add a horizontal slider
%   vert     - [0|1] add a horizontal slider
%   horizmag - magnify the width of the figure before adding the slider.
%              Default is 1.
%   vertmag  - magnify the height of the figure before adding the slider.
%              Default is 1.
%   allowsup - [0|1] allow suppression of slider by the 'x' button.
%              Default is 1.
% 
% Note:
%   clicking on the 'x' the right corner restores the original setting   
%
% Example: figure; plot(1:10); slider(gcf, 1, 1, 2, 2);
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function slider( handler, horiz, vert, horizmag, vertmag, allowsup);

if nargin < 2 
	help slider;
	return;
end
if nargin < 3
	vert = 0;
end
if nargin < 4
	horizmag = 1;
end
if nargin < 5
	vertmag = 1;
end
if nargin < 6
	allowsup = 1;
end

pos   = get(gcf, 'position');
width  = 5/pos(3)*400;
height = 5/pos(4)*400;

pos = get(gca,'position'); % plot relative to current axes
q = [0.13 0.11 0 0];
s = [0.0077 0.0081 0.0077 0.0081];
if vert
	h = uicontrol('Parent',handler, ...
 	'style', 'slider', 'Units','Normalized', 'userdata', 1, 'value', 1, ...
	'Position', [113-width -9 width 119].*s+q, 'tag', 'winslider', ...
	'string','vertslider', 'callback', ...
	['if isempty(gcbf), fig = gcf; else fig = gcbf; end;' ...
	 'h = findobj(''parent'', fig);' ...
	 'h2 = findobj(''parent'', fig, ''tag'', ''winslider'');' ...
	 'h = setdiff_bc( h, h2);' ...
	 'curobj = findobj(''parent'', fig, ''string'', ''vertslider'');' ...
	 'shift = get(curobj, ''userdata'') - get(curobj, ''value'');' ...
	 'set( curobj, ''userdata'', get(curobj, ''value''));' ... 
	 'for i = 1:length(h),' ...
	 '   curpos = get( h(i), ''position'');' ...
	 '   set( h(i), ''position'', [ curpos(1) curpos(2)+' num2str(vertmag-1) '*shift curpos(3:end)]);' ...
	 'end;' ...
	 'clear h2 h shift i curpos fig curobj;'] );
end

if horiz
	hz = uicontrol('Parent',handler, ...
 	'style', 'slider', 'Units','Normalized', 'userdata', 1, 'value', 1, ...
	'Position', [-17 -19+height 125 height].*s+q, 'tag', 'winslider', ...
	'string','horizslider', 'callback', ...
	['if isempty(gcbf), fig = gcf; else fig = gcbf; end;' ...
	 'h = findobj(''parent'', fig);' ...
	 'h2 = findobj(''parent'', fig, ''tag'', ''winslider'');' ...
	 'h = setdiff_bc( h, h2);' ...
	 'curobj = findobj(''parent'', fig, ''string'', ''horizslider'');' ...
	 'shift = get(curobj, ''userdata'') - get(curobj, ''value'');' ...
	 'set( curobj, ''userdata'', get(curobj, ''value''));' ... 
	 'for i = 1:length(h),' ...
	 '   curpos = get( h(i), ''position'');' ...
	 '   set( h(i), ''position'', [ curpos(1)+' num2str(horizmag-1) '*shift curpos(2:end)]);' ...
	 'end;' ...
	 'clear h2 h shift i curpos fig curobj;'] );
end

% button to remove the slider
% ---------------------------
but = uicontrol( 'style', 'pushbutton', 'Units','Normalized', ...
		   'string', 'x', 'position', [113-width -19+height width height].*s+q, 'tag', 'winslider', 'callback', ...
		  ['hx = findobj(''parent'', gcbf, ''tag'', ''winslider'');' ... % put slider to their extremities
		   'hx = setdiff_bc(hx, gcbo);' ...
		   'set(hx, ''value'', 1);' ...
		   'eval(get(hx(1), ''callback''));' ...
		   'if length(hx) >1, eval(get(hx(2), ''callback'')); end;' ...
		   'h = findobj(''parent'', gcbf);' ... % recompute positions
		   'for i = 1:length(h),' ...
		   '   curpos = get( h(i), ''position'');' ...
		   '   set( h(i), ''position'', [(curpos(1)+(' num2str(horizmag) '-1))/' num2str(horizmag) '  (curpos(2)+(' num2str(vertmag) '-1))/' num2str(vertmag) ' curpos(3)/' num2str(horizmag) ' curpos(4)/' num2str(vertmag) ']);' ...
		   'end;' ...
		   'clear h hx curpos;' ...
		   'delete( findobj(''parent'', gcbf, ''tag'', ''winslider'') );' ]);

if ~allowsup
	set(but, 'enable', 'off');
end
	
% magnify object in the window
% ----------------------------
h = findobj('parent', handler);
set( h, 'units', 'normalized');
h2 = findobj('parent', handler, 'tag', 'winslider');
h = setdiff_bc( h, h2);
for i = 1:length(h)
	curpos = get( h(i), 'position');
	set( h(i), 'position', [curpos(1)*horizmag-(horizmag-1)  curpos(2)*vertmag-(vertmag-1) curpos(3)*horizmag curpos(4)*vertmag]);
end

if horiz
	% set the horizontal axis to 0	
	% ----------------------------
	set(hz, 'value', 0);
	eval(get(hz, 'callback'));
end
return;

