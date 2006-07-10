% pop_comments() - edit comments
%
% Usage:
%   >> newcomments = pop_comments( oldcomments);
%   >> newcomments = pop_comments( oldcomments, title, newcomments, concat);
%
% Inputs:
%   oldcomments - old comments (string or cell array of strings)
%   title       - optional window title (string)
%   newcomments - new comments (string or cell array of strings)
%                 to assign (during commandline calls only)
%   concat      - [0|1] 1 concatenate the newcomments to the old one.
%                 Default is 0.
%
% Outputs:
%   newcomments - new comments, string
%
% Note: if new comments are given as input, there are simply
%       converted and returned by the function; otherwise a
%       window pops up.
%
% Example
%  EEG.comments = pop_comments( { 'This is the first line.' ' ' ...
%               'This is the third line.' }, 'Editing');
% EEG.comments = pop_comments(EEG.comments,'','This is the fourth line.",1);
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: not supported by cvs2svn $
% Revision 1.15  2005/11/04 17:56:02  arno
% fixing closing several windows
%
% Revision 1.14  2005/09/27 22:01:28  arno
% fix multiline text for windows
%
% Revision 1.13  2004/05/14 18:24:00  arno
% fixing double quote problem
%
% Revision 1.12  2004/05/14 17:16:18  hilit
% chnaged the help message
%
% Revision 1.11  2003/05/09 22:04:46  arno
% making the comments incremental
%
% Revision 1.10  2003/04/09 23:28:12  arno
% debuging command line call
%
% Revision 1.9  2002/08/17 01:03:34  scott
% font change, OK -> SAVE
%
% Revision 1.8  2002/08/12 01:23:16  arno
% update colors
%
% Revision 1.7  2002/08/12 01:22:20  arno
% updating colors
%
% Revision 1.6  2002/04/10 19:21:28  arno
% debuging cancel button
%
% Revision 1.5  2002/04/10 19:17:12  arno
% setting title propertie to no interpreter
%
% Revision 1.4  2002/04/08 20:52:45  arno
% removing debug msg
%
% Revision 1.3  2002/04/06 02:58:37  arno
% returning [] when no modification
%
% Revision 1.2  2002/04/06 02:49:39  arno
% reprogrammed the whole function
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 03-16-02 text interface editing -sm & ad 

function [newcomments, com] = pop_comments( comments, plottitle, newcomments, concat );

com = '';
if exist('comments') ~=1, comments = '';
elseif iscell(comments), comments = strvcat(comments{:}); 
end;

% remove trailing blanks and make multiline
comments = strmultiline( comments, 53);

if nargin < 3
    newcomments = comments;
	try, icadefs;
	catch,
		BACKCOLOR  =  [.8 .8 .8];     
		GUIBUTTONCOLOR   = [.8 .8 .8];    
	end;
	figure('menubar', 'none', 'tag', 'comment', 'color', BACKCOLOR, 'userdata', 0, ...
		   'numbertitle', 'off', 'name', 'Read/Enter comments -- pop_comments()');
	pos = get(gca,'position'); % plot relative to current axes
	q = [pos(1) pos(2) 0 0];
	s = [pos(3) pos(4) pos(3) pos(4)]./100;
	if exist('plottitle') ~=1, plottitle = ''; end;
	
	h = title(plottitle);
	set(h, 'fontname','Helvetica','fontweight', 'bold', 'interpreter', 'none');
	
	axis off;

	% create the buttons
	% ------------------
  	uicontrol('Parent',gcf, ...
  	'Units','Normalized', ...
	'Position', [0 -5 20 10].*s+q, ...
	'backgroundcolor', GUIBUTTONCOLOR, ...
	'string','CANCEL', 'callback', 'close(findobj(''tag'', ''comment''));' );
		
  	uicontrol('Parent',gcf, ...
  	'Units','Normalized', ...
	'Position', [80 -5 20 10].*s+q, ...
	'backgroundcolor', GUIBUTTONCOLOR, ...
	'string','SAVE', 'callback', ...
		[ 'set(gcbf, ''userdata'', ' ...
		'get(findobj(''parent'', gcbf, ''tag'', ''edit''), ''string''));' ]);

	%hh = text( q(1), 100*s(2)+q(2), comments, 'tag', 'edit');
	%set( hh, 'editing', 'on', 'verticalalignment', 'top');

    %hh = uicontrol('Parent',gcf, ...
  	%'Units','Normalized', ...
  	%'style', 'text', ...
	%'Position', [0 100 105 5].*s+q, ...
	%'string', 'Warning: each blank line must contain at least a ''space'' character', ...
	%'horizontalalignment', 'left', ...
    %'backgroundcolor', BACKCOLOR );

    hh = uicontrol('Parent',gcf, ...
  	'Units','Normalized', ...
  	'style', 'edit', ...
  	'tag', 'edit', ... 
	'Position', [0 10 105 85].*s+q, ...
	'string', comments, ...
	'backgroundcolor', [ 1 1 1], ...
	'horizontalalignment', 'left', ...
	'max', 3, ...
	'fontsize', 12);

    % Try to use 'courier' since it has constant character size
    lf = listfonts;
    tmppos = strmatch('Courier', lf);
    if ~isempty(tmppos)
        set(hh, 'fontname', lf{tmppos(1)}, 'fontsize', 10);
    end;
    
    waitfor(gcf, 'userdata');

    % find return mode
    if isempty(get(0, 'currentfigure')), return; end;
    tmp = get(gcf, 'userdata');
    if ~isempty(tmp) & isstr(tmp)    
        newcomments = tmp; % ok button
    else return;
    end;

	close(findobj('tag', 'comment'));
else
    if iscell(newcomments)
        newcomments = strvcat(newcomments{:});
    end;
    if nargin > 3 & concat == 1
        newcomments = strvcat(comments, newcomments);
    end;
    return;
end;

I = find( comments(:) == '''');
comments(I) = ' ';  
if nargout > 1
        if ~strcmp( comments, newcomments)
          allsame = 1;
            for index = 1:size(comments, 1)
                if ~strcmp(comments(index,:), newcomments(index,:)), allsame = 0; end;
            end;
        else
            allsame = 0;
        end;
        if allsame & ~isempty(EEG.comments)
             com =sprintf('EEG.comments = pop_comments(EEG.comments, '''', %s, 1);', vararg2str(newcomments(index+1:end,:)));
        else 
            com =sprintf('EEG.comments = pop_comments('''', '''', %s);', vararg2str(newcomments));     
        end;
end;
return;
