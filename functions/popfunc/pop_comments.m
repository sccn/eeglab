% pop_comments() - Ddd comments to an EEG dataset
%
% Usage:
%   >> EEGOUT = pop_comments( EEGIN, comments );
%
% Inputs:
%   EEGIN      - input dataset
%   comment    - character array
%
% Outputs:
%   EEGOUT     - output dataset
%
% Example
%   EEG = pop_comments( EEG, { 'This is the first line' '' ...
%               'this is the third line' });  
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

% 01-25-02 reformated help & license -ad 
% 03-16-02 text interface editing -sm & ad 

function [EEG, com] = pop_comments( EEG, comments );

com = '';
if nargin < 1
	help pop_comments;
	return;
end;	
if nargin < 2
	figure('menubar', 'none', 'numbertitle', 'off', 'name', 'About this dataset -- pop_comment()');
	pos = get(gca,'position'); % plot relative to current axes
	q = [pos(1) pos(2) 0 0];
	s = [pos(3) pos(4) pos(3) pos(4)]./100;
	set(gcf, 'userdata', 0);
	title('Comments on current dataset');

	axis off;

	try, comments = EEG.comments;
	catch, comments = [];
	end;

	% create the buttons
	% ------------------
  	uicontrol('Parent',gcf, ...
  	'Units','Normalized', ...
	'Position', [0 0 20 10].*s+q, ...
	'string','CANCEL', 'callback', ...
		[ 'set(gcbf, ''userdata'', -1);' ]);
		
  	uicontrol('Parent',gcf, ...
  	'Units','Normalized', ...
	'Position', [85 0 20 10].*s+q, ...
	'string','OK', 'callback', ...
		[ 'set(gcbf, ''userdata'', ' ...
		'get(findobj(''parent'', gcbf, ''tag'', ''edit''), ''string''));' ]);


	%hh = text( q(1), 100*s(2)+q(2), comments, 'tag', 'edit');
	%set( hh, 'editing', 'on', 'verticalalignment', 'top');
	
  	hh = uicontrol('Parent',gcf, ...
  	'Units','Normalized', ...
  	'style', 'edit', ...
  	'tag', 'edit', ... 
	'Position', [0 15 100 85].*s+q, ...
	'string', comments, ...
	'backgroundcolor', [ 1 1 1], ...
	'horizontalalignment', 'left', ...
	'max', 2, ...
	'fontsize', 12);

	waitfor(gcf, 'userdata');

	if isstr(get(gcf, 'userdata'))
		comments = get(gcf, 'userdata'); % ok button
	end;		

	close(gcf);
else
	try, comments = char(eval(comments));	catch, end; 
end;	
 
EEG.comments = comments;

I = find( comments(:) == '''');
comments(I) = ' ';  

com =sprintf('%s = pop_comments( %s, ''%s'');', inputname(1), inputname(1), array2str(comments));
return;
 
function str = array2str( array )
	str = '[';
	for index = 1:size(array,1)
		str = [ str '; [' num2str(double(array(index,:))) '] ' ];
	end;
	str = [ str ']' ];
return;
		 

