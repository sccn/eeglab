% pophelp() - Same as matlab HTHELP but does not crash under windows.
%
% Usage: >> pophelp( function );
%        >> pophelp( function, nonmatlab );
%
% Inputs:
%   function  - string for a Matlab function name 
%               (with or without the '.m' extension).
%   nonmatlab - [0|1], 1 the file is not a Matlab file
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eeglab() 

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
% Revision 1.10  2003/03/12 03:12:34  arno
% special help for pop_functions
%
% Revision 1.9  2003/03/12 02:59:32  arno
% adding pop and eponymous help
%
% Revision 1.8  2002/11/15 02:47:39  arno
% header for web
%
% Revision 1.7  2002/10/23 15:05:15  arno
% isppc -> computer
%
% Revision 1.6  2002/10/15 17:15:15  arno
% windows -> removing last function character
%
% Revision 1.5  2002/08/13 15:26:44  arno
% updating colors
%
% Revision 1.4  2002/07/29 15:53:09  arno
% same
%
% Revision 1.3  2002/07/29 15:51:06  arno
% debugging
%
% Revision 1.2  2002/07/29 15:50:15  arno
% can also read non-matlab files
%
% Revision 1.1  2002/04/05 17:39:45  jorn
% Initial revision
%
% 01-25-02 reformatted help & license -ad 
% 3/19/02 Edited help message (below) -sm 

function pophelp( funct, nonmatlab );

if nargin <1
	help pophelp;
	return;
end;
if nargin <2
	nonmatlab = 0;
end;

doc1 = readfunc(funct, nonmatlab);
if length(funct) > 4 & strcmpi(funct(1:4), 'pop_')
	try,
		doc2 = readfunc(funct(5:end), nonmatlab);
		doc1 = { doc1{:} ' _________________________________________________________________ ' ...
					   ' ' ...
                       ' The ''pop'' function above calls the eponymous Matlab function below, ' ...
                       ' which may contain more information for some parameters. '...
					   ' ' ...
					   ' _________________________________________________________________ ' ...
                       ' ' ...
				doc2{:} };
	catch, end;
end;

textgui(doc1);
h = findobj('parent', gcf, 'style', 'slider');
try, icadefs; catch, 
	GUIBUTTONCOLOR = [0.8 0.8 0.8]; 
	GUITEXTCOLOR   = 'k'; 
end;
set(h, 'backgroundcolor', GUIBUTTONCOLOR);
h = findobj('parent', gcf, 'style', 'pushbutton');
set(h, 'backgroundcolor', GUIBUTTONCOLOR);
h = findobj('parent', gca);
set(h, 'color', GUITEXTCOLOR);
set(gcf, 'color', BACKCOLOR);

return;

function [doc] = readfunc(funct, nonmatlab)

doc = {};
if nonmatlab	
	fid = fopen( funct, 'r');
else
	if findstr( funct, '.m')
		fid = fopen( funct, 'r');
	else
		fid = fopen( [funct '.m'], 'r');
	end;
end;

if fid == -1
	error('File not found');
end;

sub = 1;
try, 
    if strcmpi(computer, 'PCWIN') | strcmp(computer, 'MAC'), sub = 0; end;
catch, end;

if nonmatlab
	str = fgets( fid );
	while ~feof(fid)
		str = deblank(str(1:end-sub));
        doc = { doc{:} str(1:end) };    
        str = fgets( fid );
	end;
else
	str = fgets( fid );
	while (str(1) == '%')
		str = deblank(str(1:end-sub));
        doc = { doc{:} str(2:end) };    
		str = fgets( fid );
	end;
end;
