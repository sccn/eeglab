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

function pophelp( funct, nonmatlab );

if nargin <1
	help pophelp;
	return;
end;
if nargin <2
	nonmatlab = 0;
end;

if exist('help2html')
    if length(funct) > 3 && strcmpi(funct(end-3:end), '.txt')
        web(funct);
    else
        text1 = help2html(funct);
        if length(funct) > 4 & strcmpi(funct(1:4), 'pop_')
            try,
                text2 = help2html(funct(5:end));
                text1 = [text1 '<br><pre>___________________________________________________________________' 10 ...
                               ' ' 10 ...
                               ' The ''pop'' function above calls the eponymous Matlab function below' 10 ...
                               ' and could use some of its optional parameters' 10 ...
                               '___________________________________________________________________</pre><br><br>' text2 ];
            catch, end;
        end;

        web([ 'text://' text1 ]);
    end;
else
    if isempty(funct), return; end;
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

    textgui(doc1);1000
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
end;
return;

function [doc] = readfunc(funct, nonmatlab)

doc = {};
if iseeglabdeployed
    if isempty(find(funct == '.')), funct = [ funct '.m' ]; end;
    funct = fullfile(eeglabexefolder, 'help', funct);
end;
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
    if ~isunix, sub = 0; end;
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
fclose(fid);
