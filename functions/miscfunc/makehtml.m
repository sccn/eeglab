% makehtml() - function generating all the HTML files for files in specific 
%              directories 
%
% Usage: 
%   >> makehtml( list, outputdir);
%   >> makehtml( list, outputdir, ...
%               'key1', val1, 'key2', val2, ...);
%
% Input:
%    list        - 1) list (cell array) of filenames to convert. Ex:
%                     { 'filename1' 'filename2' 'filename3' }
%                  2) cell array of filenames to convert and the text link 
%                     on the summary page for them. For option 1, as a default
%                     the filename without the extension is used.
%                    { { 'filename1' 'link1' } ...
%                      { 'filename2' 'link2' } }
%                  3) cell array of 2-3 cell array elements containing
%                  info to generate directory HTML and matlab page. Ex
%                  { { 'directory1' 'title1' 'matlabfunc1' } ...
%                    { 'directory1'''title1' 'matlabfunc2' } }
%                  'directory' is the directory name
%                  'title' is the title that will be given the directory on
%                    the web page.
%                  'matlabfunc' is the name for a matlab function that will be
%                    created to sumarized the directory content (note that only
%                    one web page is created by several matlab files can be created
%                    It can be left empty.
%                  To scan several direcotry under the same title use convention
%                  { { { 'directory1' 'direcotory2 } 'title1' 'matlabfunc1' }
%                  List can also contains filenames
%    outputdir   - output HTML directory
%
% Optional inputs:
%   'outputfile' - output file name. Default is 'index.html'.
%   'header'     - command to insert in the header of all HTML files (i.e. javascript 
%                  declaration or meta-tags). Default: javascript 'openhelp' function. 
%   'footer'     - command to insert at the end of all HTML files (i.e. back
%                  button. Default is a reference back to the index file. 
%   'refcall'    - syntax format to call references. Default javascript function 
%                  uses 'javascript:openhelp(''%s.js'')'. Use '%s.html' for
%                  standard HTML link.     
%   'font'       - font name (default: 'Helvetica')
%   'background' - baground tag (with body tag i.e. '<body BACKGROUND="img.jpg">'). 
%   'outputlink' - command to call the generated HTML page. Default is
%                  'javascript:openhelp(''%s.js'')'. Use '%s.html' for
%                  standard HTML link.     
%   'fontindex'  - font of the HTML index file (default: 'Helvetica')
%   'backindex'  - background tag for the index file
%   'mainonly'   - ['on'|'off'] only generate main page. Default is OFF and
%                  generate all web pages.
%   
% Author: Arnaud Delorme, CNL / Salk Institute, 2002
%
% Example used for generating EEGLAB help menus:
%  makehtml({ { 'adminfunc' 'Admin functions' 'adminfunc/eeg_helpadmin.m' } ...
%             { 'popfunc', 'Interactive (pop_) functions' 'adminfunc/eeg_helppop.m' } ...
%             { { 'toolbox', 'toolbox2' } 'Signal processing functions' 'adminfunc/eeg_helpsigproc.m' }}, ...
%            '/home/www/eeglab/allfunctions');          
%
% See also: help2html()

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 2002
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
% Revision 1.6  2002/08/17 20:03:52  arno
% new background color
%
% Revision 1.5  2002/08/17 00:42:33  arno
% new modifs
%
% Revision 1.4  2002/08/16 16:29:36  arno
% new revision
%
% Revision 1.3  2002/08/16 15:12:21  arno
% after crash
%
% Revision 1.2  2002/04/06 01:58:16  arno
% debugging destination of html files
%
% Revision 1.1  2002/04/05 17:39:45  jorn
% Initial revision
%

function makehtml( directorylist, outputdir, varargin );

if nargin < 2
    help makehtml;
    return;
end;
    
if outputdir(end) ~= '/', outputdir(end+1) = '/'; end; 
    
if ~isempty( varargin )
    g = struct( varargin{:} );
else
    g = [];
end;
    
try, g.mainonly;    catch, g.mainonly = 'off'; end;
try, g.outputfile;  catch, g.outputfile = 'index.html'; end;
try, g.fontindex;   catch, g.fontindex = 'Helvetica'; end;
try, g.backindex;   catch, g.backindex =  '<body bgcolor="#fcffff">'; end;
try, g.header;      catch, g.header = [ '<script language="JavaScript"><!--' 10 'function openhelp(fnc){' 10 'self.window.location = fnc;' 10 '}' 10  '//--></script>' ]; end;
try, g.background;  catch, g.background = '<body bgcolor="#fcffff">'; end;
%try, g.background;  catch, g.background = '<body BACKGROUND="cream_stucco.jpg" bgproperties="fixed" bgcolor="#ffffe5">'; end;
try, g.refcall;     catch, g.refcall = 'javascript:openhelp(''%s.html'')'; end;
try, g.font;        catch, g.font = 'Helvetica'; end;
try, g.footer;      catch, g.footer = '<A HREF ="indexfunc.html">Back to functions</A>';  end;
try, g.outputlink;  catch, g.outputlink = [ '<tr><td VALIGN=TOP ALIGN=RIGHT NOSAVE><A HREF="javascript:openhelp(''%s.html'')">%s</A></td><td>%s</td></tr>' ];  end;

options = { 'footer', g.footer, 'background', g.background, ...
		  'refcall', g.refcall, 'font', g.font, 'header', g.header, 'outputlink', g.outputlink, 'outputonly', g.mainonly };

% ------------------------------------------- 
% scrips which generate a web page for eeglab
% ------------------------------------------- 
STYLEHEADER = '<BR><H2>%s</H2>\n'; 
OPENWIN = [ '<script language="JavaScript"><!--' 10 'function openhelp(fnc){' 10 ... 
            'self.window.location = fnc;' 10 '}' 10 '//--></script>'];
ORIGIN      = pwd;

% determine mode
% --------------
if iscell(directorylist{1}) & exist(directorylist{1}{1}) == 7
	fprintf('First cell array element is not a file\n');
	fprintf('Scanning directories...\n');
	mode = 'dir';
	% scan directories
	% ----------------
	for index = 1:length( directorylist )
		direct{ index } = scandir( directorylist{index}{1} );
	end;    
else
	fprintf('First cell array element has been identified as a file\n');
	fprintf('Scanning all files...\n');
	mode = 'files';
end;	

% remove . directory
% ------------------
rmpath('.');

% write HTML file
% ----------------
fo = fopen([ outputdir g.outputfile], 'w');
if fo == -1, error(['can not open file ''' [ outputdir g.outputfile] '''']); end;

fprintf(fo, '<HTML><HEAD>%s</HEAD>%s<FONT FACE="%s">\n', OPENWIN, g.backindex, g.fontindex);

if strcmp(mode, 'files')
	makehelphtml( directorylist, fo, 'MAIN TITLE', STYLEHEADER, outputdir, mode, options, g.mainonly );
else % direcotry
	for index = 1:length( directorylist )
		makehelphtml( direct{ index }, fo, directorylist{index}{2}, STYLEHEADER, outputdir, mode, options, g.mainonly );
	end;
end;	
fprintf( fo, '</FONT></BODY></HTML>');
fclose( fo );
if isunix
    chmodcom = sprintf('!chmod 777 %s*', outputdir);
    eval(chmodcom);
end;

% ------------------------------
% Generate help files for EEGLAB
% ------------------------------
if strcmp(mode, 'dir')
	for index = 1:length( directorylist )
		if length(directorylist{index}) > 2
			makehelpmatlab( directorylist{index}{3}, direct{ index },directorylist{index}{2}); 
		end;    
	end;
end;
addpath('.');	


return;

% scan directory list or directory
% --------------------------------
function filelist = scandir( dirlist )
    filelist = {};
    if iscell( dirlist )
        for index = 1:length( dirlist )
            tmplist = scandir( dirlist{index} );
            filelist  = { filelist{:} tmplist{:} };
        end;    
    else
        if dirlist(end) ~= '/', dirlist(end+1) = '/'; end;
        if exist(dirlist) ~= 7
            error([ dirlist ' is not a directory']);
        end;    
        tmpdir  =  dir([dirlist '*.m']); 
        filelist = { tmpdir(:).name }; 
    end;
    filelist = sort( filelist );      
return;

% ------------------------------ Function to generate help for a bunch of files -
function makehelphtml( files, fo, title, STYLEHEADER, DEST, mode, options, mainonly);
% files = cell array of string containing file names or
%         cell array of 2-strings cell array containing titles and filenames
% fo = output file
% title = title of the page or section	
	tmpdir = pwd;
	if strcmp(mode, 'files') % processing subtitle and File
		fprintf(fo, '<UL>' );
		for index = 1:length(files)
			if iscell(files{index})
				filename = files{index}{1};
			    filelink = files{index}{2};
			else
				filename = files{index};
			    filelink = '';
			end;
			fprintf('Processing %s:%s\n', filename, filelink );
			if ~isempty(filename)				
				cd(DEST); help2html( filename, [],  'outputtext', filelink, options{:}); cd(tmpdir);
				
				if strcmp(mainonly,'off')
					inputfile = which( filename);
					try, copyfile( inputfile, [ DEST filename ]); % asuming the file is in the path 
					catch, fprintf('Cannot copy file %s\n', inputfile); end;
				end;
				
				indexdot = find(filename == '.');
				if ~isempty(filelink)
					com = [ space2html(filelink)  ' -- ' space2html([ filename(1:indexdot(end)-1) '()'], ...
								 [ '<A HREF="' filename(1:indexdot(end)-1) '.html">' ], '</A><BR>')];
				else
					com = [ space2html([ filename(1:indexdot(end)-1) '()'], ...
								 [ '<A HREF="' filename(1:indexdot(end)-1) '.html">' ], '</A><BR>')];
				end;
			else 
				com = space2html(filelink, '<B>', '</B><BR>');
			end;
			fprintf( fo, '%s', com);
		end;
		fprintf(fo, '</UL>' );
	else 
		fprintf(fo, STYLEHEADER, title );
		fprintf(fo, '<table WIDTH="100%%" NOSAVE>' );
		for index = 1:length(files)
	        % Processing file only
 			fprintf('Processing %s\n', files{index});
			cd(DEST); com = help2html( files{index}, [], options{:}); cd(tmpdir);
			fprintf( fo, '%s', com);
			if strcmp(mainonly,'off')
				inputfile = which( files{index});
				try, copyfile( inputfile, [ DEST files{index} ]); % asuming the file is in the path 
				catch, fprintf('Cannot copy file %s\n', inputfile); end;
			end;
		end;
		fprintf(fo, '</table>' );
	end;	
return;

% ------------------------------ Function to pop-out a Matlab help window --------
function makehelpmatlab( filename, directory, titlewindow);
	fo = fopen( filename, 'w');
	fprintf(fo, '%%%s() - Help file for EEGLAB\n\n', filename);
	fprintf(fo, 'function noname();\n\n');
	fprintf(fo, 'command = { ...\n');
	for index = 1:length(directory)
		fprintf(fo, '''pophelp(''''%s'''');'' ...\n', directory{index});
	end;	
	fprintf(fo, '};\n');
	fprintf(fo, 'text = { ...\n');
	for index = 1:length(directory)
		fprintf(fo, '''%s'' ...\n', directory{index});
	end;	
	fprintf(fo, '};\n');
	fprintf(fo, ['textgui( text, command,' ...
	'''fontsize'', 15, ''fontname'', ''times'', ''linesperpage'', 18, ', ...
	'''title'',strvcat( ''%s'', ''(Click on blue text for help)''));\n'], titlewindow);
	fprintf(fo, 'icadefs; set(gcf, ''COLOR'', BACKCOLOR);');
	fprintf(fo, 'h = findobj(''parent'', gcf, ''style'', ''slider'');');
	fprintf(fo, 'set(h, ''backgroundcolor'', GUIBACKCOLOR);');
	fprintf(fo, 'return;\n');
	fclose( fo );
return

% convert spaces for HTML
function strout = space2html(strin, linkb, linke)
	strout = [];
	index = 1;
	while strin(index) == ' '
		strout = [ strout '&nbsp; '];
		index = index+1;
	end;
	if nargin == 3
		strout = [strout linkb strin(index:end) linke];
	else
		strout = [strout strin(index:end)];
	end;
				   