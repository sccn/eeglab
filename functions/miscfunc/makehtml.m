% makehtml() - function generating all the HTML files for files in specific 
%              directories 
%
% Usage: 
%   >> makehtml( directorylist, outputdir, titlelist, matlabhelp);
%   >> makehtml( directorylist, outputdir, titlelist, matlabhelp, ...
%               'key1', val1, 'key2', val2, ...);
%
% Input:
%    directorylist - list (cell array) of directory names. if nested lists,
%                    the functions of the directories in the nested lists
%                    are regroupeg 
%    outputdir     - output HTML directory
%    titlelist     - list of title for each directory
%    matlabhelp    - list of filenames to generate matlab help pop-up index.
%                    Default []: no help files generated   
%
% Optional inputs:
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
%   
% Author: Arnaud Delorme, CNL / Salk Institute, 2002
%
% Example: 
%  makehtml({ 'adminfunc', 'popfunc', { 'toolbox', 'toolbox2' }}, ...
%   '/home/www/eeglab/eeglab/allfunctions', ...
%   { 'Admin functions', 'Interactive (pop_) functions', 'Signal processing functions' }, ...
%   { 'adminfunc/eeg_helpadmin.m', 'adminfunc/eeg_helppop.m', 'adminfunc/eeg_helpsigproc.m' });          
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
% Revision 1.1  2002/04/05 17:39:45  jorn
% Initial revision
%

function makehtml( directorylist, outputdir, titlelist, matlabhelp, varargin );

if nargin < 3
    help makehtml;
    return;
end;
    
if nargin < 4
    matlabhelp = [];
else
    if length(directorylist) ~= length( matlabhelp ),
        error( 'the directory list and the matlab help file list must have the same length');
    end;    
end;

if length(directorylist) ~= length( titlelist ),
     error( 'the directory list and the title list must have the same length');
end;
    
if outputdir(end) ~= '/', outputdir(end+1) = '/'; end; 
    
if ~isempty( varargin )
    g = struct( varargin{:} );
else
    g = [];
end;
    
try, g.fontindex;   catch, g.fontindex = 'Helvetica'; end;
try, g.backindex;   catch, g.backindex = '<body BACKGROUND="cream_stucco.jpg" bgproperties="fixed" bgcolor="#ffffe5">'; end;
try, g.header;      catch, g.header = [ '<script language="JavaScript"><!--' 10 'function openhelp(fnc){' 10 'self.window.location = fnc;' 10 '}' 10  '//--></script>' ]; end;
try, g.background;  catch, g.background = '<body BACKGROUND="cream_stucco.jpg" bgproperties="fixed" bgcolor="#ffffe5">'; end;
try, g.refcall;     catch, g.refcall = 'javascript:openhelp(''%s.html'')'; end;
try, g.font;        catch, g.font = 'Helvetica'; end;
try, g.footer;      catch, g.footer = '<A HREF ="indexfunc.html">Back to functions</A>';  end;
try, g.outputlink;  catch, g.outputlink = [ '<tr><td VALIGN=TOP ALIGN=RIGHT NOSAVE><A HREF="javascript:openhelp(''%s.html'')">%s</A></td><td>%s</td></tr>' ];  end;

options = { 'footer', g.footer, 'background', g.background, ...
		  'refcall', g.refcall, 'font', g.font, 'header', g.header, 'outputlink', g.outputlink };

% ------------------------------------------- 
% scrips which generate a web page for eeglab
% ------------------------------------------- 
STYLEHEADER = '<BR><H2>%s</H2>\n'; 
OPENWIN = [ '<script language="JavaScript"><!--' 10 'function openhelp(fnc){' 10 ... 
            'self.window.location = fnc;' 10 '}' 10 '//--></script>'];
ORIGIN      = pwd;

% get files
% ---------
if ~isempty( matlabhelp)
    for index = 1:length( matlabhelp )
        try, delete(matlabhelp{index}); catch, end; % delete previous files   
    end;
end;   

% scan directories
% ----------------
for index = 1:length( directorylist )
    direct{ index } = scandir( directorylist{index} );
end;    

% write HTML file
% ----------------
fo = fopen([ outputdir 'indexfunc.html'], 'w');
fprintf(fo, '<HTML><HEAD>%s</HEAD>%s<FONT FACE="%s">\n', OPENWIN, g.backindex, g.fontindex);

for index = 1:length( direct )
    makehelphtml( direct{ index }, fo, titlelist{ index }, STYLEHEADER, outputdir, options );
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
if ~isempty( matlabhelp )
    for index = 1:length( matlabhelp )
        makehelpmatlab( matlabhelp{index}, direct{ index }, titlelist{ index }); 
    end;
end;    

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
function makehelphtml( files, fo, title, STYLEHEADER, DEST, options);
        tmpdir = pwd;
	fprintf(fo, STYLEHEADER, title );
	fprintf(fo, '<table WIDTH="100%%" NOSAVE>' );
	for index = 1:length(files)
		fprintf('Processing %s\n', files{index});
		cd(DEST); com = help2html( files{index}, [], options{:}); cd(tmpdir);
		fprintf( fo, '%s', com);
		inputfile = which( files{index});
		try, copyfile( inputfile, [ DEST files{index} ]); % asuming the file is in the path 
        catch, fprintf('Cannot copy file %s\n', inputfile); end;
	end;	
	fprintf(fo, '</table>' );
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
	fprintf(fo, 'return;\n');
	fclose( fo );
return

