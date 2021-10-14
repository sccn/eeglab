% makehtml() - generate .html function-index page and function help pages 
%              composed automatically from formatted Matlab function help messages
%
% Usage: 
%   >> makehtml(list, outputdir);
%   >> makehtml(list, outputdir, 'key1', val1, 'key2', val2, ...);
%
% Input:
%    list        - (1) List (cell array) of filenames to convert. 
%                     Ex: {'filename1' 'filename2' 'filename3'} 
%                     By default, the filename extension .m is understood.
%                  (2) Cell array of filenames to convert and the text link 
%                     on the summary page for them.
%                      {{'filename1' 'link1'} {'filename2' 'link2'}} ...
%                     Ex: 'link1' = 'Reject by kurtosis'. 
%                  (3) Cell array of 2 or 3 cell array elements containing
%                     info to generate .html function-index and help pages. 
%                   Ex: { {'directory1' 'heading1' 'dirindexfunc1'} ...
%                         {'directory2' 'heading2' 'dirindexfunc2'} }
%                   - 'directory': Function file directory name
%                   - 'heading': Index-file heading for the directory functions 
%                   - 'dirindexfunc': A optional Matlab pop-up help function for
%                     referencing the directory functions {default: none}
%                  (4) To scan several directories under the same heading, use 
%                  {{{'directory1' 'directory2'} 'heading1' 'dirindexfunc1' ... }
%    outputdir   - Directory for output .html help files 
%
% Optional inputs:
%   'outputfile' - Output file name. {default: 'index.html'}
%   'header'     - Command to insert in the header of all .html files (e.g., javascript 
%                  declaration or meta-tag). {default: javascript 'openhelp()' 
%                  function. See help2htm() code for details.}
%   'footer'     - Command to insert at the end of all .html files (e.g., back
%                  button. {default: reference back to the function-index file}
%   'refcall'    - Syntax format to call references. {default is  
%                  'javascript:openhelp(''%s.js'')'} Use '%s.html' for an .html link.     
%   'font'       - Font name (default: 'Helvetica')
%   'background' - Background HTML body section. Include "<" and ">". 
%   'outputlink' - Help page calling command for the function index page. {default is
%                  'javascript:openhelp(''%s.js'')'.  Use '%s.html' to use a
%                  standard .html page link instead.}    
%   'fontindex'  - Font for the .html index file (default: 'Helvetica')
%   'backindex'  - Background tag for the index file (c.f. 'background')
%   'mainheader' - Text file to insert at the beginning of the index page. Default is
%                  none.
%   'mainonly'   - ['on'|'off'] 'on' -> Generate the index page only.
%                   {default: 'off' -> generate both the index and help pages}
%   
% Author: Arnaud Delorme, CNL / Salk Institute, 2002
%
% Example: Generate EEGLAB help menus at SCCN:
%  makehtml({ { 'adminfunc' 'Admin functions' 'adminfunc/eeg_helpadmin.m' } ...
%             { 'popfunc', 'Interactive pop_functions' 'adminfunc/eeg_helppop.m' } ...
%             { { 'toolbox', 'toolbox2' } 'Signal processing functions' 'adminfunc/eeg_helpsigproc.m' }}, ...
%            '/home/www/eeglab/allfunctions', 'mainheader', '/data/common/matlab/indexfunchdr.txt');          
%
% See also: help2html2()

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 2002
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

function makehtml( directorylist, outputdir, varargin );

if nargin < 2
    help makehtml;
    return;
end
    
if outputdir(end) ~= '/', outputdir(end+1) = '/'; end; 
    
if ~isempty( varargin )
    g = struct( varargin{:} );
else
    g = [];
end
    
try, g.mainonly;    catch, g.mainonly = 'off'; end
try, g.mainheader;  catch, g.mainheader = ''; end
try, g.outputfile;  catch, g.outputfile = 'index.html'; end
try, g.fontindex;   catch, g.fontindex = 'Helvetica'; end
try, g.backindex;   catch, g.backindex =  '<body bgcolor="#fcffff">'; end
try, g.header;      catch, g.header = [ '<script language="JavaScript"><!--' 10 'function openhelp(fnc){' 10 'self.window.location = fnc;' 10 '}' 10  '//--></script>' ]; end
try, g.background;  catch, g.background = '<body bgcolor="#fcffff">'; end
%try, g.background;  catch, g.background = '<body BACKGROUND="cream_stucco.jpg" bgproperties="fixed" bgcolor="#ffffe5">'; end
try, g.refcall;     catch, g.refcall = 'javascript:openhelp(''%s.html'')'; end
try, g.font;        catch, g.font = 'Helvetica'; end
try, g.footer;      catch, g.footer = '<A HREF ="index.html">Back to functions</A>';  end
try, g.outputlink;  catch, g.outputlink = [ '<tr><td VALIGN=TOP ALIGN=RIGHT NOSAVE><A HREF="javascript:openhelp(''%s.html'')">%s</A></td><td>%s</td></tr>' ];  end

% read header text file
% ---------------------
if ~isempty(g.mainheader)
    doc = [];
    fid = fopen(g.mainheader , 'r');
    if (fid == -1), error(['Can not open file ''' g.mainheader '''' ]); end
    str = fgets( fid );
    while ~feof(fid)
        str = deblank(str(1:end-1));
        doc = [ doc str(1:end) ];
        str = fgets( fid );
    end
    g.backindex = [ g.backindex doc ];
end

options = { 'footer', g.footer, 'background', g.background, ...
		  'refcall', g.refcall, 'font', g.font, 'header', g.header, 'outputlink', g.outputlink};
if strcmpi( g.mainonly, 'on')
    options = { options{:}, 'outputonly', g.mainonly };
end

% ------------------------------------------- 
% scrips which generate a web page for eeglab
% ------------------------------------------- 
STYLEHEADER = '<BR><a name="%s"></a><H2>%s</H2>\n'; 
OPENWIN = [ '<script language="JavaScript"><!--' 10 'function openhelp(fnc){' 10 ... 
            'self.window.location = fnc;' 10 '}' 10 '//--></script>'];
ORIGIN      = pwd;

% determine mode
% --------------
if iscell(directorylist{1}) && exist(directorylist{1}{1}) == 7
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

% write .html file
% ----------------
fo = fopen([ outputdir g.outputfile], 'w');
if fo == -1, error(['cannot open file ''' [ outputdir g.outputfile] '''']); end

fprintf(fo, '<HTML><HEAD>%s</HEAD>%s<FONT FACE="%s">\n', OPENWIN, g.backindex, g.fontindex);

if strcmp(mode, 'files')
	makehelphtml( directorylist, fo, 'MAIN TITLE', STYLEHEADER, outputdir, mode, options, g.mainonly );
else % directory
	for index = 1:length( directorylist )
		makehelphtml( direct{ index }, fo, directorylist{index}{2}, STYLEHEADER, outputdir, mode, options, g.mainonly );
	end
end;	
fprintf( fo, '</FONT></BODY></HTML>');
fclose( fo );
if isunix
    chmodcom = sprintf('!chmod 777 %s*', outputdir);
    eval(chmodcom);
end

% ------------------------------
% Generate help files for EEGLAB
% ------------------------------
if strcmp(mode, 'dir')
	for index = 1:length( directorylist )
		if length(directorylist{index}) > 2
			makehelpmatlab( directorylist{index}{3}, direct{ index },directorylist{index}{2}); 
		end;    
	end
end
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
        if dirlist(end) ~= '/', dirlist(end+1) = '/'; end
        if exist(dirlist) ~= 7
            error([ dirlist ' is not a directory']);
        end;    
        tmpdir  =  dir([dirlist '*.m']); 
        filelist = { tmpdir(:).name }; 
    end
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
			end
			fprintf('Processing (mode file) %s:%s\n', filename, filelink );
			if ~isempty(filename)
                if ~exist(fullfile(DEST, [ filename(1:end-1) 'html' ]))
                    cd(DEST); 
                    try, delete([ DEST filename ]);
                    catch, end
                    help2html2( filename, [],  'outputtext', filelink, options{:}); cd(tmpdir);
                    
                    if strcmp(mainonly,'off')
                        inputfile = which( filename);
                        try, copyfile( inputfile, [ DEST filename ]); % assuming the file is in the path 
                        catch, fprintf('Cannot copy file %s\n', inputfile); end
                    end
                    
                    indexdot = find(filename == '.');
                end
                if ~isempty(filelink)
                    com = [ space2html(filelink)  ' -- ' space2html([ filename(1:indexdot(end)-1) '()'], ...
                                                                    [ '<A HREF="' filename(1:indexdot(end)-1) '.html">' ], '</A><BR>')];
                else
                    com = [ space2html([ filename(1:indexdot(end)-1) '()'], ...
                                       [ '<A HREF="' filename(1:indexdot(end)-1) '.html">' ], '</A><BR>')];
                end
			else 
				com = space2html(filelink, '<B>', '</B><BR>');
			end
			fprintf( fo, '%s', com);
		end
		fprintf(fo, '</UL>' );
	else 
		fprintf(fo, STYLEHEADER, title, title );
		fprintf(fo, '<table WIDTH="100%%" NOSAVE>' );
		for index = 1:length(files)
	        % Processing file only
            if ~exist(fullfile(DEST, [ files{index}(1:end-1) 'html' ]))
                fprintf('Processing %s\n', files{index});
                cd(DEST); com = help2html2( files{index}, [], options{:}); cd(tmpdir);
                fprintf( fo, '%s', com);
                if strcmp(mainonly,'off')
                    inputfile = which( files{index});
                    try, copyfile( inputfile, [ DEST files{index} ]); % assuming the file is in the path 
                    catch, fprintf('Cannot copy file %s\n', inputfile); end
                end
            else
                fprintf('Skipping %s\n', files{index});
                cd(DEST); 
                com = help2html2( files{index}, [], options{:}, ...
                                 'outputonly','on');
                fprintf( fo, '%s', com);
                cd(tmpdir);
            end
		end
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

% convert spaces for .html
function strout = space2html(strin, linkb, linke)
	strout = [];
	index = 1;
	while strin(index) == ' '
		strout = [ strout '&nbsp; '];
		index = index+1;
	end
	if nargin == 3
		strout = [strout linkb strin(index:end) linke];
	else
		strout = [strout strin(index:end)];
	end
				   
