% help2html() - convert a Matlab m-file help-message header 
%               into an .html help file 
% Usage:
%  >> linktext = help2html( filein, fileout, 'key1', val1, 'key2', val2 ...);
%
% Inputs:
%   filein     - input filename  (with .m extension)
%   fileout    - output filename (if empty, generated automatically)
%
% Optional inputs:
%   'header'   - command to insert in the header (i.e. javascript 
%                   declaration or meta-tags). {Default none. 
%   'footer'   - command to insert at the end of the HTML file (i.e. back
%                   button. {Default none}. 
%   'refcall'  - syntax format to call references. Default '%s.html'. For
%                   javascript function uses 'javascript:funcname(''%s.js'')'
%   'font'     - font name
%   'background'- background tag (i.e. '<body BACKGROUND="img.jpg">'). 
%                   {Default none}.
%   'outputlink'- command to call the generated HTML page. Default is
%                   standard HTML href.    
%
% Output:
%   linktext   - html text link to the output file
%
% M-file format:
%   The following lines describe the header format of an m-file function 
%   to allow html help file generation. Characters '-' and ':' are used
%   by the function for parsing.
%                                               
%%    function_name() - description line 1              
%%                      description line 2       
%%                      etc.                     
%%                                               
%%    Title1:                                     
%%      variable1  - text line 1                 
%%                   text line 2                 
%%      variable2  - text line 1                 
%%                   etc.                        
%%                                              
%%    Title2:                                    
%%      text line 1 [...](see notes)             
%%      text line 2                              
%%                                               
%%    See also:                                  
%%     function1(), function2()                     
%    
% Author:  Arnaud Delorme, Salk Institute 2001
%
% Notes: 1) In Title2, the text lines are considered as is (e.g., 
%           preserving Matlab carriage returns) unless there is a 
%           matlab continuation cue ('...'). In this case, lines are 
%           contcatenated. For Title1 with variable's name, all text lines 
%           are concatenated by default.
%        2) The pattern 'function()[,]' is detected 
%           It is printed in bold if it's the function descriptor
%           It is used as a web link to another function, but only if the html 
%           file 'function.html' exists.
%        3) If a '.jpg' file with the same name as the function exists, it is
%           inserted into the HTLM file after the function description.
%        4) Lines beginning by '%%' are not interpreted

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
% Revision 1.2  2002/04/22 22:58:30  arno
% adding extra output parameters
%
% Revision 1.1  2002/04/05 17:39:45  jorn
% Initial revision
%

function [linktext,allvars,alltext] = help2html( filename, htmlfile, varargin)

if nargin < 1
	help help2html;
	return;
end;
if nargin <3
   g = [];
else
   g = struct( varargin{:});;
end;	 	

try, g.font; 			catch, g.font		= 'Helvetica'; 	end;

g.functionname = [ '<FONT FACE="' g.font '"><FONT SIZE =+2><B>%s</B></FONT></FONT>' ];
g.description  = [ '<FONT FACE="' g.font '">%s</FONT>' ];
g.title        = [ '<FONT FACE="' g.font '"><FONT SIZE =+1><B>%s</B></FONT></FONT>' ];
g.var          = [ '<DIV ALIGN=RIGHT><FONT FACE="' g.font '"><I>%s&nbsp;&nbsp;&nbsp;</I></FONT></DIV>' ];
g.vartext      = [ '<FONT FACE="' g.font '">%s</FONT>' ];
g.text         = [ '<FONT FACE="' g.font '">%s</FONT>' ];

g.normrow      = '<tr VALIGN=TOP NOSAVE>';
g.normcol1     = '<td VALIGN=TOP NOSAVE>';
g.normcol2     = '<td VALIGN=BOTTOM NOSAVE>';
g.doublecol    = '<td ALIGN=CENTER COLSPAN="2" NOSAVE>';
g.seefile      = [ '<FONT FACE="'  g.font '">See the matlab file <A HREF="%s" target="_blank">%s</A> (may require other functions)</FONT><BR><BR>' ]; 

try, g.background; 		catch, g.background	= ''; 	end;
try, g.header; 			catch, g.header		= ''; 	end; 
try, g.footer; 			catch, g.footer		= ''; 	end; 
try, g.refcall; 		catch, g.refcall	= '%s.html'; 	end; 
try, g.outputlink; 		catch, g.outputlink = [ g.normrow g.normcol1 '<FONT FACE="' g.font '"><A HREF="%s.html">%s.html</A></td>' g.normcol2 '%s</FONT></td></tr>' ]; end; 

g.footer      = [ '<FONT FACE="'  g.font '">' g.footer '</FONT>' ];

if nargin < 1 
	help help2html;
	return;
end;	

% input file
% ---------- 
fid = fopen( filename, 'r');
if fid == -1
	error('File not found');
end;

% output file
% -----------
if nargin < 2 | isempty(htmlfile);
	index = findstr( filename, '.');
	if isempty(index), index = length(filename)+1; end;
	htmlfile = [ filename(1:index(end)-1) '.html' ];
else
	index = findstr( filename, '.');
end;
fo = fopen(htmlfile, 'w');

% write header
% ------------
fprintf(fo, '<HTML><HEAD>%s</HEAD><BODY>\n%s\n<table WIDTH="100%%" NOSAVE>\n', g.header, g.background);
	
cont = 1;
% scan file
% -----------
str = fgets( fid );
% if first line is the name of the function, reload
if str(1) ~= '%', str = fgets( fid ); end; 

% state variables
% -----------
maindescription 	 = 1;
varname = [];
oldvarname = [];
vartext = []; 	
newvar  	 = 0;
allvars = {};
alltext = {};
indexout = 1;

while (str(1) == '%')
   str = deblank(str(2:end-1));
  
   % --- DECODING INPUT
   newvar  	 = 0;
   str = deblank2(str); 
   
   if ~isempty(str)
      
	  % find a title
	  % ------------
	  i2d = findstr(str, ':');
	  newtitle = 0;
	  if ~isempty(i2d)
	  	 i2d = i2d(1);
	  	 switch lower(str(1:i2d))
	  		case { 'usage:' 'authors:' 'author:' 'notes:' 'note:' 'input:' ...
	  		'inputs:' 'outputs:' 'output' 'example:' 'examples:' 'see also:' }, newtitle = 1;
		 end;
		 if (i2d == length(str)) & (str(1) ~= '%'), newtitle = 1; end;	
   	  end;
      if newtitle
  			tilehtml = str(1:i2d); 
  			newtitle = 1;
	        oldvarname = varname;
			oldvartext = vartext;
        	if i2d < length(str)
         			vartext = formatstr(str(i2d+1:end), g.refcall);
         	else	vartext = [];
         	end;	
         	varname = [];
      else
	  % not a title
	  % ------------
         % scan lines
         [tok1 strrm] = strtok( str );
         [tok2 strrm] = strtok( strrm );

         if ~isempty(tok2) & ( tok2 == '-' | tok2 == '=') % new variable 
            newvar = 1;
            oldvarname = varname;
            oldvartext = vartext;
            varname = formatstr( tok1, g.refcall);
            strrm = deblank(strrm);            % remove tail blanks
            strrm = deblank(strrm(end:-1:1));	% remove initial blanks 
           	strrm = formatstr( strrm(end:-1:1), g.refcall);
            vartext = strrm;
         else
            % continue current text
            str = deblank(str);            % remove tail blanks
            str = deblank(str(end:-1:1));	% remove initial blanks
            str = formatstr( str(end:-1:1), g.refcall );
            if isempty(vartext)
               vartext = str;
            else	
               if ~isempty(varname) 
               	    vartext = [ vartext ' ' str]; % espace if in array
               else 
               		if all(vartext(	end-2:end) == '.')
               			vartext = [ deblank2(vartext(1:end-3)) ' ' str]; % espace if '...'
               		else
                    	vartext = [ vartext '<BR>' str];    % CR otherwise
                    end;	
               end;		
            end;
         end;	 
         newtitle = 0;		
      end;
	  % --- END OF DECODING 	
      
   	  str = fgets( fid );

      % test if last entry
      % ------------------
      if str(1) ~= '%'
  	 	 if ~newtitle
  	 	 	if ~isempty(oldvarname)
				allvars{indexout} = oldvarname;
				alltext{indexout} = oldvartext;
				indexout = indexout + 1;
           		fprintf( fo, [ '</tr>' g.normrow g.normcol1 g.var '</td>\n' ], oldvarname);
           		fprintf( fo, [ g.normcol2 g.vartext '</td></tr>\n' ], oldvartext);
       	 	else
       			if ~isempty(oldvartext)
           			fprintf( fo, [ g.normcol2 g.tabtext '</td></tr>\n' ], oldvartext);	
       			end;
         	end;	 	
         	newvar = 1;
         	oldvarname = varname;
         	oldvartext = vartext;
         end;	
      end; 

      % test if new input for an array
      % ------------------------------
      if newvar | newtitle
         if maindescription
            if ~isempty(oldvartext) % FUNCTION TITLE
               maintext = oldvartext;
               maindescription = 0;
               functioname = oldvarname( 1:findstr( oldvarname, '()' )-1);
               fprintf( fo, [g.normrow g.normcol1 g.functionname '</td>\n'],upper(functioname));
               fprintf( fo, [g.normcol2 g.description '</td></tr>\n'], [ upper(oldvartext(1)) oldvartext(2:end) ]);

               % INSERT IMAGE IF PRESENT
               imagename = [ htmlfile( 1:findstr( htmlfile, functioname )-1) functioname '.jpg' ];
			   if exist( imagename ) % do not make link if the file does not exist 
					fprintf(fo, [ g.normrow g.doublecol ...
								'<CENTER><BR><A HREF="' imagename '" target="_blank"><img SRC=' imagename ...
								' height=150 width=200></A></CENTER></td></tr>' ]);
		       end;
            end;             
   		elseif ~isempty(oldvarname)
			allvars{indexout} = oldvarname;
			alltext{indexout} = oldvartext;
			indexout = indexout + 1;
       		fprintf( fo, [ '</tr>' g.normrow g.normcol1 g.var '</td>\n' ], oldvarname);
       		fprintf( fo, [ g.normcol2 g.vartext '</td></tr>\n' ], oldvartext);
   	 	else
   			if ~isempty(oldvartext)
       			fprintf( fo, [ g.normcol2 g.text '</td></tr>\n' ], oldvartext);	
   			end;
         end;      
      end;	
      
      % print title
      % -----------
      if newtitle
		 fprintf( fo, [ g.normrow g.doublecol '<BR></td></tr>' g.normrow g.normcol1 g.title '</td>\n' ], tilehtml);
      	 if str(1) ~= '%' % last input
         	fprintf( fo, [ g.normcol2 g.text '</td></tr>\n' ], vartext);
         end;		
         oldvarname = [];
         oldvartext = [];
      end;
   else
      str = fgets( fid );
   end;
end;
fprintf( fo, [ '</table>\n<BR>' g.seefile '%s</BODY></HTML>'], lower(filename), filename, g.footer);
fclose( fid );
fclose( fo  );

% generate the output command
% ---------------------------
linktext = sprintf( g.outputlink,  filename(1:index(end)-1),  filename(1:index(end)-1), maintext ); 
return;

% -----------------
% sub-functions
% -----------------
function str = deblank2( str ); % deblank two ways
   str = deblank(str(end:-1:1));    % remove initial blanks
   str = deblank(str(end:-1:1));	% remove tail blanks 
return;

function strout = formatstr( str, refcall );
		[tok1 strrm] = strtok( str );
		strout = [];
		while ~isempty(tok1)
			tokout = functionformat( tok1, refcall );
			if isempty( strout)
				strout = tokout; 	
			else
				strout = [strout ' ' tokout ]; 	
			end;	
			[tok1 strrm] = strtok( strrm );
		end;
return;	
 
function tokout = functionformat( tokin, refcall );
	tokout = tokin;	% default
	[test, realtokin, tail] = testfunc1( tokin );
	if ~test,  [test, realtokin, tail] = testfunc2( tokin ); end;
	if test
		i1 = findstr( refcall, '%s');
		i2 = findstr( refcall(i1(1):end), '''');
		if isempty(i2) i2 = length( refcall(i1(1):end) )+1; end;
		filename  = [ realtokin refcall(i1+2:i1+i2-2)]; % concatenate filename and extension
		disp(filename)
		if exist( filename ) % do not make link if the file does not exist 
			tokout =  sprintf( [ '<A HREF="' refcall '">%s</A>' tail ' ' ], realtokin, realtokin );
		end;
	end;		
return;

function [test, realtokin, tail] = testfunc1( tokin ) % test if is string is 'function()[,]'  
	test = 0; realtokin = ''; tail = '';
	if ~isempty( findstr( tokin, '()' ) )
		realtokin = tokin( 1:findstr( tokin, '()' )-1);
		if length(realtokin) < (length(tokin)-2) tail = tokin(end); else tail = []; end;
		test = 1;
		if realtokin(1) == '(', realtokin = realtokin(2:end); end;
		if realtokin(1) == ',', realtokin = realtokin(2:end); end;
	end;
return;

function [test, realtokin, tail] = testfunc2( tokin ) % test if is string is 'FUNCTION[,]'  
	test = 0; realtokin = ''; tail = '';
	if all( upper(tokin) == tokin)
		if tokin(end) == ',' 
			realtokin = tokin(1:end-1);
			tail = ',';
		else
			realtokin = tokin;
		end;
		testokin = realtokin;
		testokin(findstr(testokin, '_')) = 'A';
		testokin(findstr(testokin, '2')) = 'A';
		if all(double(testokin) > 64) & all(double(testokin) < 91)
			test = 1;
		end;				
		realtokin = lower(realtokin);
	end;
return;	
