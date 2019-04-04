% help2html() - Convert a Matlab m-file help-message header into an .html help file 
%
% Usage:
%  >> linktext = help2html( filein, fileout, 'key1', val1, 'key2', val2 ...);
%
% Inputs:
%   filein       - input filename  (with .m extension)
%   fileout      - output filename (if empty, generated automatically)
%
% Optional inputs:
%   'header'     - command to insert in the header (i.e. javascript 
%                   declaration or meta-tags). {default: none}.
%   'footer'     - command to insert at the end of the .html file (e.g., 
%                   back button). {default: none}. 
%   'refcall'    - syntax to call references. {default: '%s.html'} For
%                  javascript function uses 'javascript:funcname(''%s.js'')'
%   'font'       - font name
%   'background' - background tag (i.e. '<body BACKGROUND="img.jpg">'). 
%                   {default: none}.
%   'outputlink' - html command text to link to the generated .html page. 
%                  Must contain two '%s' symbols to the function title 
%                  and to the function link.
%                  Ex: ... href="%s.html" ... {default: standard .html href}. 
%   'outputtext' - Text used in the outputlink. {default: the function
%                  name}
%   'outputonly' - ['on'|'off'] Only generate the linktext {default: 'off'}
%
% Output:
%   fileout      - .html file written to disk
%   linktext      - html-text link to the output file 
%
% M-file format:
%   The following lines describe the header format of an m-file function 
%   to allow .html help file generation. Characters '-' and ':' are used
%   explicitly by the function for parsing.
%                                               
%%    function_name() - description line 1              
%%                      description line 2       
%%                      etc.                     
%%                                               
%%    Title1:                                     
%%      descriptor1  - text line 1                 
%%                     text line 2                 
%%      descriptor2  - [type] text line 1                 
%%      "descriptor 3" - text line 1 (see notes)                 
%%                       etc.                        
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
% Notes: 1) The text lines below Title2 are considered as is (i.e.,
%           preserving Matlab carriage returns) unless there is a 
%           Matlab continuation cue ('...'). In this case, lines are 
%           contcatenated. As below 'Title1', all text lines following
%           each  descriptor (i.e., single_word followed by '-' or '='
%           or multiple quoted words followed by a '-' or '=') 
%           are concatenated. 
%        2) The pattern 'function()' is detected and is printed in bold 
%           if it is the first function descriptor. Otherwise,
%           it is used as a web link to the .html function file 
%           'function.html' if this exists.
%        3) If a 'function.jpg' image file (with same 'function' name) exists, 
%           the image is inserted into the function .html file following 
%           the function description. If the .jpg file is absent, the function
%           checks for the presence of a .gif file.
%        4) Lines beginning by '%%' are not interpreted and will be printed as is.
%        5) if [type] is present in a "descriptor2  - [type] text line 1"
%           type is bolded.

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

function [linktext,allvars,alltext] = help2html( filename, htmlfile, varargin)

if nargin < 1
	help help2html;
	return;
end
if nargin <3
   g = [];
else
   g = struct( varargin{:});;
end;	 	

try, g.font; 			catch, g.font		= 'Helvetica'; 	end

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

try, g.outputonly; 		catch, g.outputonly	= 'off'; end
try, g.background; 		catch, g.background	= ''; 	end
try, g.header; 			catch, g.header		= ''; 	end; 
try, g.footer; 			catch, g.footer		= ''; 	end; 
try, g.refcall; 		catch, g.refcall	= '%s.html'; 	end; 
try, g.outputlink; 		catch, g.outputlink = [ g.normrow g.normcol1 '<FONT FACE="' g.font '"><A HREF="%s.html">%s.html</A></td>' g.normcol2 '%s</FONT></td></tr>' ]; end; 

g.footer      = [ '<FONT FACE="'  g.font '">' g.footer '</FONT>' ];

if nargin < 1 
	help help2html;
	return;
end;	

% output file
% -----------
if nargin < 2 || isempty(htmlfile);
	indexdot = findstr( filename, '.');
	if isempty(indexdot), indexdot = length(filename)+1; end
	htmlfile = [ filename(1:indexdot(end)-1) '.html' ];
else
	indexdot = findstr( filename, '.');
end

% open files
% ---------- 
fid = fopen( filename, 'r');
if fid == -1
	error('Input file not found');
end
if ~strcmp(g.outputonly, 'on')
	fo = fopen(htmlfile, 'w');
	if fo == -1
		error('Cannot open output file');
	end
	% write header
	% ------------
	fprintf(fo, '<HTML><HEAD>%s</HEAD><BODY>\n%s\n<table WIDTH="100%%" NOSAVE>\n', g.header, g.background);
end

	
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
	  		'inputs:' 'outputs:' 'output:' 'example:' 'examples:' 'see also:' }, newtitle = 1;
		 end
		 if (i2d == length(str)) && (str(1) ~= '%'), newtitle = 1; end;	
   	  end
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
         [tok1 strrm] = mystrtok( str );
         [tok2 strrm] = strtok( strrm );

         if ~isempty(tok2) && ( tok2 == '-' || tok2 == '=') % new variable 
            newvar = 1;
            oldvarname = varname;
            oldvartext = vartext;
			if ~maindescription
				varname = formatstr( tok1, g.refcall);
            else 
				varname = tok1;
			end
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
            end
         end;	 
         newtitle = 0;		
      end
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
           		fprintf( fo, [ g.normcol2 g.vartext '</td></tr>\n' ], finalformat(oldvartext));
       	 	else
       			if ~isempty(oldvartext)
           			fprintf( fo, [ g.normcol2 g.text '</td></tr>\n' ], finalformat(oldvartext));	
       			end
         	end;	 	
         	newvar = 1;
         	oldvarname = varname;
         	oldvartext = vartext;
         end;	
      end; 

      % test if new input for an array
      % ------------------------------
      if newvar || newtitle
         if maindescription
            if ~isempty(oldvartext) % FUNCTION TITLE
               maintext = oldvartext;

			   % generate the output command
			   % ---------------------------
			   try, g.outputtext; 		catch, g.outputtext	= ''; 	end; 
			   if isempty(g.outputtext),  g.outputtext	=  filename(1:indexdot(end)-1); end
			   linktext = sprintf( g.outputlink, g.outputtext,  filename(1:indexdot(end)-1), maintext ); 
			   if strcmp(g.outputonly, 'on')
                   fclose(fid);
				   return;
			   end
			   
               maindescription = 0;
               functioname = oldvarname( 1:findstr( oldvarname, '()' )-1);
               fprintf( fo, [g.normrow g.normcol1 g.functionname '</td>\n'],upper(functioname));
               fprintf( fo, [g.normcol2 g.description '</td></tr>\n'], [ upper(oldvartext(1)) oldvartext(2:end) ]);

               % INSERT IMAGE IF PRESENT
               imagename = [];
               imagename1 = [ htmlfile( 1:findstr( htmlfile, functioname )-1) functioname '.jpg' ];
               imagename2 = [ htmlfile( 1:findstr( htmlfile, functioname )-1) functioname '.gif' ];
			   if exist( imagename2 ), imagename = imagename2; end; 
			   if exist( imagename1 ), imagename = imagename1; end; 
               if ~isempty(imagename) % do not make link if the file does not exist 
                   imageinfo = imfinfo(imagename);
                   if imageinfo.Width < 600
                       fprintf(fo, [ g.normrow g.doublecol ...
                                     '<CENTER><BR><A HREF="' imagename '" target="_blank"><img SRC=' imagename ...
                                     '></A></CENTER></td></tr>' ]);
                   else
                       fprintf(fo, [ g.normrow g.doublecol ...
                                     '<CENTER><BR><A HREF="' imagename '" target="_blank"><img SRC=' imagename ...
                                     ' width=600></A></CENTER></td></tr>' ]);
                   end
                   
               end
            end;             
   		elseif ~isempty(oldvarname)
			allvars{indexout} = oldvarname;
			alltext{indexout} = oldvartext;
			indexout = indexout + 1;
       		fprintf( fo, [ '</tr>' g.normrow g.normcol1 g.var '</td>\n' ], oldvarname);
       		fprintf( fo, [ g.normcol2 g.vartext '</td></tr>\n' ], finalformat(oldvartext));
   	 	else
   			if ~isempty(oldvartext)
       			fprintf( fo, [ g.normcol2 g.text '</td></tr>\n' ], finalformat(oldvartext));	
   			end
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
      end
   else
      str = fgets( fid );
   end
end
fprintf( fo, [ '</table>\n<BR>' g.seefile '%s</BODY></HTML>'], lower(filename), filename, g.footer);
fclose( fid );
fclose( fo  );

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
		end
return;	

% final formating
function str = finalformat(str); % bold text in bracket if just in the beginning
    tmploc = sort(union(find(str == '['), find(str == ']')));
    if ~isempty(tmploc) && str(1) == '['
        if mod(length(tmploc),2) ~= 0, str, error('Opening but no closing bracket'); end
        tmploc = tmploc(1:2);
        str = [ str(1:tmploc(1)) '<b>' str(tmploc(1)+1:tmploc(2)-1) '</b>' str(tmploc(2):end) ];
    end

    %tmploc = find(str == '"');
    %if ~isempty(tmploc)
    %    if mod(length(tmploc),2) ~= 0, str, error('Opening but no closing parenthesis'); end
    %    for index = length(tmploc):-2:1
    %        str = [ str(1:tmploc(index-1)-1) '<b>' str(tmploc(index-1)+1:tmploc(index)-1) '</b>' str(tmploc(index)+1:end) ];
    %    end
    %end
    
function tokout = functionformat( tokin, refcall );
	tokout = tokin;	% default
	[test, realtokin, tail, beg] = testfunc1( tokin );
	if ~test,  [test, realtokin, tail] = testfunc2( tokin ); end
	if test
		i1 = findstr( refcall, '%s');
		i2 = findstr( refcall(i1(1):end), '''');
		if isempty(i2) i2 = length( refcall(i1(1):end) )+1; end
		filename  = [ realtokin refcall(i1+2:i1+i2-2)]; % concatenate filename and extension
		%disp(filename)
		if exist( filename ) % do not make link if the file does not exist 
			tokout =  sprintf( [ beg '<A HREF="' refcall '">%s</A>' tail ' ' ], realtokin, realtokin );
		end
	end;		
return;

function [test, realtokin, tail, beg] = testfunc1( tokin ) % test if is string is 'function()[,]'  
	test = 0; realtokin = ''; tail = ''; beg = '';
	if ~isempty( findstr( tokin, '()' ) )
		if length(tokin)<3, return; end
		realtokin = tokin( 1:findstr( tokin, '()' )-1);
		if length(realtokin) < (length(tokin)-2) tail = tokin(end); else tail = []; end
		test = 1;
		if realtokin(1) == '(', realtokin = realtokin(2:end); beg = '('; end
		if realtokin(1) == ',', realtokin = realtokin(2:end); beg = ','; end
	end
return;

function [test, realtokin, tail] = testfunc2( tokin ) % test if is string is 'FUNCTION[,]'  
	test = 0; realtokin = ''; tail = '';
	if all( upper(tokin) == tokin)
		if tokin(end) == ',' 
			realtokin = tokin(1:end-1);
			tail = ',';
		else
			realtokin = tokin;
		end
		testokin = realtokin;
		testokin(findstr(testokin, '_')) = 'A';
		testokin(findstr(testokin, '2')) = 'A';
		if all(double(testokin) > 64) && all(double(testokin) < 91)
			test = 1;
		end;				
		realtokin = lower(realtokin);
	end
return;	

function [tok, str] = mystrtok(str)
    
    [tok str] = strtok(str);
    if tok(1) == '"'
        while tok(end) ~= '"'
            [tok2 str] = strtok(str);
            if isempty(tok2), tok, error('can not find closing quote ''"'' in previous text'); end
            tok = [tok ' ' tok2];
        end
    end
    
