% gethelpvar() - convert a Matlab m-file help-message header 
%               into out variables 
% Usage:
%  >> [vartext, varnames] = gethelpvar( filein, varlist);
%
% Inputs:
%   filein     - input filename  (with .m extension)
%   varlist    - optional list of variable to return. If absent
%                retrun all variables.
%
% Outputs:
%   vartext     - variable text
%   varnames    - name of the varialbe returned;
%
% Notes: see help2html() for file format
%
% Example: >> gethelpvar( 'gethelpvar.m', 'filein')
%
% Author:  Arnaud Delorme, Salk Institute 20 April 2002
%
% See also: help2html()

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

function [alltext, allvars] = gethelpvar( filename, varlist )

if nargin < 1
	help gethelpvar;
	return;
end

% input file
% ---------- 
fid = fopen( filename, 'r');
if fid == -1
	error('File not found');
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
	  		'inputs:' 'outputs:' 'output' 'example:' 'examples:' 'see also:' }, newtitle = 1;
		 end
		 if (i2d == length(str)) && (str(1) ~= '%'), newtitle = 1; end;	
   	  end
      if newtitle
  			tilehtml = str(1:i2d); 
  			newtitle = 1;
	        oldvarname = varname;
        		oldvartext = vartext;
        	if i2d < length(str)
         			%vartext = formatstr(str(i2d+1:end), g.refcall);
         			vartext = str(i2d+1:end);
         	else	vartext = [];
         	end;	
         	varname = [];
      else
	  % not a title
	  % ------------
         % scan lines
         [tok1 strrm] = strtok( str );
         [tok2 strrm] = strtok( strrm );

         if ~isempty(tok2) && ( isequal(tok2,'-') || isequal(tok2,'=')) % new variable 
            newvar = 1;
            oldvarname = varname;
            oldvartext = vartext;
            varname = tok1;
            strrm = deblank(strrm);            % remove tail blanks
            strrm = deblank(strrm(end:-1:1));	% remove initial blanks 
           	%strrm = formatstr( strrm(end:-1:1), g.refcall);
           	strrm = strrm(end:-1:1);
            vartext = strrm;
         else
            % continue current text
            str = deblank(str);            % remove tail blanks
            str = deblank(str(end:-1:1));	% remove initial blanks
            %str = formatstr( str(end:-1:1), g.refcall );
            str = str(end:-1:1);
            if isempty(vartext)
               vartext = str;
            else	
               if ~isempty(varname) 
               	    vartext = [ vartext 10 str]; % espace if in array
               else 
               		if length(vartext)>3 && all(vartext(	end-2:end) == '.')
               			vartext = [ deblank2(vartext(1:end-3)) 10 str]; % espace if '...'
               		else
                    	vartext = [ vartext 10 str];    % CR otherwise
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
           		%fprintf( fo, [ '</tr>' g.normrow g.normcol1 g.var '</td>\n' ], oldvarname);
           		%fprintf( fo, [ g.normcol2 g.vartext '</td></tr>\n' ], oldvartext);
       	 	else
       			if ~isempty(oldvartext)
           			%fprintf( fo, [ g.normcol2 g.tabtext '</td></tr>\n' ], oldvartext);	
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
               maindescription = 0;
               functioname = oldvarname( 1:findstr( oldvarname, '()' )-1);
               %fprintf( fo, [g.normrow g.normcol1 g.functionname '</td>\n'],upper(functioname));
               %fprintf( fo, [g.normcol2 g.description '</td></tr>\n'], [ upper(oldvartext(1)) oldvartext(2:end) ]);

               % INSERT IMAGE IF PRESENT
               %imagename = [ htmlfile( 1:findstr( htmlfile, functioname )-1) functioname '.jpg' ];
			   %if exist( imagename ) % do not make link if the file does not exist 
					%fprintf(fo, [ g.normrow g.doublecol ...
					%			'<CENTER><BR><A HREF="' imagename '" target="_blank"><img SRC=' imagename ...
					%			' height=150 width=200></A></CENTER></td></tr>' ]);
		       %end
            end;             
   		elseif ~isempty(oldvarname)
			allvars{indexout} = oldvarname;
			alltext{indexout} = oldvartext;
			indexout = indexout + 1;
       		%fprintf( fo, [ '</tr>' g.normrow g.normcol1 g.var '</td>\n' ], oldvarname);
       		%fprintf( fo, [ g.normcol2 g.vartext '</td></tr>\n' ], oldvartext);
   	 	else
   			if ~isempty(oldvartext)
       			%fprintf( fo, [ g.normcol2 g.text '</td></tr>\n' ], oldvartext);	
   			end
         end;      
      end;	
      
      % print title
      % -----------
      if newtitle
		 %fprintf( fo, [ g.normrow g.doublecol '<BR></td></tr>' g.normrow g.normcol1 g.title '</td>\n' ], tilehtml);
      	 if str(1) ~= '%' % last input
         	%fprintf( fo, [ g.normcol2 g.text '</td></tr>\n' ], vartext);
         end;		
         oldvarname = [];
         oldvartext = [];
      end
   else
      str = fgets( fid );
   end
end
fclose( fid );

% remove quotes of variables
% --------------------------
for index = 1:length(allvars)
	if allvars{index}(1) == '''', allvars{index} = eval( allvars{index} ); end
end

if exist('varlist') == 1
	if ~iscell(varlist), varlist = { varlist }; end
        newtxt = mat2cell(zeros(length(varlist), 1), length(varlist), 1); % preallocation
	for index = 1:length(varlist)
		loc = strmatch( varlist{index}, allvars);
		if ~isempty(loc)
			newtxt{index} = alltext{loc(1)};
		else
			disp([ 'warning: variable ''' varlist{index} ''' not found']);
			newtxt{index} = '';
		end
	end
	alltext = newtxt;
end
	
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
 
function tokout = functionformat( tokin, refcall );
	tokout = tokin;	% default
	[test, realtokin, tail] = testfunc1( tokin );
	if ~test,  [test, realtokin, tail] = testfunc2( tokin ); end
	if test
		i1 = findstr( refcall, '%s');
		i2 = findstr( refcall(i1(1):end), '''');
		if isempty(i2) i2 = length( refcall(i1(1):end) )+1; end
		filename  = [ realtokin refcall(i1+2:i1+i2-2)];
		if exist( filename ) % do not make link if the file does not exist 
			tokout =  sprintf( [ '<A HREF="' refcall '">%s</A>' tail ' ' ], realtokin, realtokin );
		end
	end;		
return;

function [test, realtokin, tail] = testfunc1( tokin ) % test if is string is 'function()[,]'  
	test = 0; realtokin = ''; tail = '';
	if ~isempty( findstr( tokin, '()' ) )
		realtokin = tokin( 1:findstr( tokin, '()' )-1);
		if length(realtokin) < (length(tokin)-2) tail = tokin(end); else tail = []; end
		test = 1;
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
