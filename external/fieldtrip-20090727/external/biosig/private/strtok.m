%% Copyright (C) 2000 Paul Kienzle
%% Copyright (C) 2008 Alois Schloegl
%% $Id: strtok.m,v 1.1 2009-07-07 02:23:50 arno Exp $
%% This function is part of BioSig http://biosig.sf.net 
%% Originally, it was part of Octave. It was modified for the use with FreeMat
%%
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%

%% usage: [tok, rem] = strtok(str, delim)
%% 
%% Find all characters up to but not including the first character which
%% is in the string delim.  If rem is requested, it contains the
%% remainder of the string, starting at the first deliminator. Leading
%% delimiters are ignored.  If delim is not specified, space is assumed.

%% TODO: check what to do for a null delimiter
function [tok, rem] = strtok(str, delim)

  if nargin<1 || nargin > 2
    usage('[tok, rem] = strtok(str, delim)');
  end
  if nargin < 2 || isempty(delim), delim = " "; end

  if isempty(str)
    tok = '';rem = '';
  elseif length(delim) > 3
    start = 1;
    len = length(str);
    while start<=len
      if all(str(start) ~= delim), break; end
      start=start+1;
    end
    stop = start;
    while stop<=len
      if any(str(stop) == delim), break; end
      stop=stop+1;
    end
    tok = str(start:stop-1);
    rem = str(stop:len);
  else
    if length(delim)==1
      idx = find(str == delim);
    elseif length(delim)==2
      idx = find(str == delim(1) | str==delim(2));
    else
      idx = find(str == delim(1) | str==delim(2) | str==delim(3));
    end
    if isempty(idx)
      tok = str;
      rem = "";
    else
      skip = find(idx(:)' ~= 1:length(idx)); %% find first non-leading delimiter
      if isempty(skip)
      	tok = str(idx(length(idx))+1:length(str));
      	rem = '';
      else
      	tok = str(skip(1):idx(skip(1))-1);
      	rem = str(idx(skip(1)):length(str));
      end
    end
  end

end

%!demo
%! strtok("this is the life")
%! % split at the first space, returning "this"

%!demo
%! s = "14*27+31"
%! while 1
%!   [t,s] = strtok(s, "+-*/");
%!   printf("<%s>", t);
%!   if isempty(s), break; end
%!   printf("<%s>", s(1));
%! end
%! printf("\n");
%! % ----------------------------------------------------
%! % Demonstrates processing of an entire string split on
%! % a variety of delimiters. Tokens and delimiters are 
%! % printed one after another in angle brackets.  The
%! % string is:

%!# test the tokens for all cases
%!assert(strtok(""), "");             # no string
%!assert(strtok("this"), "this");     # no delimiter in string
%!assert(strtok("this "), "this");    # delimiter at end
%!assert(strtok("this is"), "this");  # delimiter in middle
%!assert(strtok(" this"), "this");    # delimiter at start
%!assert(strtok(" this "), "this");   # delimiter at start and end
%!assert(strtok(" "), ""(1:0));            # delimiter only

%!# test the remainder for all cases
%!test [t,r] = strtok(""); assert(r, "");
%!test [t,r] = strtok("this"); assert(r, "");
%!test [t,r] = strtok("this "); assert(r, " ");
%!test [t,r] = strtok("this is"); assert(r, " is");
%!test [t,r] = strtok(" this"); assert(r, "");
%!test [t,r] = strtok(" this "); assert(r, " ");
%!test [t,r] = strtok(" "); assert(r, "");

%!# simple check with 2 and 3 delimeters
%!assert(strtok("this is", "i "), "th");
%!assert(strtok("this is", "ij "), "th");

%!# test all cases for 4 delimiters since a different 
%!# algorithm is used when more than 3 delimiters
%!assert(strtok("","jkl "), "");
%!assert(strtok("this","jkl "), "this");
%!assert(strtok("this ","jkl "), "this");
%!assert(strtok("this is","jkl "), "this");
%!assert(strtok(" this","jkl "), "this");
%!assert(strtok(" this ","jkl "), "this");
%!assert(strtok(" ","jkl "), ""(1:0));

%!# test 'bad' string orientations
%!assert(strtok(" this "'), "this"');   # delimiter at start and end
%!assert(strtok(" this "',"jkl "), "this"');
