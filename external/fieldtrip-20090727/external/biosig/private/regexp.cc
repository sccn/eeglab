/*

Copyright (C) 2005 David Bateman
Copyright (C) 2002-2005 Paul Kienzle

Octave is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2, or (at your option) any
later version.

Octave is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; see the file COPYING.  If not, write to the
Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
Boston, MA 02110-1301, USA.

*/


#include <octave/config.h>
#undef HAVE_PCRE
#define HAVE_PCRE 1
typedef int octave_idx_type;




#ifdef HAVE_CONFIG_H
# include <config.h>
#else
# include <octave/config.h>
#endif

#include <algorithm>
#include <sstream>

#include "defun-dld.h"
#include "error.h"
#include "gripes.h"
#include "oct-obj.h"
#include "utils.h"

#include "Cell.h"
#include "oct-map.h"
#include "str-vec.h"
#include "quit.h"

#ifdef HAVE_PCRE
#include <pcre.h>
#else
#ifdef HAVE_REGEX
#ifdef __MINGW32__
#define __restrict
#endif
#include <regex.h>
#endif
#endif

// FIXME rewrite so that regexp builds up a list of matches in a linked
// structure then converts this structure to a contiguous array.  Resizing
// the contiguous arrays for each match takes quadratic time.  Probably should
// have in the order of 2k blocks in linked structure.  Don't bother collecting
// and composing return values the user does not want.

// FIXME refactor and incorporate regexprep so that it works directly off the
// link structure rather than constructing a cell array of matrices for it to
// work from.

static octave_value_list
octregexp (const octave_value_list &args, int nargout, const std::string &nm,
	   bool case_insensitive)
{
  octave_value_list retval;
#if defined (HAVE_REGEX) || defined (HAVE_PCRE) 
  int nargin = args.length();
  int nopts = nargin - 2;
  bool once = false;
  bool lineanchors = false;
  bool dotexceptnewline = false;
  bool freespacing = false;

  if (nargin < 2)
    {
      print_usage(nm);
      return retval;
    }

  std::string buffer = args(0).string_value ();
  if (error_state)
    {
      gripe_wrong_type_arg (nm.c_str(), args(0));
      return retval;
    }

  std::string pattern = args(1).string_value ();
  if (error_state)
    {
      gripe_wrong_type_arg (nm.c_str(), args(1));
      return retval;
    }

  for (int i = 2; i < nargin; i++)
    {
      std::string str = args(i).string_value();
      if (error_state)
	{
	  error ("%s: optional arguments must be strings", nm.c_str());
	  break;
	}
      std::transform (str.begin (), str.end (), str.begin (), tolower);
      if (str.find("once", 0) == 0)
	{
	  once = true;
	  nopts--;
	}
      else if (str.find("ignorecase", 0) == 0)
	{
	  case_insensitive = true;
	  nopts--;
	}
      else if (str.find("dotall", 0) == 0)
	{
	  dotexceptnewline = false;
	  nopts--;
	}
      else if (str.find("stringanchors", 0) == 0)
	{
	  lineanchors = false;
	  nopts--;
	}
      else if (str.find("matchcase", 0) == 0)
	{
	  case_insensitive = false;
	  nopts--;
	}
      else if (str.find("literalspacing", 0) == 0)
	{
	  freespacing = false;
	  nopts--;
	}
#if HAVE_PCRE
      // Only accept these options with pcre
      else if (str.find("dotexceptnewline", 0) == 0)
	{
	  dotexceptnewline = true;
	  nopts--;
	}
      else if (str.find("lineanchors", 0) == 0)
	{
	  lineanchors = true;
	  nopts--;
	}
      else if (str.find("freespacing", 0) == 0)
	{
	  freespacing = true;
	  nopts--;
	}
      else if (str.find("start", 0) && str.find("end", 0) &&
	       str.find("tokenextents", 0) && str.find("match", 0) &&
	       str.find("tokens", 0) && str.find("names", 0))
	error ("%s: unrecognized option", nm.c_str());
#else
      else if (str.find("names", 0) == 0 ||
	       str.find("dotexceptnewline", 0) == 0 ||
	       str.find("lineanchors", 0) == 0 ||
	       str.find("freespacing", 0) == 0)
	error ("%s: %s not implemented in this version", str.c_str(), nm.c_str());
      else if (str.find("start", 0) && str.find("end", 0) &&
	       str.find("tokenextents", 0) && str.find("match", 0) &&
	       str.find("tokens", 0))
	error ("%s: unrecognized option", nm.c_str());
#endif
    }

  if (!error_state)
    {
      Octave_map nmap;
      Cell t, m, te;
      NDArray s, e;

      // named tokens "(?<name>...)" are only treated with PCRE not regex.
#if HAVE_PCRE
      // The syntax of named tokens in pcre is "(?P<name>...)" while we need
      // a syntax "(?<name>...)", so fix that here. Also an expression like
      // "(?<first>\w+)\s+(?<last>\w+)|(?<last>\w+),\s+(?<first>\w+)" should
      // be perfectly legal, while pcre does not allow the same named token
      // name on both sides of the alternative. Also fix that here by replacing
      // name tokens by dummy names, and dealing with the dummy names later.
      
      size_t pos = 0;
      size_t new_pos;
      string_vector named;
      int nnames = 0;
      int inames = 0;
      std::ostringstream buf;
      Array<int> named_idx;

      // Add mode flags
      while ((new_pos = pattern.find ("(?<",pos)) != NPOS)
	{
	  size_t tmp_pos = pattern.find_first_of ('>',new_pos);

	  if (tmp_pos == NPOS)
	    {
	      error ("syntax error in pattern");
	      break;
	    }

	  std::string tmp_name = pattern.substr(new_pos+3,tmp_pos-new_pos-3);
	  bool found = false;

	  for (int i = 0; i < nnames; i++)
	    if (named(i) == tmp_name)
	      {
		named_idx.resize(inames+1);
		named_idx(inames) = i;
		found = true;
		break;
	      }
	  if (! found)
	    {
	      named_idx.resize(inames+1);
	      named_idx(inames) = nnames;
	      named.append(tmp_name);
	      nnames++;
	    }

	  if (new_pos - pos > 0)
	    buf << pattern.substr(pos,new_pos-pos);
	  if (inames < 10)
	    buf << "(?P<n00" << inames++;
	  else if (inames < 100)
	    buf << "(?P<n0" << inames++;
	  else
	    buf << "(?P<n" << inames++;
	  pos = tmp_pos;
	}

      buf << pattern.substr(pos);

      if (error_state)
	return retval;

      // Compile expression
      pcre *re;
      const char *err;
      int erroffset;
      std::string buf_str = buf.str ();
      re = pcre_compile (buf_str.c_str (),
			 (case_insensitive ? PCRE_CASELESS : 0) |
			 (dotexceptnewline ? 0 : PCRE_DOTALL) |
			 (lineanchors ? PCRE_MULTILINE : 0) |
			 (freespacing ? PCRE_EXTENDED : 0),
			 &err, &erroffset, NULL);
    
      if (re == NULL) {
	error("%s: %s at position %d of expression", nm.c_str(), 
	      err, erroffset);
	return retval;
      }

      int subpatterns;
      int namecount;
      int nameentrysize;
      char *nametable;
      int idx = 0;
      int sz = 0;

      pcre_fullinfo(re, NULL, PCRE_INFO_CAPTURECOUNT,  &subpatterns);
      pcre_fullinfo(re, NULL, PCRE_INFO_NAMECOUNT, &namecount);
      pcre_fullinfo(re, NULL, PCRE_INFO_NAMEENTRYSIZE, &nameentrysize);
      pcre_fullinfo(re, NULL, PCRE_INFO_NAMETABLE, &nametable);

      OCTAVE_LOCAL_BUFFER(int, ovector, (subpatterns+1)*3);
      OCTAVE_LOCAL_BUFFER(int, nidx, namecount);

      for (int i = 0; i < namecount; i++)
	{
	  // Index of subpattern in first two bytes MSB first of name.
	  // Extract index.
	  nidx[i] = (static_cast<int>(nametable[i*nameentrysize])) << 8 |
	    static_cast<int>(nametable[i*nameentrysize+1]);
	}

      Cell named_tokens(dim_vector(nnames,1));

      while(true)
	{
	  OCTAVE_QUIT;

	  int matches = pcre_exec(re, NULL, buffer.c_str(), 
				  buffer.length(), idx, 
				  (idx ? PCRE_NOTBOL : 0),
				  ovector, (subpatterns+1)*3);

	  if (matches < 0 && matches != PCRE_ERROR_NOMATCH)
	    {
	      error ("%s: internal error calling pcre_exec", nm.c_str());
	      pcre_free(re);
	      return retval;
	    }
	  else if (matches == PCRE_ERROR_NOMATCH)
	    break;
	  else if (ovector[1] <= ovector[0])
	    break;
	  else
	    {
	      int pos_match = 0;
	      Matrix mat_te(matches-1,2);
	      for (int i = 1; i < matches; i++)
		{
		  if (ovector[2*i] >= 0 && ovector[2*i+1] > 0)
		    {
		      mat_te(pos_match,0) = double (ovector[2*i]+1);
		      mat_te(pos_match++,1) = double (ovector[2*i+1]);
		    }
		}
	      mat_te.resize(pos_match,2);
#if 1
	      s.resize (dim_vector(1, sz+1));
	      s(sz) = double (ovector[0]+1);
	      e.resize (dim_vector(1, sz+1));
	      e(sz) = double (ovector[1]);
	      te.resize(dim_vector(1, sz+1));
	      te(sz) = mat_te;
#endif

	      const char **listptr;
	      int status = pcre_get_substring_list(buffer.c_str(), ovector, 
						   matches, &listptr);

	      if (status == PCRE_ERROR_NOMEMORY) {
		error("%s: cannot allocate memory in pcre_get_substring_list",
		      nm.c_str());
		pcre_free(re);
		return retval;
	      }


	      Cell cell_t (dim_vector(1,pos_match));
	      pos_match = 0;
	      for (int i = 1; i < matches; i++)
		if (ovector[2*i] >= 0 && ovector[2*i+1] > 0)
		  cell_t(pos_match++) = std::string(*(listptr+i));

#if 1
	      m.resize (dim_vector(1, sz+1));
	      m(sz) =  std::string(*listptr);
	      t.resize (dim_vector(1, sz+1));
	      t(sz) = cell_t;
#endif

	      if (namecount > 0)
		for (int i = 1; i < matches; i++)
		  {
		    if (ovector[2*i] >= 0 && ovector[2*i+1] > 0)	
		      {
			if (sz == 0)
			  {
			    named_tokens(named_idx(i-1)) = 
			      std::string(*(listptr+nidx[i-1]));
			  }
			else
			  {
			    Cell tmp = named_tokens(named_idx(i-1));
			    tmp.resize(dim_vector(1,sz+1));
			    tmp(sz) = std::string(*(listptr+nidx[i-1]));
			    named_tokens(named_idx(i-1)) = tmp;
			  }
		      }
		  }

	      pcre_free_substring_list(listptr);

	      if (once)
		break;

	      idx = ovector[1];
	      sz++;
	    }
	}

      for (int i = 0; i < nnames; i++)
	nmap.assign (named(i), named_tokens(i));

      pcre_free(re);
#else
      regex_t compiled;
      int err=regcomp(&compiled, pattern.c_str(), REG_EXTENDED | 
		      (case_insensitive ? REG_ICASE : 0));
      if (err)
	{
	  int len = regerror(err, &compiled, NULL, 0);
	  OCTAVE_LOCAL_BUFFER (char, errmsg, len);
	  regerror(err, &compiled, errmsg, len);
	  error("%s: %s in pattern (%s)", nm.c_str(), errmsg, 
		pattern.c_str());
	  regfree(&compiled);
	  return retval;
	}

      int subexpr = 1;
      int idx = 0;
      int sz = 0;
      for (unsigned int i=0; i < pattern.length(); i++)
	  subexpr += ( pattern[i] == '(' ? 1 : 0 );
      OCTAVE_LOCAL_BUFFER (regmatch_t, match, subexpr );

      while(true)
	{
	  OCTAVE_QUIT;
	  if (regexec(&compiled, buffer.c_str() + idx, subexpr, 
		      match, (idx ? REG_NOTBOL : 0)) == 0) 
	    {
	      // Count actual matches
	      int matches = 0;
	      while (matches < subexpr && match[matches].rm_so >= 0) 
		matches++;

	      s.resize (dim_vector(1, sz+1));
	      s(sz) = double (match[0].rm_so+1+idx);
	      e.resize (dim_vector(1, sz+1));
	      e(sz) = double (match[0].rm_eo+idx);
	      te.resize(dim_vector(1, sz+1));
	      Matrix mat_te(matches-1,2);
	      for (int i = 1; i < matches; i++)
		{
		  mat_te(i-1,0) = double (match[i].rm_so+1+idx);
		  mat_te(i-1,1) = double (match[i].rm_eo+idx);
		}
	      te(sz) = mat_te;

	      m.resize (dim_vector(1, sz+1));
	      m(sz) =  buffer.substr (match[0].rm_so+idx, 
					 match[0].rm_eo-match[0].rm_so);

	      t.resize (dim_vector(1, sz+1));
	      Cell cell_t (dim_vector(1,matches-1));
	      for (int i = 1; i < matches; i++)
		cell_t(i-1) = buffer.substr (match[i].rm_so+idx, 
					     match[i].rm_eo-match[i].rm_so);
	      t(sz) = cell_t;

	      idx += match[0].rm_eo;
	      sz++;

	      if (once)
		break;
	    }
	  else
	    break;
	}
      regfree(&compiled);
#endif

      retval(5) = nmap;
      retval(4) = t;
      retval(3) = m;
      retval(2) = te;
      retval(1) = e;
      retval(0) = s;

      // Alter the order of the output arguments
      if (nopts > 0)
	{
	  int n = 0;
	  octave_value_list new_retval;
	  new_retval.resize(nargout);

	  OCTAVE_LOCAL_BUFFER (int, arg_used, 6);
	  for (int i = 0; i < 6; i++)
	    arg_used[i] = false;
	  
	  for (int i = 2; i < nargin; i++)
	    {
	      int k = 0;
	      std::string str = args(i).string_value();
	      std::transform (str.begin (), str.end (), str.begin (), tolower);
	      if (str.find("once", 0) == 0
		  || str.find("stringanchors", 0) == 0
		  || str.find("lineanchors", 0) == 0
		  || str.find("matchcase", 0) == 0
		  || str.find("ignorecase", 0) == 0
		  || str.find("dotall", 0) == 0
		  || str.find("dotexceptnewline", 0) == 0
		  || str.find("literalspacing", 0) == 0
		  || str.find("freespacing", 0) == 0
	      )
		continue;
	      else if (str.find("start", 0) == 0)
		k = 0;
	      else if (str.find("end", 0) == 0)
		k = 1;
	      else if (str.find("tokenextents", 0) == 0)
		k = 2;
	      else if (str.find("match", 0) == 0)
		k = 3;
	      else if (str.find("tokens", 0) == 0)
		k = 4;
	      else if (str.find("names", 0) == 0)
		k = 5;

	      new_retval(n++) = retval(k);
	      arg_used[k] = true;

	      if (n == nargout)
		break;
	    }

	  // Fill in the rest of the arguments
	  if (n < nargout)
	    {
	      for (int i = 0; i < 6; i++)
		{
		  if (! arg_used[i])
		    new_retval(n++) = retval(i);
		}
	    }

	  retval = new_retval;
	}
    }

#else
  error ("%s: not available in this version of Octave", nm.c_str());
#endif
  return retval;
}

DEFUN_DLD (regexp, args, nargout,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {[@var{s}, @var{e}, @var{te}, @var{m}, @var{t}, @var{nm}] =} regexp (@var{str}, @var{pat})\n\
@deftypefnx {Loadable Function} {[@dots{}] =} regexp (@var{str}, @var{pat}, @var{opts}, @dots{})\n\
\n\
Regular expression string matching. Matches @var{pat} in @var{str} and\n\
returns the position and matching substrings or empty values if there are\n\
none.\n\
\n\
The matched pattern @var{pat} can include any of the standard regex\n\
operators, including:\n\
\n\
@table @code\n\
@item .\n\
Match any character\n\
@item * + ? @{@}\n\
Repetition operators, representing\n\
@table @code\n\
@item *\n\
Match zero or more times\n\
@item +\n\
Match one or more times\n\
@item ?\n\
Match zero or one times\n\
@item @{@}\n\
Match range operator, which is of the form @code{@{@var{n}@}} to match exactly\n\
@var{n} times, @code{@{@var{m},@}} to match @var{m} or more times,\n\
@code{@{@var{m},@var{n}@}} to match between @var{m} and @var{n} times.\n\
@end table\n\
@item [@dots{}] [^@dots{}]\n\
List operators, where for example @code{[ab]c} matches @code{ac} and @code{bc}\n\
@item ()\n\
Grouping operator\n\
@item |\n\
Alternation operator. Match one of a choice of regular expressions. The\n\
alternatives must be delimited by the grouoing operator @code{()} above\n\
@item ^ $\n\
Anchoring operator. @code{^} matches the start of the string @var{str} and\n\
@code{$} the end\n\
@end table\n\
\n\
In addition the following escaped characters have special meaning. It should\n\
be noted that it is recommended to quote @var{pat} in single quotes rather\n\
than double quotes, to avoid the escape sequences being interpreted by octave\n\
before being passed to @code{regexp}.\n\
\n\
@table @code\n\
@item \\b\n\
Match a word boundary\n\
@item \\B\n\
Match within a word\n\
@item \\w\n\
Matches any word character\n\
@item \\W\n\
Matches any non word character\n\
@item \\<\n\
Matches the beginning of a word\n\
@item \\>\n\
Matches the end of a word\n\
@item \\s\n\
Matches any whitespace character\n\
@item \\S\n\
Matches any non whitespace character\n\
@item \\d\n\
Matches any digit\n\
@item \\D\n\
Matches any non-digit\n\
@end table\n\
\n\
The outputs of @code{regexp} by default are in the order as given below\n\
\n\
@table @asis\n\
@item @var{s}\n\
The start indices of each of the matching substrings\n\
\n\
@item @var{e}\n\
The end indices of each matching substring\n\
\n\
@item @var{te}\n\
The extents of each of the matched token surrounded by @code{(@dots{})} in\n\
@var{pat}.\n\
\n\
@item @var{m}\n\
A cell array of the text of each match.\n\
\n\
@item @var{t}\n\
A cell array of the text of each token matched.\n\
\n\
@item @var{nm}\n\
A structure containing the text of each matched named token, with the name\n\
being used as the fieldname. A named token is denoted as\n\
@code{(?<name>@dots{})}\n\
@end table\n\
\n\
Particular output arguments or the order of the output arguments can be\n\
selected by additional @var{opts} arguments. These are strings and the\n\
correspondence between the output arguments and the optional argument\n\
are\n\
\n\
@multitable @columnfractions 0.2 0.3 0.3 0.2\n\
@item @tab 'start'        @tab @var{s}  @tab\n\
@item @tab 'end'          @tab @var{e}  @tab\n\
@item @tab 'tokenExtents' @tab @var{te} @tab\n\
@item @tab 'match'        @tab @var{m}  @tab\n\
@item @tab 'tokens'       @tab @var{t}  @tab\n\
@item @tab 'names'        @tab @var{nm}  @tab\n\
@end multitable\n\
\n\
A further optional argument is 'once', that limits the number of returned\n\
matches to the first match. Additional arguments are\n\
\n\
@table @asis\n\
@item matchcase\n\
Make the matching case sensitive.\n\
@item ignorecase\n\
Make the matching case insensitive.\n\
@item stringanchors\n\
Match the anchor characters at the beginning and end of the string.\n\
@item lineanchors\n\
Match the anchor characters at the beginning and end of the line.\n\
@item dotall\n\
The character @code{.} matches the newline character.\n\
@item dotexceptnewline\n\
The character @code{.} matches all but the newline character.\n\
@item freespacing\n\
The pattern can include arbitrary whitespace and comments starting with\n\
@code{#}.\n\
@item literalspacing\n\
The pattern is taken literally.\n\
@end table\n\
@end deftypefn")
{
  return octregexp (args, nargout, "regexp", false);
}

/*

## seg-fault test
%!assert(regexp("abcde","."),[1,2,3,4,5])

## Check that anchoring of pattern works correctly
%!assert(regexp('abcabc','^abc'),1);
%!assert(regexp('abcabc','abc$'),4);
%!assert(regexp('abcabc','^abc$'),[]);

%!test
%! [s, e, te, m, t] = regexp(' No Match ', 'f(.*)uck');
%! assert (s,[])
%! assert (e,[])
%! assert (te,{})
%! assert (m, {})
%! assert (t, {})

%!test
%! [s, e, te, m, t] = regexp(' FiRetrUck ', 'f(.*)uck');
%! assert (s,[])
%! assert (e,[])
%! assert (te,{})
%! assert (m, {})
%! assert (t, {})

%!test
%! [s, e, te, m, t] = regexp(' firetruck ', 'f(.*)uck');
%! assert (s,2)
%! assert (e,10)
%! assert (te{1},[3,7])
%! assert (m{1}, 'firetruck')
%! assert (t{1}{1}, 'iretr')

%!test
%! [s, e, te, m, t] = regexp('short test string','\w*r\w*');
%! assert (s,[1,12])
%! assert (e,[5,17])
%! assert (size(te), [1,2])
%! assert (isempty(te{1}))
%! assert (isempty(te{2}))
%! assert (m{1},'short')
%! assert (m{2},'string')
%! assert (size(t), [1,2])
%! assert (isempty(t{1}))
%! assert (isempty(t{2}))

%!test
%! [s, e, te, m, t] = regexp('short test string','\w*r\w*','once');
%! assert (s,1)
%! assert (e,5)
%! assert (size(te), [1,1])
%! assert (isempty(te{1}))
%! assert (m{1},'short')
%! ## Matlab gives [1,0] here but that seems wrong.
%! assert (size(t), [1,1])

%!test
%! [m, te, e, s, t] = regexp('short test string','\w*r\w*','once', 'match', 'tokenExtents', 'end', 'start', 'tokens');
%! assert (s,1)
%! assert (e,5)
%! assert (size(te), [1,1])
%! assert (isempty(te{1}))
%! assert (m{1},'short')
%! ## Matlab gives [1,0] here but that seems wrong.
%! assert (size(t), [1,1])

%!test
%! ## This test is expected to fail if PCRE is not installed
%! if (!isempty(findstr(octave_config_info ("DEFS"),"HAVE_PCRE")))
%!   [s, e, te, m, t, nm] = regexp('short test string','(?<word1>\w*t)\s*(?<word2>\w*t)');
%!   assert (s,1)
%!   assert (e,10)
%!   assert (size(te), [1,1])
%!   assert (te{1}, [1 5; 7, 10])
%!   assert (m{1},'short test')
%!   assert (size(t),[1,1])
%!   assert (t{1}{1},'short')
%!   assert (t{1}{2},'test')
%!   assert (size(nm), [1,1])
%!   assert (!isempty(fieldnames(nm)))
%!   assert (sort(fieldnames(nm)),{'word1';'word2'})
%!   assert (nm.word1,'short')
%!   assert (nm.word2,'test')
%! endif

%!test
%! ## This test is expected to fail if PCRE is not installed
%! if (!isempty(findstr(octave_config_info ("DEFS"),"HAVE_PCRE")))
%!   [nm, m, te, e, s, t] = regexp('short test string','(?<word1>\w*t)\s*(?<word2>\w*t)', 'names', 'match', 'tokenExtents', 'end', 'start', 'tokens');
%!   assert (s,1)
%!   assert (e,10)
%!   assert (size(te), [1,1])
%!   assert (te{1}, [1 5; 7, 10])
%!   assert (m{1},'short test')
%!   assert (size(t),[1,1])
%!   assert (t{1}{1},'short')
%!   assert (t{1}{2},'test')
%!   assert (size(nm), [1,1])
%!   assert (!isempty(fieldnames(nm)))
%!   assert (sort(fieldnames(nm)),{'word1';'word2'})
%!   assert (nm.word1,'short')
%!   assert (nm.word2,'test')
%! endif

%!test
%! ## This test is expected to fail if PCRE is not installed
%! if (!isempty(findstr(octave_config_info ("DEFS"),"HAVE_PCRE")))
%!   [t, nm] = regexp("John Davis\nRogers, James",'(?<first>\w+)\s+(?<last>\w+)|(?<last>\w+),\s+(?<first>\w+)','tokens','names');
%!   assert (size(t), [1,2]);
%!   assert (t{1}{1},'John');
%!   assert (t{1}{2},'Davis');
%!   assert (t{2}{1},'Rogers');
%!   assert (t{2}{2},'James');
%!   assert (size(nm), [1,1]);
%!   assert (nm.first{1},'John');
%!   assert (nm.first{2},'James');
%!   assert (nm.last{1},'Davis');
%!   assert (nm.last{2},'Rogers');
%! endif

%!assert(regexp("abc\nabc",'.'),[1:7])
%!assert(regexp("abc\nabc",'.','dotall'),[1:7])
%!test
%! if (!isempty(findstr(octave_config_info ("DEFS"),"HAVE_PCRE")))
%!   assert(regexp("abc\nabc",'(?s).'),[1:7])
%!   assert(regexp("abc\nabc",'.','dotexceptnewline'),[1,2,3,5,6,7])
%!   assert(regexp("abc\nabc",'(?-s).'),[1,2,3,5,6,7])
%! endif

%!assert(regexp("caseCaSe",'case'),1)
%!assert(regexp("caseCaSe",'case',"matchcase"),1)
%!assert(regexp("caseCaSe",'case',"ignorecase"),[1,5])
%!test
%! if (!isempty(findstr(octave_config_info ("DEFS"),"HAVE_PCRE")))
%!   assert(regexp("caseCaSe",'(?-i)case'),1)
%!   assert(regexp("caseCaSe",'(?i)case'),[1,5])
%! endif

%!assert (regexp("abc\nabc",'c$'),7)
%!assert (regexp("abc\nabc",'c$',"stringanchors"),7)
%!test
%! if (!isempty(findstr(octave_config_info ("DEFS"),"HAVE_PCRE")))
%!   assert (regexp("abc\nabc",'(?-m)c$'),7)
%!   assert (regexp("abc\nabc",'c$',"lineanchors"),[3,7])
%!   assert (regexp("abc\nabc",'(?m)c$'),[3,7])
%! endif

%!assert (regexp("this word",'s w'),4)
%!assert (regexp("this word",'s w','literalspacing'),4)
%!test
%! if (!isempty(findstr(octave_config_info ("DEFS"),"HAVE_PCRE")))
%!   assert (regexp("this word",'(?-x)s w','literalspacing'),4)
%!   assert (regexp("this word",'s w','freespacing'),[])
%!   assert (regexp("this word",'(?x)s w'),[])
%! endif

%!error regexp('string', 'tri', 'BadArg');
%!error regexp('string');

*/

DEFUN_DLD(regexpi, args, nargout,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {[@var{s}, @var{e}, @var{te}, @var{m}, @var{t}, @var{nm}] =} regexpi (@var{str}, @var{pat})\n\
@deftypefnx {Loadable Function} {[@dots{}] =} regexpi (@var{str}, @var{pat}, @var{opts}, @dots{})\n\
\n\
Case insensitive regular expression string matching. Matches @var{pat} in\n\
@var{str} and returns the position and matching substrings or empty values\n\
if there are none. See @code{regexp} for more details\n\
@end deftypefn")
{
  return octregexp (args, nargout, "regexp", true);
}

/*

## seg-fault test
%!assert(regexpi("abcde","."),[1,2,3,4,5])

## Check that anchoring of pattern works correctly
%!assert(regexpi('abcabc','^abc'),1);
%!assert(regexpi('abcabc','abc$'),4);
%!assert(regexpi('abcabc','^abc$'),[]);

%!test
%! [s, e, te, m, t] = regexpi(' No Match ', 'f(.*)uck');
%! assert (s,[])
%! assert (e,[])
%! assert (te,{})
%! assert (m, {})
%! assert (t, {})

%!test
%! [s, e, te, m, t] = regexpi(' FiRetrUck ', 'f(.*)uck');
%! assert (s,2)
%! assert (e,10)
%! assert (te{1},[3,7])
%! assert (m{1}, 'FiRetrUck')
%! assert (t{1}{1}, 'iRetr')

%!test
%! [s, e, te, m, t] = regexpi(' firetruck ', 'f(.*)uck');
%! assert (s,2)
%! assert (e,10)
%! assert (te{1},[3,7])
%! assert (m{1}, 'firetruck')
%! assert (t{1}{1}, 'iretr')

%!test
%! [s, e, te, m, t] = regexpi('ShoRt Test String','\w*r\w*');
%! assert (s,[1,12])
%! assert (e,[5,17])
%! assert (size(te), [1,2])
%! assert (isempty(te{1}))
%! assert (isempty(te{2}))
%! assert (m{1},'ShoRt')
%! assert (m{2},'String')
%! assert (size(t), [1,2])
%! assert (isempty(t{1}))
%! assert (isempty(t{2}))

%!test
%! [s, e, te, m, t] = regexpi('ShoRt Test String','\w*r\w*','once');
%! assert (s,1)
%! assert (e,5)
%! assert (size(te), [1,1])
%! assert (isempty(te{1}))
%! assert (m{1},'ShoRt')
%! ## Matlab gives [1,0] here but that seems wrong.
%! assert (size(t), [1,1])

%!test
%! [m, te, e, s, t] = regexpi('ShoRt Test String','\w*r\w*','once', 'match', 'tokenExtents', 'end', 'start', 'tokens');
%! assert (s,1)
%! assert (e,5)
%! assert (size(te), [1,1])
%! assert (isempty(te{1}))
%! assert (m{1},'ShoRt')
%! ## Matlab gives [1,0] here but that seems wrong.
%! assert (size(t), [1,1])

%!test
%! ## This test is expected to fail if PCRE is not installed
%! if (!isempty(findstr(octave_config_info ("DEFS"),"HAVE_PCRE")))
%!   [s, e, te, m, t, nm] = regexpi('ShoRt Test String','(?<word1>\w*t)\s*(?<word2>\w*t)');
%!   assert (s,1)
%!   assert (e,10)
%!   assert (size(te), [1,1])
%!   assert (te{1}, [1 5; 7, 10])
%!   assert (m{1},'ShoRt Test')
%!   assert (size(t),[1,1])
%!   assert (t{1}{1},'ShoRt')
%!   assert (t{1}{2},'Test')
%!   assert (size(nm), [1,1])
%!   assert (!isempty(fieldnames(nm)))
%!   assert (sort(fieldnames(nm)),{'word1';'word2'})
%!   assert (nm.word1,'ShoRt')
%!   assert (nm.word2,'Test')
%! endif

%!test
%! ## This test is expected to fail if PCRE is not installed
%! if (!isempty(findstr(octave_config_info ("DEFS"),"HAVE_PCRE")))
%!   [nm, m, te, e, s, t] = regexpi('ShoRt Test String','(?<word1>\w*t)\s*(?<word2>\w*t)', 'names', 'match', 'tokenExtents', 'end', 'start', 'tokens');
%!   assert (s,1)
%!   assert (e,10)
%!   assert (size(te), [1,1])
%!   assert (te{1}, [1 5; 7, 10])
%!   assert (m{1},'ShoRt Test')
%!   assert (size(t),[1,1])
%!   assert (t{1}{1},'ShoRt')
%!   assert (t{1}{2},'Test')
%!   assert (size(nm), [1,1])
%!   assert (!isempty(fieldnames(nm)))
%!   assert (sort(fieldnames(nm)),{'word1';'word2'})
%!   assert (nm.word1,'ShoRt')
%!   assert (nm.word2,'Test')
%! endif

%!assert(regexpi("abc\nabc",'.'),[1:7])
%!assert(regexpi("abc\nabc",'.','dotall'),[1:7])
%!test
%! if (!isempty(findstr(octave_config_info ("DEFS"),"HAVE_PCRE")))
%!   assert(regexpi("abc\nabc",'(?s).'),[1:7])
%!   assert(regexpi("abc\nabc",'.','dotexceptnewline'),[1,2,3,5,6,7])
%!   assert(regexpi("abc\nabc",'(?-s).'),[1,2,3,5,6,7])
%! endif

%!assert(regexpi("caseCaSe",'case'),[1,5])
%!assert(regexpi("caseCaSe",'case',"matchcase"),1)
%!assert(regexpi("caseCaSe",'case',"ignorecase"),[1,5])
%!test
%! if (!isempty(findstr(octave_config_info ("DEFS"),"HAVE_PCRE")))
%!   assert(regexpi("caseCaSe",'(?-i)case'),1)
%!   assert(regexpi("caseCaSe",'(?i)case'),[1,5])
%! endif

%!assert (regexpi("abc\nabc",'c$'),7)
%!assert (regexpi("abc\nabc",'c$',"stringanchors"),7)
%!test
%! if (!isempty(findstr(octave_config_info ("DEFS"),"HAVE_PCRE")))
%!   assert (regexpi("abc\nabc",'(?-m)c$'),7)
%!   assert (regexpi("abc\nabc",'c$',"lineanchors"),[3,7])
%!   assert (regexpi("abc\nabc",'(?m)c$'),[3,7])
%! endif

%!assert (regexpi("this word",'s w'),4)
%!assert (regexpi("this word",'s w','literalspacing'),4)
%!test
%! if (!isempty(findstr(octave_config_info ("DEFS"),"HAVE_PCRE")))
%!   assert (regexpi("this word",'(?-x)s w','literalspacing'),4)
%!   assert (regexpi("this word",'s w','freespacing'),[])
%!   assert (regexpi("this word",'(?x)s w'),[])
%! endif

%!error regexpi('string', 'tri', 'BadArg');
%!error regexpi('string');

*/

/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/