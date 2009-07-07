// This code is public domain.
// Author: Paul Kienzle
#ifdef HAVE_CONFIG_H
# include <config.h>
#else
# include <octave/config.h>
#endif

#include "defun-dld.h"
#include "error.h"
#include "parse.h"
#include "quit.h"
#include "Cell.h"

DEFUN_DLD(regexprep,args,nargout,"\
-*- texinfo -*-\n\
@deftypefn {Function File}  @var{string} = regexprep(@var{string}, @var{pat}, @var{repstr}, @var{options})\n\
Replace matches of @var{pat} in  @var{string} with @var{repstr}.\n\
\n\
\n\
The replacement can contain @code{$i}, which subsubstitutes\n\
for the ith set of parentheses in the match string.  E.g.,\n\
@example\n\
\n\
   regexprep(\"Bill Dunn\",'(\\w+) (\\w+)','$2, $1')\n\
\n\
@end example\n\
returns \"Dunn, Bill\"\n\
\n\
@var{options} may be zero or more of\n\
@table @samp\n\
\n\
@item once\n\
Replace only the first occurance of @var{pat} in the result.\n\
\n\
@item warnings\n\
This option is present for compatibility but is ignored.\n\
\n\
@item ignorecase or matchcase\n\
Ignore case for the pattern matching (see @code{regexpi}).\n\
Alternatively, use (?i) or (?-i) in the pattern.\n\
\n\
@item lineanchors and stringanchors\n\
Whether characters ^ and $ match the beginning and ending of lines.\n\
Alternatively, use (?m) or (?-m) in the pattern.\n\
\n\
@item dotexceptnewline and dotall\n\
Whether . matches newlines in the string.\n\
Alternatively, use (?s) or (?-s) in the pattern.\n\
\n\
@item freespacing or literalspacing\n\
Whether whitespace and # comments can be used to make the regular expression more readable.\n\
Alternatively, use (?x) or (?-x) in the pattern.\n\
\n\
@end table\n\
@seealso{regexp,regexpi}\n\
@end deftypefn")
{
  octave_value_list retval;

  int nargin = args.length();
  int nopts = nargin - 3;

  if (nargin < 3)
    {
      print_usage("regexprep");
      return retval;
    }

  // Make sure we have string,pattern,replacement
  const std::string buffer = args(0).string_value ();
  if (error_state) return retval;
  const std::string pattern = args(1).string_value ();
  if (error_state) return retval;
  const std::string replacement = args(2).string_value ();
  if (error_state) return retval;
  
  // Pack options excluding 'tokenize' and various output
  // reordering strings into regexp arg list
  octave_value_list regexpargs(nargin-1,octave_value());
  regexpargs(0) = args(0);
  regexpargs(1) = args(1);
  int len=2;
  for (int i = 3; i < nargin; i++) 
    {
      const std::string opt = args(i).string_value();
      if (opt != "tokenize" && opt != "start" && opt != "end"
	  && opt != "tokenextents" && opt != "match" && opt != "tokens"
	  && opt != "names"  && opt != "warnings") 
	{
	  regexpargs(len++) = args(i);
	}
    }
  regexpargs.resize(len);
  
  //std::cout << "Buffer " << buffer << std::endl;
  //std::cout << "Pattern " << pattern << std::endl;
  //std::cout << "Replacement " << replacement << std::endl;

  // Identify replacement tokens; build a vector of group numbers in
  // the replacement string so that we can quickly calculate the size 
  // of the replacement.
  int tokens = 0;
  for (size_t i=1; i < replacement.size(); i++) 
    {
      if (replacement[i-1]=='$' && isdigit(replacement[i])) 
	{
	  tokens++, i++;
	}
    }
  std::vector<int> token(tokens);
  int k=0;
  for (size_t i=1; i < replacement.size(); i++) 
    {
      if (replacement[i-1]=='$' && isdigit(replacement[i])) 
	{
	  token[k++] = replacement[i]-'0';
	  i++;
	}
    }

  // Perform replacement
  std::string rep;
  if (tokens > 0) 
    {
      
      // Call regexp, asking for start, end, and capture start/end
      octave_value_list regexpret = feval("regexp", regexpargs, 3);
      if (regexpret(0).is_empty()) 
	{
	  retval(0) = args(0);
	  return retval;
	}
      const ColumnVector s(regexpret(0).vector_value());
      const ColumnVector e(regexpret(1).vector_value());
      const Cell te(regexpret(2).cell_value());
      if (error_state) return retval;
      
      // Determine replacement length
      const size_t replen = replacement.size() - 2*tokens;
      int delta = 0;
      for (int i=0; i < s.length(); i++) 
	{
	  OCTAVE_QUIT;

	  const Matrix pairs(te(i).matrix_value());
	  size_t pairlen = 0;
	  for (int j=0; j < tokens; j++) 
	    {
	      if (token[j] == 0) 
		pairlen += size_t(e(i)-s(i))+1;
	      else if (token[j] <= pairs.rows()) 
		pairlen += size_t(pairs(token[j]-1,1)-pairs(token[j]-1,0))+1;
	    }
	  delta += int(replen + pairlen) - int(e(i)-s(i)+1);
	}
      
      // std::cout << "replacement delta is " << delta << std::endl;
      
      // Build replacement string
      rep.reserve(buffer.size()+delta);
      size_t from = 0;
      for (int i=0; i < s.length(); i++) 
	{
	  OCTAVE_QUIT;

	  const Matrix pairs(te(i).matrix_value());
	  rep.append(&buffer[from], size_t(s(i)-1)-from);
	  from = size_t(e(i)-1)+1;
	  for (size_t j=1; j < replacement.size(); j++) 
	    {
	      if (replacement[j-1]=='$' && isdigit(replacement[j])) 
		{
		  int k = replacement[j]-'0';
		  if (k==0) 
		    { 
		      // replace with entire match
		      rep.append(&buffer[size_t(e(i)-1)],
				 size_t(e(i)-s(i))+1);
		    } 
		  else if (k <= pairs.rows()) 
		    {
		      // replace with group capture
		      // std::cout << "k=" << k << " [" << pairs(k-1,0) << "," << pairs(k-1,1) << "]" << std::endl;
		      rep.append(&buffer[size_t(pairs(k-1,0)-1)],
				 size_t(pairs(k-1,1)-pairs(k-1,0))+1);
		    } 
		  else 
		    {
		      // replace with nothing
		    }
		  j++;
		} 
	      else 
		{
		  rep.append(1,replacement[j-1]);
		  if (j == replacement.size()-1) 
		    {
		      rep.append(1,replacement[j]);
		    }
		}
	    }
	  // std::cout << "->" << rep << std::endl;
	}
      rep.append(&buffer[from],buffer.size()-from);
      
    } 
  else 
    {

      // Call regexp, asking for start, end
      octave_value_list regexpret = feval("regexp", regexpargs, 2);
      if (regexpret(0).is_empty()) 
	{
	  retval(0) = args(0);
	  return retval;
	}

      const ColumnVector s(regexpret(0).vector_value());
      const ColumnVector e(regexpret(1).vector_value());
      if (error_state) return retval;

      // Determine replacement length
      const size_t replen = replacement.size();
      int delta = 0;
      for (int i=0; i < s.length(); i++) 
	{
          OCTAVE_QUIT;
	  delta += int(replen) - int(e(i)-s(i)+1);
	}

      // std::cout << "replacement delta is " << delta << std::endl;
      
      // Build replacement string
      rep.reserve(buffer.size()+delta);
      size_t from = 0;
      for (int i=0; i < s.length(); i++) 
	{
          OCTAVE_QUIT;
	  rep.append(&buffer[from],size_t(s(i)-1)-from);
	  from = size_t(e(i)-1)+1;
	  rep.append(replacement);
	  // std::cout << "->" << rep << std::endl;
	}
      rep.append(&buffer[from],buffer.size()-from);
      
    }
  
  retval(0) = rep;
  return retval;
}

/*
%!test  # Replace with empty
%! xml = '<!-- This is some XML --> <tag v="hello">some stuff<!-- sample tag--></tag>';
%! t = regexprep(xml,'<[!?][^>]*>','');
%! assert(t,' <tag v="hello">some stuff</tag>')

%!test  # Replace with non-empty
%! xml = '<!-- This is some XML --> <tag v="hello">some stuff<!-- sample tag--></tag>';
%! t = regexprep(xml,'<[!?][^>]*>','?');
%! assert(t,'? <tag v="hello">some stuff?</tag>')

%!test  # Check that 'tokenize' is ignored
%! xml = '<!-- This is some XML --> <tag v="hello">some stuff<!-- sample tag--></tag>';
%! t = regexprep(xml,'<[!?][^>]*>','','tokenize');
%! assert(t,' <tag v="hello">some stuff</tag>')

%!test  # Capture replacement
%! if (!isempty(findstr(octave_config_info ("DEFS"),"HAVE_PCRE")))
%!   data = "Bob Smith\nDavid Hollerith\nSam Jenkins";
%!   result = "Smith, Bob\nHollerith, David\nJenkins, Sam";
%!   t = regexprep(data,'(?m)^(\w+)\s+(\w+)$','$2, $1');
%!   assert(t,result)
%! end

# Return the original if no match
%!assert(regexprep('hello','world','earth'),'hello')

## Test a general replacement
%!assert(regexprep("a[b]c{d}e-f=g", "[^A-Za-z0-9_]", "_"), "a_b_c_d_e_f_g");

## Make sure it works at the beginning and end
%!assert(regexprep("a[b]c{d}e-f=g", "a", "_"), "_[b]c{d}e-f=g");
%!assert(regexprep("a[b]c{d}e-f=g", "g", "_"), "a[b]c{d}e-f=_");

## Options
%!assert(regexprep("a[b]c{d}e-f=g", "[^A-Za-z0-9_]", "_", "once"), "a_b]c{d}e-f=g");
%!assert(regexprep("a[b]c{d}e-f=g", "[^A-Z0-9_]", "_", "ignorecase"), "a_b_c_d_e_f_g");

## Option combinations
%!assert(regexprep("a[b]c{d}e-f=g", "[^A-Z0-9_]", "_", "once", "ignorecase"), "a_b]c{d}e-f=g");

*/