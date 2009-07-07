/*

Copyright (C) 2006 David Bateman

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; see the file COPYING.  If not, write to the Free
Software Foundation, 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

In addition to the terms of the GPL, you are permitted to link
this program with any Open Source program, as defined by the
Open Source Initiative (www.opensource.org)

*/

#include <octave/config.h>
#include <octave/oct.h>
#include <octave/Cell.h>
#include <octave/ov-str-mat.h>
#include <octave/ov-colon.h>

DEFUN_DLD (mat2cell, args, ,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {@var{b} =} num2cell (@var{a}, @var{m}, @var{n})\n\
@deftypefnx {Loadable Function} {@var{b} =} num2cell (@var{a}, @var{d1}, @var{d2}, @dots{})\n\
@deftypefnx {Loadable Function} {@var{b} =} num2cell (@var{a}, @var{r})\n\
Converts the matrix @var{a} to a cell array If @var{a} is 2-D, then\n\
it is required that @code{sum (@var{m}) == size (@var{a}, 1)} and\n\
@code{sum (@var{n}) == size (@var{a}, 2)}. Similarly, if @var{a} is\n\
a multi-dimensional and the number of dimensional arguments is equal\n\
to the dimensions of @var{a}, then it is required that @code{sum (@var{di})\n\
== size (@var{a}, i)}.\n\
\n\
Given a single dimensional argument @var{r}, the other dimensional\n\
arguments are assumed to equal @code{size (@var{a},@var{i})}.\n\
\n\
An example of the use of mat2cell is\n\
\n\
@example\n\
@group\n\
mat2cell (reshape(1:16,4,4),[3,1],[3,1])\n\
@result{} @{\n\
  [1,1] =\n\
\n\
     1   5   9\n\
     2   6  10\n\
     3   7  11\n\
\n\
  [2,1] =\n\
\n\
     4   8  12\n\
\n\
  [1,2] =\n\
\n\
    13\n\
    14\n\
    15\n\
\n\
  [2,2] = 16\n\
@}\n\
@end group\n\
@end example\n\
@end deftypefn\n\
@seealso{num2cell,cell2mat}")
{
  int nargin = args.length();
  octave_value retval;

  if (nargin < 2)
    usage ("mat2cell");
  else
    {
      dim_vector dv = args(0).dims();
      dim_vector new_dv;
      new_dv.resize(dv.length());
      
      if (nargin > 2)
	{
	  octave_idx_type nmax = -1;

	  if (nargin - 1 != dv.length())
	    error ("mat2cell: Incorrect number of dimensions");
	  else
	    {
	      for (octave_idx_type j = 0; j < dv.length(); j++)
		{
		  ColumnVector d = ColumnVector (args(j+1).vector_value 
						 (false, true));

		  if (d.length() < 1)
		    {
		      error ("mat2cell: dimension can not be empty");
		      break;
		    }
		  else
		    {
		      if (nmax < d.length())
			nmax = d.length();

		      for (octave_idx_type i = 1; i < d.length(); i++)
			{
			  OCTAVE_QUIT;

			  if (d(i) >= 0)
			    d(i) += d(i-1);
			  else
			    {
			      error ("mat2cell: invalid dimensional argument");
			      break;
			    }
			}

		      if (d(0) < 0)
			error ("mat2cell: invalid dimensional argument");
		      
		      if (d(d.length() - 1) != dv(j))
			error ("mat2cell: inconsistent dimensions");

		      if (error_state)
			break;

		      new_dv(j) = d.length();
		    }
		}
	    }

	  if (! error_state)
	    {
	      // Construct a matrix with the index values
	      Matrix dimargs(nmax, new_dv.length());
	      for (octave_idx_type j = 0; j < new_dv.length(); j++)
		{
		  OCTAVE_QUIT;

		  ColumnVector d = ColumnVector (args(j+1).vector_value 
						 (false, true));

		  dimargs(0,j) = d(0);
		  for (octave_idx_type i = 1; i < d.length(); i++)
		    dimargs(i,j) = dimargs(i-1,j) + d(i);
		}


	      octave_value_list lst (new_dv.length(), octave_value());
	      Cell ret (new_dv);
	      octave_idx_type nel = new_dv.numel();
	      octave_idx_type ntot = 1;

	      for (int j = 0; j < new_dv.length()-1; j++)
		ntot *= new_dv(j);

	      for (octave_idx_type i = 0; i <  nel; i++)
		{
		  octave_idx_type n = ntot;
		  octave_idx_type ii = i;
		  for (octave_idx_type j =  new_dv.length() - 1;  j >= 0; j--)
		    {
		      OCTAVE_QUIT;
		  
		      octave_idx_type idx = ii / n;
		      lst (j) = Range((idx == 0 ? 1. : dimargs(idx-1,j)+1.),
				      dimargs(idx,j));
		      ii = ii % n;
		      if (j != 0)
			n /= new_dv(j-1);
		    }
		  ret(i) = args(0).do_index_op(lst, 0);
		  if (error_state)
		    break;
		}
	  
	      if (!error_state)
		retval = ret;
	    }
	}
      else
	{
	  ColumnVector d = ColumnVector (args(1).vector_value 
					 (false, true));

	  double sumd = 0.;
	  for (octave_idx_type i = 0; i < d.length(); i++)
	    {
	      OCTAVE_QUIT;

	      if (d(i) >= 0)
		sumd += d(i);
	      else
		{
		  error ("mat2cell: invalid dimensional argument");
		  break;
		}
	    }

	  if (sumd != dv(0))
	    error ("mat2cell: inconsistent dimensions");

	  new_dv(0) = d.length();
	  for (octave_idx_type i = 1; i < dv.length(); i++)
	    new_dv(i) = 1;

	  if (! error_state)
	    {
	      octave_value_list lst (new_dv.length(), octave_value());
	      Cell ret (new_dv);

	      for (octave_idx_type i = 1; i < new_dv.length(); i++)
		lst (i) = Range (1., static_cast<double>(dv(i)));
	      
	      double idx = 0.;
	      for (octave_idx_type i = 0; i <  new_dv(0); i++)
		{
		  OCTAVE_QUIT;

		  lst(0) = Range(idx + 1., idx + d(i));
		  ret(i) = args(0).do_index_op(lst, 0);
		  idx += d(i);
		  if (error_state)
		    break;
		}
	  
	      if (!error_state)
		retval = ret;
	    }
	}
    }

  return retval;
}

/*

%!test
%! x = reshape(1:20,5,4);
%! c = mat2cell(x,[3,2],[3,1]);
%! assert(c,{[1,6,11;2,7,12;3,8,13],[16;17;18];[4,9,14;5,10,15],[19;20]})

%!test
%! x = 'abcdefghij';
%! c = mat2cell(x,1,[0,4,2,0,4,0]);
%! empty1by0str = resize('',1,0);
%! assert(c,{empty1by0str,'abcd','ef',empty1by0str,'ghij',empty1by0str})

*/
	  
/*
;;; Local Variables: ***
;;; mode: C++ ***
;;; End: ***
*/
