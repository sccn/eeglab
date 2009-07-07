/*

Copyright (C) Peter Rydesäter 2002, Mitthögskolan, SWEDEN
Copyright (C) 2005 David Bateman

This file is part of Octave.

Octave is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2, or (at your option) any
later version.

Octave is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with Octave; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.

*/

//
//   Please, send me an e-mail if you use this and/or improve it. 
//   If you think this is usful fore your work, please remember me and my
//   hard work with it!
//
//   Peter.Rydesater@mh.se
//
//   Peter Rydesäter
//   Mitthögskolan (Mid Sweden University)
//   SE-831 25 ÖSTERSUND
//   SWEDEN
//

#include <string>
#include <octave/config.h>
#include <octave/ov.h>
#include <octave/defun-dld.h>
#include <octave/error.h>
#include <octave/ov-scalar.h>
#include <octave/oct-map.h>

static std::string matstr;
static int matlen =0;
static int matpos = 0;
#define STREND (-1)

static int iseof (void)
{
  return (matpos >= matlen || matpos < 0);
}

static int nextch (void)
{
  if (iseof())
    return STREND;
  else
    return matstr[matpos];
}

static int readch (void)
{
  if (iseof())
    return STREND;
  else
    return matstr[matpos++];
}

static void unreadch (int ch)
{
  if (matpos > 0)
    matpos--;
}

static void skipblank (void)
{
  while (!iseof())
    {
      int ch = readch();
      if (!isspace(ch))
	{
	  unreadch(ch);
	  break;
	}
    }
}

static void readlabel (std::string &str)
{
  str = std::string ();
  int n = 0;
  while (!iseof())
    {
      int ch = readch();
      if(ch == '=' || isspace(ch) || ch == '>' || ch == '/' || iseof())
	{
	  unreadch(ch);
	  break;
	}
      if (!isalnum(ch)) 
	{
	  ch='_';    
	  if (n == 0) 
	    {
	      str += 'x';
	      n++;
	    }
	}
      str += ch;
      n++;
    }
}

static void readcom (std::string &str)
{
  while (!iseof())
    {
      int ch = readch();
      str += ch;
      if (ch == '>')
	break;	
    }
}

static void readtext (std::string &str)
{
  str = std::string ();
  int n = 0;
  skipblank();                 //remove blank;
  while (!iseof())
    {
      if (nextch() == '<')
	break;
      str += readch ();
      n++;
    }

  int m = n;
  while (m > 0 && isspace(str[m-1])) // remove blank;
    m--;
  if (m != n)
    str = str.substr(0, m);
}

static void readparvalue (std::string &str)
{
  str = std::string ();
  skipblank();
  if (nextch() == '"')
    {
      readch();
      while (!iseof())
	{
	  int ch = readch();
	  if(ch == '"')
	    break;
	  str += ch;
	}
    }
  else
    {
      while (!iseof()) 
	{
	  int ch = nextch();
	  if (isspace(ch) || ch == '/' || ch == '>')
	    break;
	  str += readch();
	}
    }
}

DEFUN_DLD (__xmldata__, args, nargout,
  "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {[@var{name}, @var{data}, @var{type}, @var{p}} = xmldata (@var{str}, @var{p})\n\
\n\
Read a XML tag or text at a given position in a XML string. This is an\n\
internal function and should not be called directly.\n\
\n\
@seealso{xmlstruct}\n\
@end deftypefn")
{
  int nargin = args.length ();
  octave_value_list retval;

  retval(3) = octave_value (double (0));
  retval(2) = octave_value (double (0));
  retval(1) = octave_value (double (0));
  retval(0) = octave_value (double (0));

  if (nargin < 2)
    return retval;
  if (!args(0).is_string())
    return retval;
  matstr = args(0).string_value();
  matlen = matstr.length();
  matpos = static_cast<int>(args(1).scalar_value()) - 1;
  if (nargout < 3)
    return retval;
  if (!iseof())
    {
      std::string str;
      readtext(str);
      if (str.length () == 0)
	{
	  std::string buff;
	  int tagtype = 1;
	  Octave_map map (dim_vector (1,1));

	  while (readch() != '<')
	    if (iseof()) 
	      goto finished;
	  skipblank();
	  if (nextch() == '/')
	    {
	      tagtype = 3;
	      readch();
	      skipblank();
	    }
	  if (nextch() == '!')
	    {
	      buff = '<';
	      readcom (buff);
	      retval(0) = octave_value (std::string(buff));
	      retval(1) = octave_value (std::string());
	      retval(2) = octave_value (double (0));
	      goto finished;
	    }
	  readlabel(buff);
	  if (tagtype == 3)
	    {
	      retval(0) = octave_value (std::string(buff));
	      retval(1) = octave_value (std::string());
	    }
	  else
	    {
	      retval(0) = octave_value (std::string(buff));
	      retval(1) = map;
	    }
	  while (!iseof())
	    {
	      skipblank();
	      if (nextch() == '/')
		{
		  readch();
		  tagtype = 2;
		}
	      skipblank();
	      if (nextch() == '>')
		{
		  readch();
		  retval(2) = octave_value ((double)tagtype);
		  break;
		}
	      std::string struct_el;
	      readlabel(struct_el);
	      map.contents(struct_el);
	      retval(1) = map;
	      skipblank();
	      if (nextch() == '=')
		{
		  readch();
		  readparvalue(buff);
		  map.assign(struct_el, buff);
		  retval(1) = map;
		}
	    }
	}
      else
	{
	  retval(0) = octave_value (std::string(str));
	  retval(1) = octave_value (std::string());
	  retval(2) = octave_value (double (0));
	}
    }

 finished:
  retval(3) = octave_value (double ((matpos+1)*(!iseof())));

  return retval;
}
