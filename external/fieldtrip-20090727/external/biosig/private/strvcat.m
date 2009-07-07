%% Copyright (C) 1996 Kurt Hornik
%%
%% This program is free software; you can redistribute it and/or modify it
%% under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2, or (at your option)
%% any later version.
%%
%% This program is distributed in the hope that it will be useful, but
%% WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%% General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; see the file COPYING.  If not, write to the Free
%% Software Foundation, 59 Temple Place - Suite 330, Boston, MA
%% 02111-1307, USA.

%% -*- texinfo -*-
%% @deftypefn {Function File} {} strvcat (@var{s_1}, @dots{}, @var{s_n})
%% Return a matrix containing the strings (and cell-strings) 
%% @var{s_1}, @dots{}, @var{s_n} as
%% its rows.  Each string is padded with blanks in order to form a valid
%% matrix.  Unlike @var{str2mat}, empty strings are ignored.
%%
%% @end deftypefn

%% Author: Kurt Hornik <Kurt.Hornik@ci.tuwien.ac.at>
%% Adapted-By: jwe
%% Modified: Paul Kienzle <pkienzle@kienzle.powernet.co.uk> converted
%%           str2mat to strvcat.  Same function except that strvcat
%%           ignores empty strings.
%% Modified by Alois Schloegl <a.schloegl@ieee.org> Mar 2005
%%	     added support for cell-strings 


function retval = strvcat (varargin)

  if (nargin == 0)
    usage ('strvcat (s1, ...)');
  end

  va_arg_cnt = 1;

  nr = zeros(nargin, 1);
  nc = zeros(nargin, 1);
  K = 0; 
  for k = 1:nargin,
    s = nth(varargin, va_arg_cnt);
    va_arg_cnt = va_arg_cnt+1;
    if iscell(s),
	for k1 = 1:length(s)
	    K = K+1;
	    [nr(K), nc(K)] = size (s{k1});
	end;
    else
        K = K+1;
        [nr(K), nc(K)] = size (s);
    end;	
  end

  retval_nr = sum (nr);
  retval_nc = max (nc);

  retval = setstr (ones(retval_nr, retval_nc) * 32);

  va_arg_cnt = 1;

  row_offset = 0;
  for k = 1 : nargin
    s = nth (varargin, va_arg_cnt);
    va_arg_cnt = va_arg_cnt+1;
    if iscell(s),
	for k1 = 1:length(s)
	    retval(row_offset+k1,1:length(s{k1})) = char(s{k1});
	end;
        row_offset = row_offset + length(s);
    else
	if (nc(k) > 0)
    	    retval ((row_offset + 1) : (row_offset + nr(k)), 1:nc(k)) = char(s);
	end
        row_offset = row_offset + size(s,1);
    end;
  end

end
