%% Copyright (C) 2000 Paul Kienzle
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

%% -*- texinfo -*-
%% @deftypefn {Function File} {} now
%% Returns the current local time as a number of days since Jan 1, 0000.
%% By this reckoning, Jan 1, 1970 is day number 719529.  The fractional
%% portion, @code{rem(now,1)} corresponds to the portion of the current
%% day.
%%
%% @seealso{date,clock,datenum,datestr,datevec,calendar,weekday}
%% @end deftypefn

function t = now
  t = datenum(clock);
  %% The following doesn't work (e.g., one hour off on 2005-10-04):
  %%   seconds since 1970-1-1 corrected by seconds from GMT to local time
  %%   divided by 86400 sec/day plus day num for 1970-1-1
  %%   t = (time - mktime(gmtime(0)))/86400 + 719529;
  %% mktime(gmtime(0)) does indeed return the offset from Greenwich to the
  %% local time zone, but we need to account for daylight savings time
  %% changing by an hour the offset from CUT for part of the year.
