## Copyright (C) 2001 Bill Denney <bill@givebillmoney.com>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

## -*- texinfo -*-
## @deftypefn {Function File} {Y =} datesplit(date, P)
## @deftypefnx {Function File} {[Y,M,D,h,m,s] =} datesplit(date, P)
## Split a date string into the Year, Month, Day, hour, minute, and
## second.  This routine tries to be as forgiving as possible to the
## date input while requiring that the date is not ambiguous.
##
## Anywhere possible where it would not be ambiguous, efforts were made
## to make times possible with seconds and AM/PM as optional.  Also,
## along the same lines, where possible, commas were allowed with
## spaces, and the year/month/day separators were allowed as period (.),
## slash (/), and dash (-).  Not all format possibilities are shown in
## the following table, but a date like @code{dd-mmm-yyyy HH:MM:SS} is
## parsed just as well as @code{d/mmm.yyyy,  ,H:MM, AM}.
##
## Supported @code{date} formats include (the same as datestr):
## @multitable @columnfractions 0.1 0.45 0.45
## @item @strong{Code} @tab @strong{Format} @tab @strong{Example}
## @item  0  @tab dd-mmm-yyyy HH:MM:SS    @tab 07-Sep-2000 15:38:09
## @item  1  @tab dd-mmm-yyyy             @tab 07-Sep-2000 
## @item  2  @tab mm/dd/yy                @tab 09/07/00 
## @item  3  @tab mmm                     @tab Sep 
## @item  6  @tab mm/dd                   @tab 09/13 
## @item 10  @tab yyyy                    @tab 2000 
## @item 12  @tab mmmyy                   @tab Sep00 
## @item 13  @tab HH:MM:SS                @tab 15:38:09 
## @item 14  @tab HH:MM:SS PM             @tab 03:38:09 PM
## @item 15  @tab HH:MM                   @tab 15:38 
## @item 16  @tab HH:MM PM                @tab 03:38 PM 
## @item 17  @tab QQ-YY                   @tab Q3-00
## @item 19  @tab dd/mm                   @tab 13/03
## @item 20  @tab dd/mm/yy                @tab 13/03/95
## @item 21  @tab mmm.dd.yyyy HH:MM:SS    @tab Mar.03.1962 13:53:06
## @item 22  @tab mmm.dd.yyyy             @tab Mar.03.1962
## @item 23  @tab mm/dd/yyyy              @tab 03/13/1962
## @item 24  @tab dd/mm/yyyy              @tab 12/03/1962
## @item 25  @tab yy/mm/dd                @tab 95/03/13
## @item 26  @tab yyyy/mm/dd              @tab 1995/03/13
## @item 27  @tab QQ-YYYY                 @tab Q4-2132
## @item 28  @tab mmmyyyy                 @tab Mar2047
## @item 29  @tab yyyymmdd                @tab 20470313
## @item 30  @tab yyyymmddTHHMMSS         @tab 20470313T132603
## @item 31  @tab yyyy-mm-dd HH:MM:SS     @tab 1047-03-13 13:26:03
## @end multitable
##
## The parameter @code{P} is needed to convert date strings with 2 digit
## years into dates with 4 digit years.  2 digit years are assumed to be
## between @code{P} and @code{P+99}. If @code{P} is not given then the 
## current year - 50 is used, so that dates are centered on the present.
## For birthdates, you would want @code{P} to be current year - 99.  For
## appointments, you would want @code{P} to be current year.
##
## This function makes no strong attempt to verify the accuracy of the
## numbers that it returns in that it doesn't (currently) check to see
## that you're not trying to use the date Feb 30.  When applicable, it
## tries to make your input work, though.  It will try to determine if
## you're using the date "03/13/95" that the date is "March 13, 1995",
## but if there is doubt, datesplit will return an error instead of
## trying to guess the wrong value.
##
## @seealso{date,clock,now,datestr,datenum,calendar,weekday} 
## @end deftypefn

## TODO:
##  * Some formats are ambiguous.  Allow the user to specify the format
##    to remove ambiguity.
##  * Validate the dates.
##  * Possible bug (after dates are validated): There are times where
##    the year is assumed, Feb 29 may be a valid date, but with the
##    assumed year, it may become invalid.
##  * Internationalize.  Not all months use the English system.
##  * Vectorize.  That requires vectorization of regexp though...

## Author: Bill Denney <bill@givebillmoney.com>

function [y, m, d, h, mi, s] = datesplit(ds, P)

  if nargin < 2
    P = [];
  endif

  today = datevec(now);

  if (isempty(P))
    P = today(1)-50;
  endif

  global __month_names = ["Jan";"Feb";"Mar";"Apr";"May";"Jun";...
			  "Jul";"Aug";"Sep";"Oct";"Nov";"Dec"];
  global __day_names   = ["Sun";"Mon";"Tue";"Wed";"Thu";"Fri";"Sat"];
  global __time_names  = ["AM";"PM"];

  if (iscellstr(ds))
    ds = ds{1};
  endif
  ds = tolower(deblank(ds));

  if (nargin < 1)
    error("datesplit: no input arguments");
  elseif (nargin == 1)
    fmt = [];
  endif
  %% we have to determine the format, this could be error prone

  ## format  0  dd-mmm-yyyy HH:MM:SS    e.g. 07-Sep-2000 15:38:09
  [match, d, m, y, h, mi, s, ap] = \
      regexp("^(3[01]|[0-2]?[0-9])[-./]([a-z]{3})[-./]([0-9]{4})[, ]+(2[0-3]|[01]?[0-9]):([0-5][0-9])(:[0-5][0-9]|)[, ]*([ap]m|)$", ds);

  ## format 21  mmm.dd.yyyy HH:MM:SS    e.g. Mar.03.1962 13:53:06
  if (isempty(match))
    [match, m, d, y, h, mi, s, ap] = \
	regexp("^([a-z]{3})[-./](3[01]|[0-2]?[0-9])[-./]([0-9]{4})[, ]+(2[0-3]|[01]?[0-9]):([0-5][0-9])(:[0-5][0-9]|)[, ]*([ap]m|)$", ds);
  endif

  ## format 31  yyyy-mm-dd HH:MM:SS     e.g. 2004-03-13 13:26:03
  if (isempty(match))
    [match, y, m, d, h, mi, s, ap] = \
	regexp("^([0-9]{4})[-./](1[012]|0?[0-9])[-./](3[01]|[0-2]?[0-9])[, ]+(2[0-3]|[01]?[0-9]):([0-5][0-9])(:[0-5][0-9]|)[, ]*([ap]m|)$", ds);
  endif

  ## format 30  yyyymmddTHHMMSS         e.g. 20470313T132603
  if (isempty(match))
    [match, y, m, d, h, mi, s] = \
	regexp("^([0-9]{4})(1[012]|0[0-9])(3[01]|[012][0-9])t(2[0-3]|[01][0-9])([0-5][0-9])([0-5][0-9])$", ds);
    ap = "NA";
  endif

  ## format 13  HH:MM:SS                e.g. 15:38:09
  ## format 14  HH:MM:SS PM             e.g. 03:38:09 PM
  ## format 15  HH:MM                   e.g. 15:38
  ## format 16  HH:MM PM                e.g. 03:38 PM
  if (isempty(match))
    [match, h, mi, s, ap] = \
	regexp("^(2[0-3]|[01]?[0-9]):([0-5][0-9])(:[0-5][0-9]|)[, ]*([ap]m|)$", ds);

    if (! isempty(match))
      %% assume that it is as of today
      y = today(1);
      m = today(2);
      d = today(3);
    endif
  endif

  ## format  1  dd-mmm-yyyy             e.g. 07-Sep-2000
  if (isempty(match))
    [match, d, m, y] = \
	regexp("^(3[01]|[012]?[0-9])[-./]([a-z]{3})[-./]([0-9]{4})$", ds);

    if (! isempty(match))
      %% assume the beginning of the day
      h = 0;
      mi = 0;
      s = 0;
      ap = "NA";
    endif
  endif

  ## format 22  mmm.dd.yyyy             e.g. Mar.03.1962
  if (isempty(match))
    [match, m, d, y] = \
	regexp("^([a-z]{3})[-./](3[01]|[012]?[0-9])[-./]([0-9]{4})$", ds);

    if (! isempty(match))
      %% assume the beginning of the day
      h = 0;
      mi = 0;
      s = 0;
      ap = "NA";
    endif
  endif

  ## format  2  mm/dd/yy                e.g. 09/07/00
  ## format 23  mm/dd/yyyy              e.g. 03/13/1962
  ## format 20  dd/mm/yy                e.g. 13/03/95
  ## format 24  dd/mm/yyyy              e.g. 12/03/1962
  ## format 25  yy/mm/dd                e.g. 95/03/13
  ## format 26  yyyy/mm/dd              e.g. 1995/03/13
  if (isempty(match))
    [match, d, m, y] = \
	regexp("^([0-9]{1,2}|[0-9]{4})[-./](3[01]|[012]?[0-9])[-./]([0-9]{1,2}|[0-9]{4})$", ds);

    if (! isempty(match))
      %% we have to determine if the date is unambiguous
      d = str2num(d);
      m = str2num(m);
      y = str2num(y);

      if ((y == 0) || (y > 31))
	%% we've got the year correct
	if ((m > 12) && (d < 13))
	  %% we're operating on mm/dd/yyyy
	  tmp = m;
	  m = d;
	  d = tmp;
	elseif ((m < 13) && (d > 12))
	  %% it's fine
	else
	  %% it's ambiguous
	  error(["datesplit: ambiguous date " ds]);
	endif
      elseif ((d == 0) || (d > 31))
	%% the day and the year need to be switched
	tmp = y;
	y = d;
	d = tmp;
      else
	%% it's ambiguous
	error(["datesplit: ambiguous date " ds]);
      endif

      %% assume the beginning of the day
      h = 0;
      mi = 0;
      s = 0;
      ap = "NA";
    endif

  endif

  ## format 29  yyyymmdd                e.g. 20470313
  if (isempty(match))
    [match, y, m, d] = \
	regexp("^([0-9]{4})(1[012]|0?[0-9])(3[01]|[012][0-9])$", ds);
    %% I've never seen a date that has the year first that was not
    %% yyyy/mm/dd, so I'm going to assume that it's unambiguous.

    if (! isempty(match))
      %% assume the beginning of the day
      h = 0;
      mi = 0;
      s = 0;
      ap = "NA";
    endif
  endif

  ## format 17  QQ-YY                   e.g. Q3-00
  ## format 27  QQ-YYYY                 e.g. Q4-2132
  if (isempty(match))
    [match, q, y] = \
	regexp("^q([1-4])[-./]([0-9]{2}|[0-9]{4})$", ds);
    if (! isempty(match))
      %% Assume that it's the end of the quarter
      q = str2num(q);
      m = 3*q;
      dayopts = [31 30 30 31];
      d = dayopts(q);
    
      %% assume the end of the day
      h = 23;
      mi = 59;
      s = 59;
      ap = "NA";
    endif
  endif

  ## format 28  mmmyyyy                 e.g. Mar2047
  ## format 12  mmmyy                   e.g. Sep00
  if (isempty(match))
    [match, m, y] = \
	regexp("^([a-z]{3})([0-9]{2}|[0-9]{4})$", ds);
    if (! isempty(match))
      %% assume the beginning of the month
      d = 1;
      h = 0;
      mi = 0;
      s = 0;
      ap = "NA";
    endif
  endif

  ## format  6  mm/dd                   e.g. 09/07
  ## format 19  dd/mm                   e.g. 13/03
  if (isempty(match))
    [match, m, d] = \
	regexp("^(3[01]|[012]?[0-9])[-./](3[01]|[012][0-9])$", ds);

    if (! isempty(match))
      m = str2num(m);
      d = str2num(d);

      %% we have to determine if the date is unambiguous
      if ((m > 12) && (d < 13))
	%% we're operating on mm/dd/yyyy
	tmp = m;
	m = d;
	d = tmp;
      elseif ((m < 13) && (d > 12))
	%% it's fine
      else
	%% it's ambiguous
	error(["datesplit: ambiguous date " ds]);
      endif
      %% assume this year and the beginning of the day
      y = today(1);
      h = 0;
      mi = 0;
      s = 0;
      ap = "NA";
    endif
  endif

  ## format 10  yyyy                    e.g. 2000
  if (isempty(match))
    [match, y] = regexp("^([0-9]{4})$", ds);

    if (! isempty(match))
      %% assume the beginning of the year
      m = 1;
      d = 1;
      h = 0;
      mi = 0;
      s = 0;
      ap = "NA";
    endif
  endif

  ## format  3  mmm                     e.g. Sep
  if (isempty(match))
    m = strmatch(ds, tolower(__month_names));

    if (! isempty(m))
      match = 1;
      %% assuming the beginning of the month, this year
      y = today(1);
      d = 1;
      h = 0;
      mi = 0;
      s = 0;
      ap = "NA";
    endif
  endif

  ## format  8  ddd                     e.g. Thu
  %% People shouldn't use this function for something like this

  if (isempty(match))
    %% you mean I did all that work, and you still didn't use a valid
    %% date?  Darn you!
    error(["datesplit: Unknown date format " ds]);
  endif

  if (! isempty(match))
    if isempty(s)
      s = 0;
    elseif (ischar(s) && (1 == findstr(s,":")))
      s = s(2:3);
    endif
    if isempty(ap)
      ap = "NA";
    endif
  endif

  %% start converting the date from characters to numbers
  if (ischar(y))
    y = str2num(y);
    if (isempty(y))
      error(["datesplit: Invalid year specification " y]);
    endif
  endif
  %% Handle Y2K issues...
  if (y < 100)
    y = y + 1900;
    if (y < P)
      y = y + 100;
    endif
  endif

  if (ischar(m))
    m_num = str2num(m);
    if (isempty(m_num))
      m = strmatch(m, tolower(__month_names));
    else
      m = m_num;
    endif
    if (isempty(m))
      error(["datesplit: Invalid month specification"]);
    endif
  endif

  if (ischar(d))
    d = str2num(d);
    if (isempty(d))
      error(["datesplit: Invalid day specification " d]);
    endif
  endif

  if (ischar(h))
    h = str2num(h);
    if (isempty(h))
      error(["datesplit: Invalid hour specification " h]);
    elseif ((ap(2) == "M") && (h > 12))
      error(["datesplit: Invalid hour specification, AM or PM specified but"
	     "hour is greater than 12."]);
    endif
    
    if (strcmpi(ap, "PM") && (h < 12))
      h = h + 12;
    elseif (strcmpi(ap, "AM") && (h == 12))
      h = 0;
    endif
  endif

  if (ischar(mi))
    mi = str2num(mi);
    if (isempty(mi) || (mi > 59))
      error(["datesplit: Invalid minute specification " mi]);
    endif
  endif

  if (ischar(s))
    s = str2num(s);
    if (isempty(s) || (s > 59))
      error(["datesplit: Invalid second specification " s]);
    endif
  endif

  if (nargout <= 1)
    y = [y, m, d, h, mi, s];
  endif

endfunction

%!shared nowvec
%! nowvec=datevec(now); % Some tests could fail around midnight!
%!assert (datevec("07-Sep-2000 15:38:09"),[2000,9,7,15,38,9]);
%!assert (datevec("07-Sep-2000"),[2000,9,7,0,0,0]);
%!#ambiguous assert (datevec("09/07/00"),[2000,9,7,0,0,0]);
%!assert (datevec("Sep"),[nowvec(1),9,1,0,0,0]);
%!#ambiguous assert (datevec("09/13"),[nowvec(1),9,13,0,0,0]);
%!assert (datevec("2000"),[2000,1,1,0,0,0]);
%!assert (datevec("Sep00"),[2000,9,1,0,0,0]);
%!assert (datevec("15:38:09"),[nowvec(1:3),15,38,9]);
%!assert (datevec("03:38:09 PM"),[nowvec(1:3),15,38,9]);
%!assert (datevec("15:38"),[nowvec(1:3),15,38,0]);
%!assert (datevec("3:38 PM"),[nowvec(1:3),15,38,0]);
%!assert (datevec("Q3-00"),[2000,9,30,23,59,59]);
%!assert (datevec("Mar.03.1962 13:53:06"),[1962,3,3,13,53,6]);
%!assert (datevec("03/13/1962"),[1962,3,13,0,0,0]);
%!assert (datevec("1995/03/13"),[1995,3,13,0,0,0]);
%!assert (datevec("Q4-2132"),[2132,12,31,23,59,59]);
%!assert (datevec("Mar2047"),[2047,3,1,0,0,0]);
%!assert (datevec("20470313"),[2047,3,13,0,0,0]);
%!assert (datevec("20470313T132603"),[2047,3,13,13,26,3]);
%!assert (datevec("1047-03-13 13:26:03"),[1047,3,13,13,26,3]);
