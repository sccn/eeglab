## Copyright (C) 2000 Paul Kienzle
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
## @deftypefn {Function File} {} datestr(date,code,P)
## Format the given date/time according to the format @code{code}.  The date
## 730736.65149 (2000-09-07 15:38:09.0934) would be formated as follows:
## @multitable @columnfractions 0.1 0.45 0.45
## @item @strong{Code} @tab @strong{Format} @tab @strong{Example}
## @item  0 @tab dd-mmm-yyyy HH:MM:SS   @tab 07-Sep-2000 15:38:09
## @item  1 @tab dd-mmm-yyyy            @tab 07-Sep-2000 
## @item  2 @tab mm/dd/yy               @tab 09/07/00 
## @item  3 @tab mmm                    @tab Sep 
## @item  4 @tab m                      @tab S 
## @item  5 @tab mm                     @tab 09
## @item  6 @tab mm/dd                  @tab 09/07 
## @item  7 @tab dd                     @tab 07 
## @item  8 @tab ddd                    @tab Thu 
## @item  9 @tab d                      @tab T 
## @item 10 @tab yyyy                   @tab 2000 
## @item 11 @tab yy                     @tab 00
## @item 12 @tab mmmyy                  @tab Sep00 
## @item 13 @tab HH:MM:SS               @tab 15:38:09 
## @item 14 @tab HH:MM:SS PM            @tab 03:38:09 PM
## @item 15 @tab HH:MM                  @tab 15:38 
## @item 16 @tab HH:MM PM               @tab 03:38 PM 
## @item 17 @tab QQ-YY                  @tab Q3-00
## @item 18 @tab QQ                     @tab Q3
## @item 19 @tab dd/mm                  @tab 13/03
## @item 20 @tab dd/mm/yy               @tab 13/03/95
## @item 21 @tab mmm.dd.yyyy HH:MM:SS   @tab Mar.03.1962 13:53:06
## @item 22 @tab mmm.dd.yyyy            @tab Mar.03.1962
## @item 23 @tab mm/dd/yyyy             @tab 03/13/1962
## @item 24 @tab dd/mm/yyyy             @tab 12/03/1962
## @item 25 @tab yy/mm/dd               @tab 95/03/13
## @item 26 @tab yyyy/mm/dd             @tab 1995/03/13
## @item 27 @tab QQ-YYYY                @tab Q4-2132
## @item 28 @tab mmmyyyy                @tab Mar2047
## @item 29 @tab yyyymmdd               @tab 20470313
## @item 30 @tab yyyymmddTHHMMSS        @tab 20470313T132603
## @item 31 @tab yyyy-mm-dd HH:MM:SS    @tab 1047-03-13 13:26:03
## @end multitable
##
## If no code is given or code is -1, then use 0, 1 or 13 as the
## default, depending on whether the date portion or the time portion 
## of the date is empty.
##
## If a vector of dates is given, a vector of date strings is returned.
##
## The parameter @code{P} is needed by @code{datevec} to convert date strings
## with 2 digit years into dates with 4 digit years.  See @code{datevec}
## for more information.
##
## @seealso{date,clock,now,datestr,datenum,calendar,weekday} 
## @end deftypefn

## TODO: parse arbitrary code strings.
## TODO: e.g., for  Wednesday 2001-03-05 09:04:06 AM, use
## TODO:     yy    01
## TODO:     yyyy  2001
## TODO:     m     M
## TODO:     mm    03
## TODO:     mmm   Mar
## TODO:     d     W
## TODO:     dd    05
## TODO:     ddd   Wed
## TODO:     HH    09
## TODO:     MM    04
## TODO:     SS    06
## TODO:     PM    AM
## TODO: Vectorize.  It is particularly easy since all the codes are
## TODO:    fixed width.  Just generate the parts in separate arrays and
## TODO:    concatenate.
## TODO: Invent format codes which don't have leading zeros or uppercase?

function retval = datestr(date,code,P)
  if (nargin == 0 || nargin > 3 )
    usage("datestr(date [, code]) or datestr('date' [, code [, P]])");
  endif
  if (nargin < 3) P = []; endif
  if (nargin < 2) code = []; endif
  V = datevec(date, P);

  if (isempty(code))
    if (all(V(:,1)==0) && all(all(V(:,2:3) == 1)))
      code = 13;
    elseif (all(V(:,4:6)==0))
      code = 1; 
    else
      code = 0;
    endif
  endif

  global __month_names = ["Jan";"Feb";"Mar";"Apr";"May";"Jun";...
			  "Jul";"Aug";"Sep";"Oct";"Nov";"Dec"];
  global __time_names = ["AM";"PM"];
  for i=1:rows(V)
    Y=V(i,1); M=V(i,2); D=V(i,3);
    h=V(i,4); m=V(i,5); s=V(i,6);
    Y2 = rem(Y,100);
    switch (code)
      case { 0, 'dd-mmm-yyyy HH:MM:SS' }
	str = sprintf("%02d-%s-%04d %02d:%02d:%02d",...
			    D,__month_names(M,:),Y,h,m,floor(s));
      case { 1, 'dd-mmm-yyyy' }
	str = sprintf("%02d-%s-%04d",D,__month_names(M,:),Y);
      case { 2, 'mm/dd/yy' }
	str = sprintf("%02d/%02d/%02d",M,D,Y2);
      case { 3, 'mmm' } 
	str = sprintf("%s",__month_names(M,:));
      case { 4, 'm' }
	str = sprintf("%s",__month_names(M,1));
      case { 5, 'mm' } 
	str = sprintf("%02d",M);
      case { 6, 'mm/dd' } 
	str = sprintf("%02d/%02d",M,D);
      case { 7, 'dd' } 
	str = sprintf("%02d",D);
      case { 8, 'ddd' } 
	[d,str] = weekday(datenum(Y,M,D));
      case { 9, 'd' } 
	[d,str] = weekday(datenum(Y,M,D));
	str = str(1);
      case { 10, 'yyyy' } 
	str = sprintf("%04d", Y);
      case { 11, 'yy' } 
	str = sprintf("%02d", Y2);
      case { 12, 'mmmyy' } 
	str = sprintf("%s%02d", __month_names(M,:),Y2);
      case { 13, 'HH:MM:SS' } 
	str = sprintf("%02d:%02d:%02d", h, m, floor(s));
      case { 14, 'HH:MM:SS PM' }
	str = sprintf("%02d:%02d:%02d %s", rem(h,12), m, floor(s), \
			     __time_names(floor(h/12)+1,:));
      case { 15, 'HH:MM' } 
	str = sprintf("%02d:%02d", h, m);
      case { 16, 'HH:MM PM' }
	str = sprintf("%02d:%02d %s", rem(h,12), m, \
			     __time_names(floor(h/12)+1,:));
      case { 17, 'QQ-YY' } 
	str = sprintf("Q%d-%02d", floor((M+2)/3),Y2);
      case { 18, 'QQ' }
	str = sprintf("Q%d", floor((M+2)/3));
      case { 19, 'dd/mm' }
	str = sprintf('%02d/%02d', D, M);
      case { 20, 'dd/mm/yy' }
	str = sprintf('%02d/%02d/%02d', D, M, Y2);
      case { 21, 'mmm.dd.yyyy HH:MM:SS' }
	str = sprintf('%s.%02d.%04d %02d:%02d:%02d', 
		__month_names(M,:), D, Y, h, m, s);
      case { 22, 'mmm.dd.yyyy' }
	str = sprintf('%s.%02d.%04d', __month_names(M,:), D, Y);
      case { 23, 'mm/dd/yyyy' }
	str = sprintf('%02d/%02d/%04d', M, D, Y);
      case { 24, 'dd/mm/yyyy' }
	str = sprintf('%02d/%02d/%04d', D, M, Y);
      case { 25, 'yy/mm/dd' }
	str = sprintf('%02d/%02d/%02d', Y2, M, D);
      case { 26, 'yyyy/mm/dd' }
	str = sprintf('%04d/%02d/%02d', Y, M, D);
      case { 27, 'QQ-YYYY' }
	str = sprintf('Q%d-%04d', floor((M+2)/3), Y);
      case { 28, 'mmmyyyy' }
	str = sprintf('%s%04d', __month_names(M,:), Y);
      case { 29, 'yyyymmdd' }
	str = sprintf('%04d%02d%02d', Y, M, D);
      case { 30, 'yyyymmddTHHMMSS' }
	str = sprintf('%04d%02d%02dT%02d%02d%02d',Y,M,D,h,m,s);
      case { 31, 'yyyy-mm-dd HH:MM:SS' }
	str = sprintf('%04d-%02d-%02d %02d:%02d:%02d',Y,M,D,h,m,s);
    endswitch
    if i == 1
      retval = str;
    else 
      retval = [ retval ; str ] ;
    endif
  endfor
endfunction
