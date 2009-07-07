%% -*- texinfo -*-
%% @deftypefn {Function File} {} datenum(Y, M, D [, h , m [, s]])
%% @deftypefnx {Function File} {} datenum('date' [, P])
%% Returns the specified local time as a day number, with Jan 1, 0000
%% being day 1. By this reckoning, Jan 1, 1970 is day number 719529.  
%% The fractional portion, corresponds to the portion of the specified day.
%%
%% Years can be negative and/or fractional.
%% Months below 1 are considered to be January.
%% Days of the month start at 1.
%% Days beyond the end of the month go into subsequent months.
%% Days before the beginning of the month go to the previous month.
%% Days can be fractional.
%%
%% XXX WARNING XXX this function does not attempt to handle Julian
%% calendars so dates before Octave 15, 1582 are wrong by as much
%% as eleven days.  Also be aware that only Roman Catholic countries
%% adopted the calendar in 1582.  It took until 1924 for it to be 
%% adopted everywhere.  See the Wikipedia entry on the Gregorian 
%% calendar for more details.
%%
%% XXX WARNING XXX leap seconds are ignored.  A table of leap seconds
%% is available on the Wikipedia entry for leap seconds.
%%
%% @seealso{date,clock,now,datestr,datevec,calendar,weekday}
%% @end deftypefn

%% Algorithm: Peter Baum (http://vsg.cape.com/~pbaum/date/date0.htm)
%% Author: Paul Kienzle
%% This program is granted to the public domain.

function n = datenum(Y,M,D,h,m,s)
  monthstart = [306,337,0,31,61,92,122,153,184,214,245,275];

  if nargin == 0 || (nargin > 2  && isstr(Y)) || nargin > 6
    usage('n=datenum(''date'' [, P]) or n=datenum(Y, M, D [, h, m [, s]])');
  end
  if ischar(Y)
    if nargin < 2, M=[]; end
    [Y,M,D,h,m,s] = datevec(Y,M);
  else
    if nargin < 6, s = 0; end
    if nargin < 5, m = 0; end
    if nargin < 4, h = 0; end
    if nargin == 1
      nc = columns(Y);
      if nc > 6 || nc < 3,
        error('expected date vector containing [Y,M,D,h,m,s]');
      end
      s = 0;
      m = 0; 
      h = 0;
      if nc >= 6, s = Y(:,6); end
      if nc >= 5, m = Y(:,5); end
      if nc >= 4, h = Y(:,4); end
      D = Y(:,3);
      M = Y(:,2);
      Y = Y(:,1);
    end 
  end

  M(M<1) = 1; %% For compatibility.  Otherwise allow negative months.

  %% Set start of year to March by moving Jan. and Feb. to previous year.
  %% Correct for months > 12 by moving to subsequent years.
  z = Y + fix((M-14)/12);

  %% Lookup number of days to beginning of the month.
  M = mod(M-1,12)+1;
  f = M;
  f(:) = monthstart(M);

  %% Add Y+M+D+h+m+s and correct for leap years.
  n = D + f + 365*z+floor(z/4)-floor(z/100)+floor(z/400) + 60 + ...
	(h+(m+s/60)/60)/24;


%!assert(datevec(datenum(2003,11,28)),[2003,11,28,0,0,0])
