% movav() - Perform a moving average of data indexed by xvals.
%           Supports use of a moving non-rectangular window.
%           Can be used to resample a data matrix.
% Usage:
%     >> [outdata,outx] = movav(data,xvals,xwidth,xadv,firstx,lastx,xwin);
%
% Inputs:
%   data   = input data (chans,frames)
%   xvals  = index for each data frame (column) {def|0 -> 1:frames}
%            Note that default is fastest, assumes equal x-spacing.
%   xwidth = smoothing-window width in xvals {def|0 -> (lastx-firstx)/4}
%   xadv   = xvals step size {default|0 -> 1}
%            Note that to reduce yyy frames to xxx, use yyy/(xxx+2)
%   firstx = low xval of first averaging window {def|0 -> low xvals}
%   lastx  = high xval of last averaging window {def|0 -> high xvals}
%   xwin   = vector of window values {def|0 -> ones() = square window}
%            May be long; linear interp. is NOT used between values.
%            An example is >> gauss(1001,2) ->  [0.018 ... 1.0 ... 0.018]
%
% Outputs:
%   outdata = smoothed data (chans,
%   outx    = xval midpoints of successive output data windows
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 10-25-97 

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 10-25-97 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: not supported by cvs2svn $

% 3-20-98 fixed bug in multi-channel windowed averaging -sm
% 6-10-98 changed mean() and sum() to nanmean() and nansum() -sm
% 2-16-99 tested for stat toolbox functions nanmean() and nansum() -sm
% 9-03-01 fixed gauss() example -sm
% 01-25-02 reformated help & licenses -ad 

function [outdata,outx] = movav(data,xvals,xwidth,xadv,firstx,lastx,xwin)

MAXPRINT = 1; % max outframe numbers to print on tty
NEARZERO = 1e-22;
verbose = 0;

nanexist = 0;  
if exist('nanmean') % test for stat toolbox nan routines
   nanexist = 1;
end
if nargin<1
   help movav
   return
else
  [chans,frames]=size(data);
end
if chans>1 & frames == 1,
  data   = data';   % make row vector
  tmp    = chans;
  chans  = frames;
  frames = tmp;
end
if frames < 4
  fprintf('movav(): data size (%d,%d) is too short.]\n',chans,frames);
  return
end

fastave = 0;
if nargin<2 
  xvals = 0;
end
if length(xvals)==1 
  if size(xvals,1)>1
    xvals=xvals';    % make row vector
  elseif xvals(1)==0,
    fastave =1;
  else
    help movav
    return
  end
end
if size(xvals,1)>1 & size(xvals,2)>1
  help movav
  return
end
if fastave == 0 & frames ~= length(xvals)
    fprintf('movav(): columns in (%d) xvals vector and (%d) in data matrix must be equal.\n',length(xvals),size(data,2));
    return
end

if nargin < 7
  xwin = 0;
end

if nargin < 6
  lastx = 0;
end
if lastx == 0,
  if fastave
    lastx = frames;
  else
    lastx = max(xvals);
  end
end

if nargin<5,
  firstx = 0;
end
if firstx==0,
  if fastave
    firstx = 1;
  else
    firstx = min(xvals);
  end
end

if nargin<4,
  xadv = 0;
end
if isempty(xadv) | xadv == 0,
  xadv = 1.0;
end

if nargin<3,
  xwidth = 0;
end
if xwidth==0,
  xwidth = (lastx-firstx)/4;
end
wlen = 1;  % default;
if fastave==0
  if length(xwin)==1 & xwin ~=0,  % should be a vector or 0
    help movav
    return
  elseif size(xwin,1)>1 & size(xwin,2)>1 % not a matrix
    help movav
    return
  end
  if size(xwin,1)>1
    xwin = xwin';   % make row vector
  end

  if xwin~=0
    if abs(sum(xwin)) < NEARZERO
      fprintf('movav(): abs(sum(xwin)) too small. Cannot normalize.\n');
    else
      xwin = xwin/abs(sum(xwin)); % make xwin values sum to 1;
    end
    wlen = length(xwin);
  end
end

outframes = floor(((lastx-firstx+xadv)-xwidth)/xadv)+1;
if verbose
  fprintf('movav() will output %d frames.\n',outframes);
end
if outframes < 1,
   outframes = 1;
end
outdata = zeros(chans,outframes);
outx = zeros(1,outframes);
outxval = firstx+xwidth/2;
%
%%%%%%%%%%%%%%%%%%%%%% Print header %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if verbose,
  fprintf('Performing moving averaging:\n')
  fprintf('Output will be %d chans by %d frames',chans,outframes);
  if wlen>1,
    fprintf(' using the specified %d-point window.\n',wlen);
  else
    fprintf(' using a square window.\n');
  end
end
 %fprintf('   firstx = %g, lastx= %g, xwidth = %g xadv = %g\n',...
 %                 firstx,lastx,xwidth,xadv);
%
%%%%%%%%%%%%%%%%%%% Perform averaging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
lox = firstx;
i = 0; %fastave default
for f=1:outframes
    hix = lox+xwidth;
    outx(1,f)=outxval;
    outxval = outxval + xadv;
    if fastave==0
      i = find(xvals>=lox & xvals <= hix);
    end
    if length(i)==0,
      if f>1,
       outdata(:,f) = outdata(:,f-1); % If no data, replicate
if verbose,
 fprintf('r');
end
      else
       outdata(:,f) = zeros(chans,1); %  or else output zeros
if verbose
  fprintf('0');
end
      end
    elseif length(xwin)==1,
      if fastave
        if nanexist
         outdata(:,f) = nanmean(data(:,round(lox):round(hix))')'; % Else average
        else
         outdata(:,f) = mean(data(:,round(lox):round(hix))')'; % Else average
        end
      else
        if nanexist
          outdata(:,f) = nanmean(data(:,i)')'; % Else average
        else
          outdata(:,f) = mean(data(:,i)')'; % Else average
        end
      end
if verbose
  fprintf('.');
end
%
%%%%%%%%%%%%%%%%% Windowed averaging %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    else                                 
       wadv=(hix-lox)/wlen;
       ix = ceil((xvals(i)-lox)/wadv);
       zs = find(ix==0);
       ix(zs) = ones(1,zs);
       if length(xwin)>1
          sumx = sum(xwin(ix));
       else
          sumx=1;
       end
       if abs(sumx) < NEARZERO  % cannot normalize
         if f>1,
          outdata(:,f) = outdata(:,f-1); % if no data, replicate
if verbose,
  fprintf('R');
end
         else
          outdata(:,f) = zeros(chans,1); % or output zeros
if verbose,
  fprintf('0');
end
         end
       else
        if nanexist
         outdata(:,f) = nansum((((ones(chans,1)*xwin(ix)).*data(:,i))/sumx)')'; 
        else                % perform normalized windowed smoothing
         outdata(:,f) = sum((((ones(chans,1)*xwin(ix)).*data(:,i))/sumx)')'; 
        end               
       end 
    end
    lox = lox+xadv;
    if (outframes<MAXPRINT) 
      fprintf('%d ',f);
    end
end
if verbose,
  fprintf('\n');
end

