% eventlock() - DEPRECATED: Please use eegalign() instead.
% eventlock() - Time lock data epochs to specified event frames or event x-values.
%               Use to timelock existing data epochs to reaction times or other events
%
% Usage:
%  >> [dataout,medval,shiftframes] = eventlock(data,frames,eventframes,medval);
%  >> [dataout,medval,shiftframes] = eventlock(data,xvals,eventvals,medval);
%
% Inputs for multi-channel data:
%   data        = input data, size(chans,frames*epochs) 
%   frames      = scalar, frames per epoch {0 -> data length}
%   eventframes = time locking event frames, size(1,epochs)
%   medval      = median eventframe to align to {default: median(eventvals)}
%
% Inputs for single-channel data:
%   data         = input data, size(frames, epochs)
%   xvals        = vector of epoch time-values, size(1,frames) 
%                  OR xvals = [startval nframes srate] (ms N Hz)
%   eventvals    = x-values of time locking events, size(1,epochs)
%   medval       = median event time to align to {default: median(eventvals)}
%
% Outputs:
%   dataout      = shifted/sorted data out; size(dataout) = size(data)
%   medval       = median eventval (time or frame). Data is aligned to this.
%   shifts       = number of frames shifted for each epoch
%                 (Note: shift > 0 means shift forward|right)
%
% Note: Missing values are filled with NaN. Some toolbox functions
%       ( timef(), crossf(), crosscoher() ) handle NaNs correctly. 
%       To truncate to non-NaN values, use
%       >> dataout = matsel(dataout,frames,...
%              1+max(shifts(find(shifts>=0))):frames+min(shifts(find(shifts<=0))));
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 8/20/99

% Copyright (C) 8/20/99 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 8/23/99 debugged right shifts -sm
% 8/25/99 allowed single- or multiple-channel data -sm
% 3/9/00 fixed order of dims in help msg -sm & ev
% 01-25-02 reformated help & license -ad 

function [data,medval,shifts] = eventlock(data,arg2,eventvals,medval)

VERBOSE=0;
%
%%%%%%%%%%%%%% Reshape data to (frames,ntrials) %%%%%%%%%%%%%
%
if nargin < 3
  help eventlock
  return 
end
MULTICHANNEL = 0;
if length(arg2) == 1 % then assume multi-channel: size(chans,epochs*frames) 
   MULTICHANNEL = 1;
   if arg2==0
    frames = length(eventvals); % default 1 epoch
   else
    frames = arg2;
   end
   xvals = 1:frames;
else                 % assume single-channel: size(frames,epochs) 
   MULTICHANNEL = 0;
   xvals = arg2;
   if length(xvals) ~= 3
     frames = length(xvals);
   else
     frames = xvals(2);
   end
end
rows = size(data,1); 
if MULTICHANNEL
   ntrials = size(data,2)/frames;
else
   ntrials = size(data,2);
end

if length(xvals) == 3 % assume [startval nvals srate] format
  srate = xvals(3);
  endval = xvals(1)+(xvals(2)-1)*xvals(3)+1e-16;
  xvals = xvals(1):1000/xvals(3):endval; 
  xvals = xvals(1:frames);  % make sure of length!
else
  srate = 1000/(xvals(2)-xvals(1));
end

if MULTICHANNEL
  if frames*ntrials ~= size(data,2)
   fprintf(...
 'Input data length (%d) is not a multiple of the eventframes length (%d).\n',...
                   size(data,2),                                    frames);
   return
  end 
elseif frames~= size(data,1) 
   fprintf(...
 'Input data frames (%d) is not equal to xvals length (%d).\n',...
                   size(data,1),                                    length(xvals));
   return
end 

if length(eventvals) ~= ntrials
   fprintf(...
 'eventlock(): Number of eventvals (%d) must equal number of trials (%d).\n',...
                     length(eventvals), ntrials);
   return
end 
  
if MULTICHANNEL
  fprintf('Input data is %d channels by %d frames * %d trials.\n',...
                             rows,frames,ntrials);
else
  fprintf('Input data is %d frames by %d trials.\n',...
                             frames,ntrials);
end

%
%%%%%%%%%%%%%%%%%%% Align data to eventvals %%%%%%%%%%%%%%%%%%%
%

if MULTICHANNEL
    frames = rows*frames;
    data = reshape(data,frames,ntrials);
end % NB: data is now size((chans*)frames,ntrials)

aligndata=repmat(nan,frames,ntrials); % begin with matrix of NaN's
shifts = zeros(1,ntrials);

if ~exist('medval') | isempty(medval) 
   medval= median(eventvals);
end
[medval medx] = min(abs(medval-xvals)); 
 medval = xvals(medx);

if MULTICHANNEL
  fprintf('Aligning data to median event frame %g.\n',medval); 
else % SINGLE CHANNEL
  fprintf('Aligning data to median event time %g.\n',medval); 
end

for i=1:ntrials, %%%%%%%%% foreach trial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% idx = eventvals(i);
% shft  = round(medx-idx);
 [tmpx idx] = min(abs(xvals-eventvals(i)));
 shft  = medx-idx;

if VERBOSE;fprintf('Trial %d -shift = %d\n',i,shft);end

 shifts(i) = shft;
 if MULTICHANNEL
   shft = shft*rows; % shift in multiples of chans
 end

 if shft>0, % shift right %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if frames > shft 
     aligndata((1+shft):end,i) = data(1:(end-shft),i);
  else
     fprintf('No eventlocked data for epoch %d - required shift too large!\n',i);
  end

 elseif shft < 0 % shift left %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if frames > -shft
     aligndata(1:(end+shft),i) = data((1-shft):end,i);
  else
     fprintf('No eventlocked data for epoch %d - required shift too large.\n',i);
  end
 else % shft == 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     aligndata(:,i) = data(:,i);
 end 
end % end trial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Shifted trials by %d to %d frames.\n',min(shifts),max(shifts));
if ~MULTICHANNEL
   fprintf('                  %g to %g msec.\n',...
                        1000/srate*min(shifts),1000/srate*max(shifts));
end
data = aligndata;                       % now data is aligned to sortvar

if MULTICHANNEL
  data = reshape(data,rows,frames/rows*ntrials); % demultiplex channels
end
