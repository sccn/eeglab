% imagescloglog() - make an imagesc(0) plot with log y-axis and
%                   x-axis values
%
% Usage:  >> imagescloglog(times,freqs,data);
% Usage:  >> imagescloglog(times,freqs,data,clim,xticks,yticks,'key','val',...);
%
% Inputs:
%   times = vector of x-axis values (LOG spaced)
%   freqs = vector of y-axis values (LOG spaced)
%   data  = matrix of size (freqs,times)
%
% Optional inputs:
%   clim   = optional color limit
%   xticks = graduation for x axis
%   yticks = graduation for y axis
%   ...    = 'key', 'val' properties for figure
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 4/2003 

% Copyright (C) 4/2003 Arnaud Delorme, SCCN/INC/UCSD, scott@sccn.ucsd.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function imagescloglog(times,freqs,data,clim, xticks, yticks, varargin)

  if size(data,1) ~= length(freqs)
      fprintf('logfreq(): data matrix must have %d rows!\n',length(freqs));
      return
  end
  if size(data,2) ~= length(times)
      fprintf('logfreq(): data matrix must have %d columns!\n',length(times));
      return
  end
  if min(freqs)<= 0
      fprintf('logfreq(): frequencies must be > 0!\n');
      return
  end
  try, icadefs; catch, warning('Using MATLAB default colormap'); end
  
  steplog = log(times(2))-log(times(1)); % same for all points
  realborders = [exp(log(times(1))-steplog/2) exp(log(times(end))+steplog/2)];
  newtimes    = linspace(realborders(1), realborders(2), length(times));
  
  % regressing 3 times
  border  = mean(newtimes(2:end)-newtimes(1:end-1))/2; % automatically added to the borders in imagesc
  newtimes = linspace(realborders(1)+border, realborders(2)-border, length(times));
  border  = mean(newtimes(2:end)-newtimes(1:end-1))/2; % automatically added to the borders in imagesc
  newtimes = linspace(realborders(1)+border, realborders(2)-border, length(times));
  border  = mean(newtimes(2:end)-newtimes(1:end-1))/2; % automatically added to the borders in imagesc
  newtimes = linspace(realborders(1)+border, realborders(2)-border, length(times));

  % problem with log images in Matlab: border are automatically added
  % to account for half of the width of a line: but they are added as
  % if the data was linear. The commands below compensate for this effect
  
  steplog = log(freqs(2))-log(freqs(1)); % same for all points
  realborders = [exp(log(freqs(1))-steplog/2) exp(log(freqs(end))+steplog/2)];
  newfreqs    = linspace(realborders(1), realborders(2), length(freqs));
  
  % regressing 3 times
  border  = mean(newfreqs(2:end)-newfreqs(1:end-1))/2; % automatically added to the borders in imagesc
  newfreqs = linspace(realborders(1)+border, realborders(2)-border, length(freqs));
  border  = mean(newfreqs(2:end)-newfreqs(1:end-1))/2; % automatically added to the borders in imagesc
  newfreqs = linspace(realborders(1)+border, realborders(2)-border, length(freqs));
  border  = mean(newfreqs(2:end)-newfreqs(1:end-1))/2; % automatically added to the borders in imagesc
  newfreqs = linspace(realborders(1)+border, realborders(2)-border, length(freqs));
  
  if nargin == 4 && ~isempty(clim)
      imagesc(newtimes,newfreqs,data,clim);
  else 
      imagesc(newtimes,newfreqs,data);
  end
  
  set(gca, 'yscale', 'log', 'xscale', 'log');
  try colormap(DEFAULT_COLORMAP); catch, end
  
  % puting ticks
  % ------------
  if nargin >= 5
      divs = xticks;
  else 
      divs = linspace(log(times(1)), log(times(end)), 10);
      divs = ceil(exp(divs)); divs = unique_bc(divs); % ceil is critical here, round might misalign
                                               % out-of border label with within border ticks
  end
  set(gca, 'xtickmode', 'manual');
  set(gca, 'xtick', divs);
  if nargin >= 6
      divs = yticks;
  else 
      divs = linspace(log(freqs(1)), log(freqs(end)), 10);
      divs = ceil(exp(divs)); divs = unique_bc(divs); % ceil is critical here, round might misalign
                                               % out-of border label with within border ticks
  end
  set(gca, 'ytickmode', 'manual');
  set(gca, 'ytick', divs);
  
  % additional properties
  % ---------------------
  set(gca, 'yminortick', 'off', 'xaxislocation', 'bottom', 'box', 'off', 'ticklength', [0.03 0], 'tickdir','out', 'color', 'none');  
  if ~isempty(varargin)
      set(gca, varargin{:});
  end
  
