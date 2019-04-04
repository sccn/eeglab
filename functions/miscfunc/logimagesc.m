% logimagesc() - make an imagesc(0) plot with log y-axis values (ala semilogy())
%
% Usage:  >> [logfreqs,dataout] = logimagesc(times,freqs,data);
%
% Input:
%   times = vector of x-axis values
%   freqs = vector of y-axis values
%   data  = matrix of size (freqs,times)
%
% Optional Input:
%   plot = ['on'|'off'] plot image or return output (default 'on').
%
% Note: Entering text() onto the image requires specifying (x,log(y)).

% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 4/2000 

% Copyright (C) 4/2000 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 08-07-00 made ydir normal -sm
% 01-25-02 reformated help & license -ad 

function [lgfreqs,datout, h, yt, yl] = logimagesc(times,freqs,data,varargin)

  if nargin < 1
      help logimagesc;
      return
  end
  if size(data,1) ~= length(freqs)
      fprintf('logfreq(): data matrix must have %d rows!\n',length(freqs));
      datout = data;
      return
  end
  if size(data,2) ~= length(times)
      fprintf('logfreq(): data matrix must have %d columns!\n',length(times));
      datout = data;
      return
  end
  if min(freqs)<= 0
      fprintf('logfreq(): frequencies must be > 0!\n');
      datout = data;
      return
  end
  
  try, icadefs; catch, warning('Using MATLAB default colormap'); end
  
  lfreqs = log(freqs);
  lgfreqs = linspace(lfreqs(1),lfreqs(end),length(lfreqs));
  lgfreqs = lgfreqs(:);
  lfreqs = lfreqs(:);
  [mesht meshf] = meshgrid(times,lfreqs);
  try
      datout = griddata(mesht,meshf,double(data),times,lgfreqs);
  catch
      fprintf('error in logimagesc.m calling griddata.m, trying v4 method.');
      datout = griddata(mesht,meshf,data,times,lgfreqs,'v4');
  end
  datout(find(isnan(datout(:)))) = 0;
  
  if ~isempty(varargin)
      plot = varargin{2};
  else
      plot = 'on';
  end
  
  if strcmp(plot, 'on')
      imagesc(times,freqs,data);
      try colormap(DEFAULT_COLORMAP); catch, end
      nt = ceil(min(freqs)); % new tick - round up min y to int
      ht = floor(max(freqs)); % high freq - round down

      yt=get(gca,'ytick');
      yl=get(gca,'yticklabel');
      
      h=imagesc(times,lgfreqs,datout); % plot the image
      set(gca,'ydir','normal')

      i = 0; yt = [];
      yl = cell(1,100);

      tickscale = 1.618; % log scaling power for frequency ticks
      while (nt*tickscale^i < ht )
        yt = [yt log(round(nt*tickscale^i))];
        yl{i+1}=int2str(round(nt*tickscale^i));
        i=i+1;
      end

      if ht/(nt*tickscale^(i-1)) > 1.35
         yt = [yt log(ht)];
         yl{i+1} = ht;
      else
         i=i-1;
      end
     yl = {yl{1:i+1}};
     set(gca,'ytick',yt);
     set(gca,'yticklabel',yl);

%     if nt > min(yt),
%         set(gca,'ytick',log([nt yt]));
%         set(gca,'yticklabel',{int2str(nt) yl});
%      else
%         set(gca,'ytick',log([yt]));
%         set(gca,'yticklabel',{yl});
%      end

  end 
