% logimagesc() - make an imagesc(0) plot with log y-axis values (ala semilogy())
%
% Usage:  >> [logfreqs,dataout] = logimagesc(times,freqs,data);
%
% Input:
%   times = vector of x-axis values
%   freqs = vector of y-axis values
%   data  = matrix of size (freqs,times)
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 4/2000 

% Copyright (C) 4/2000 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% 08-07-00 made ydir normal -sm
% 01-25-02 reformated help & license -ad 

function [lgfreqs,datout] = logimagesc(times,freqs,data)

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
  lfreqs = log(freqs);
  lgfreqs = linspace(lfreqs(1),lfreqs(end),length(lfreqs));
  lgfreqs = lgfreqs(:);
  lfreqs = lfreqs(:);
  [mesht meshf] = meshgrid(times,lfreqs);
  datout = griddata(mesht,meshf,data,times,lgfreqs);

  imagesc(times,freqs,data);
  nt = ceil(min(freqs)); % new tick - round up min freq to int-hz
  yt=get(gca,'ytick');
  yl=get(gca,'yticklabel');

  imagesc(times,lgfreqs,datout);
  set(gca,'ydir','normal')
  set(gca,'ytick',log([nt yt]));
  set(gca,'yticklabel',{int2str(nt) yl});
