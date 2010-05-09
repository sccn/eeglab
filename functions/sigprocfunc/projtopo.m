% projtopo() - plot projections of one or more ICA components along with 
%              the original data in a 2-d topographic array. Returns 
%              the data plotted. Click on subplot to examine separately.
% Usage:
%       >> [projdata] = projtopo(data,weights,[compnums],'chan_locs',...
%                                 'title',[limits],colors,chans);
% Inputs:
%   data       = single epoch of runica() input data (chans,frames) 
%   weights    = unimxing weight matrix (runica() weights*sphere)
%  [compnums]  = vector of component numbers to project and plot 
%  'chan_locs' = channel locations file. Example: >> topoplot example  
%                 Else [rows cols] for rectangular grid array
% Optional:
%  'title'     = (short) plot title {default|0 -> none}
%  [limits]    = [xmin xmax ymin ymax]  (x's in msec) 
%                          {default|0|both y's 0 -> use data limits}
%   colors     = file of color codes, 3 chars per line  ('.' = space)
%                          {default|0 -> default color order (black/white first)}
%   chans      = vector of channel numbers to plot {default|0: all}
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 04-02-98 
%
% See also: plotproj(), topoplot()

% Copyright (C) 04-02-98 from plotproj() Scott Makeig, SCCN/INC/UCSD,
% scott@sccn.ucsd.edu
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

% Without color arg, reads filename for PROJCOLORS from icadefs.m

% 03-15-00 added plotchans arg -sm
% 03-16-00 added axcopy() feature -sm & tpj
% 12-19-00 adjusted new icaproj() args -sm
% 12-20-00 removed argument "sphere" -sm
% 07-12-01 fixed nargin test to allow 8 args -sm & cmk
% 01-25-02 reformated help & license, added links -ad 

function [projdata] = projtopo(data,weights,compnums,chan_locs,titl,limits,colors,plotchans)

icadefs      % read default PROJCOLORS & MAXPLOTDATACHANS variables from icadefs.m
axsize = 0;  % use plottopo() default
DEFAULT_TITLE = '';
%
% Substitute for missing arguments
%
if nargin < 7,
    colors = 'white1st.col';
elseif colors==0 | isempty(colors)
    colors = 'white1st.col';
end

if nargin < 6,
    limits = 0;
end
if nargin < 5,
    titl = 0;
end
if titl==0,
    titl = DEFAULT_TITLE;
end

if nargin < 4 | nargin > 8
    help projtopo
    fprintf('projtopo(): requires 4-8 arguments.\n\n');
    return
end
%
% Test data size
%
[chans,framestot] = size(data);
if ~exist('plotchans') | isempty(plotchans) | plotchans==0
   plotchans = 1:chans; % default
end
frames = framestot; % assume one epoch

[wr,wc]           = size(weights);
%
% Substitute for 0 arguments
%
if compnums == 0,
    compnums = [1:wr];
end;
if size(compnums,1)>1,        % handle column of compnums !
    compnums = compnums';
end;
if length(compnums) > MAXPLOTDATACHANS,
    fprintf(...
  'projtopo(): cannot plot more than %d channels of data at once.\n',...
         MAXPLOTDATACHANS);
    return
end;

if max(compnums)>wr,
   fprintf(...
 '\n    projtopo(): Component index (%d) > number of components (%d).\n', ...
                                 max(compnums),wr);
   return
end

fprintf('Reconstructing (%d chan, %d frame) data summing %d components.\n', ...
                   chans,frames,length(compnums));
%
% Compute projected data for single components
%
projdata = data;
fprintf('projtopo(): Projecting component(s) ');
for s=compnums,            % for each component 
   fprintf('%d ',s);
   proj = icaproj(data,weights,s); % let offsets distribute 
   projdata = [projdata proj];  % append projected data onto projdata
end;
fprintf('\n');
%
% Make the plot
%
% >> plottopo(data,'chan_locs',frames,limits,title,channels,axsize,colors,ydir) 

plottopo(projdata,chan_locs,size(data,2),limits,titl,plotchans,axsize,colors);
                                                % make the plottopo() plot
axcopy(gcf);
