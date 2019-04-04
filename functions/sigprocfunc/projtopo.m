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
elseif colors==0 || isempty(colors)
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

if nargin < 4 || nargin > 8
    help projtopo
    fprintf('projtopo(): requires 4-8 arguments.\n\n');
    return
end
%
% Test data size
%
[chans,framestot] = size(data);
if ~exist('plotchans') || isempty(plotchans) || plotchans==0
   plotchans = 1:chans; % default
end
frames = framestot; % assume one epoch

[wr,wc]           = size(weights);
%
% Substitute for 0 arguments
%
if compnums == 0,
    compnums = [1:wr];
end
if size(compnums,1)>1,        % handle column of compnums !
    compnums = compnums';
end
if length(compnums) > MAXPLOTDATACHANS,
    fprintf(...
  'projtopo(): cannot plot more than %d channels of data at once.\n',...
         MAXPLOTDATACHANS);
    return
end

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
end
fprintf('\n');
%
% Make the plot
%
% >> plottopo(data,'chan_locs',frames,limits,title,channels,axsize,colors,ydir) 

plottopo(projdata,chan_locs,size(data,2),limits,titl,plotchans,axsize,colors);
                                                % make the plottopo() plot
axcopy(gcf);
