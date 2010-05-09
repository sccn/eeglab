% plotproj() - plot projections of one or more ICA components along with 
%              the original data (returns the data plotted)
%
% Usage:
%   >> [projdata] = plotproj(data,weights,compnums);
%   >> [projdata] = plotproj(data,weights,compnums, ...
%                                 title,limits,chanlist,channames,colors);
%
% Inputs:
%   data        = single epoch of runica() input data (chans,frames) 
%   weights     = unmixing matrix (=weights*sphere)
%   compnums    = vector of component numbers to project and plot 
%
% Optional inputs:
%   title       = 'fairly short plot title' {0 -> 'plotproj()'}
%   limits      = [xmin xmax ymin ymax]  (x's in msec) 
%                          {0, or both y's 0 -> data limits}
%   chanlist    = list of data channels to plot {0 -> all}
%   channames   = channel location file or structure (see readlocs())
%   colors      = file of color codes, 3 chars per line  ('.' = space)
%                          {0 -> default color order (black/white first)}
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 05-01-96 
%
% See also: plotdata()

% Without color arg, reads filename for PROJCOLORS from icadefs.m

% Copyright (C) 05-01-96 from plotdata() Scott Makeig, SCCN/INC/UCSD,
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

% 05-25-96 added chanlist, nargin tests, rearranged variable order -sm
% 07-29-96 debugged, added chanlist channames option -sm
% 10-26-96 added test for column of compnums  -sm
% 02-18-97 improved usage message -sm
% 02-19-97 merged versions -sm
% 03-19-97 changed var() to diag(cov()), use datamean arg instead of frames/baseframes -sm
% 04-24-97 tested datamean for 1-epoch; replaced cov() with mean-squares -sm
% 05-20-97 read icadefs for PROJCOLORS & MAXPLOTDATACHANS -sm
% 06-05-97 use arbitrary chanlist as default channames -sm 
% 06-07-97 changed order of args to conform to runica -sm
% 06-12-97 made sumdata(chanlist in line 159 below -sm
% 07-23-97 removed datamean from args; let mean distribut4e over components -sm
% 09-09-97 corrected write out line " summing " -sm
% 10-31-97 removed errcode var -sm
% 11-05-97 added test for channames when chanlist ~= 1:length(chanlist) -sm & ch
% 12-19-00 adjusted new icaproj() args -sm
% 01-12-01 removed sphere arg -sm
% 01-25-02 reformated help & license, added links -ad 

function [projdata] = plotproj(data,weights,compnums,titl,limits,chanlist,channels,colors);

icadefs       % read default PROJCOLORS & MAXPLOTDATACHANS variables from icadefs.m
DEFAULT_TITLE = 'plotproj()';

%
% Substitute for missing arguments
%

if nargin < 8,
    colors = 'white1st.col';
elseif colors==0,
    colors = 'white1st.col';
end

if nargin < 7,
    channels = 0;
end
if nargin < 6
    chanlist = 0;
end
if nargin < 5,
    limits = 0;
end
if nargin < 4,
    titl = 0;
end
if titl==0,
    titl = DEFAULT_TITLE;
end

if nargin < 3,
    fprintf('plotproj(): must have at least four arguments.\n\n');
    help plotproj
    return
end
%
% Test data size
%
[chans,framestot] = size(data);

frames = framestot; % assume one epoch

[wr,wc]           = size(weights);

if wc ~= chans
    fprintf('plotproj(): sizes of weights and data incompatible.\n\n');
    return
end
%
% Substitute for 0 arguments
%
if chanlist == 0,
    chanlist = [1:chans];
end;
if compnums == 0,
    compnums = [1:wr];
end;
if size(compnums,1)>1,        % handle column of compnums !
    compnums = compnums';
end;
if length(compnums) > 256,
    fprintf('plotproj(): cannot plot more than %d channels of data at once.\n',256);
    return
end;

if channels ~= 0  % if chan name file given
   if ~all(chanlist == [1:length(chanlist)])
     fprintf('plotproj(): Cannot read an arbitrary chanlist of channel names.\n');
     return
   end
end
if channels==0,
    channels = chanlist;
end

if max(compnums)>wr,
   fprintf(...
 '\n    plotproj(): Component index (%d) > number of components (%d).\n', ...
                                 max(compnums),wr);
   return
end

fprintf('Reconstructing (%d chan, %d frame) data summing %d components.\n', ...
                   chans,frames,length(compnums));
%
% Compute projected data for single components
%
projdata = data(chanlist,:);
fprintf('plotproj(): Projecting component(s) ');
for s=compnums,            % for each component 
   fprintf('%d ',s);
   proj = icaproj(data,weights,s);   % let offsets distribute 
   projdata = [projdata proj(chanlist,:)];  % append projected data onto projdata
   % size(projdata) = [length(chanlist)  framestot*(length(compnums)+1)]
end;
fprintf('\n');
%     
% Compute percentage of variance accounted for
%
sumdata = icaproj(data,weights,compnums);% let offsets distribute 
sigmssq = mean(sum(data(chanlist,:).*data(chanlist,:)));
   data(chanlist,:) = data(chanlist,:) - sumdata(chanlist,:);
difmssq = mean(sum(data(chanlist,:).*data(chanlist,:)));

pvaf   = round(100.0*(1.0-difmssq/sigmssq)); % percent variance accounted for

rtitl = ['(p.v.a.f. ' int2str(pvaf) '%)'];
%
% Make the plot
%
plotdata(projdata,length(data),limits,titl,channels,colors,rtitl);
                                                % make the plot
