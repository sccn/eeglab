% gradplot() - Compute the gradient of EEG scalp map(s) on a square grid
%
% Usage:
%            >> [gradX, gradY] = gradplot(maps,eloc_file,draw)
% Inputs:
%    maps      - Activity levels, size (nelectrodes,nmaps)
%    eloc_file - Electrode location filename (.loc file) containing electrode
%              - coordinates (For format, see >> topoplot example)
%    draw      - If defined, draw the gradient map(s) {default: no}
%
% Outputs:
%   gradX  - Gradients in X direction
%   gradY  - Gradients in Y directions 
%
% Note: Use cart2pol() to convert to polar (amp, direction) coordinates).
%
% Authors: Marissa Westerfield & Arnaud Delorme, CNL/Salk Institute, La Jolla 3/10/01
%
% See also: topoplot(), lapplot()

% Copyright (C) 3/10/01 Marissa Westerfield & Arnaud Delorme, CNL/Salk Institute, La Jolla
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

% 01-25-02 reformated help & license, added links -ad 

function [gradx, grady] = gradplot( map, locs_file, draw ) 

if nargin < 2
	help gradplot;
	return;
end

NCHANS = size(map,1);
GRID_SCALE = 2*NCHANS+5;
MAX_RADIUS = 0.5;

% --------------------------------
% Read the electrode location file
% --------------------------------
if ischar(locs_file) % a locs file
        [tmpeloc labels Th Rd ind] = readlocs(locs_file,'filetype', ...
                                              'loc');
        [x,y] = pol2cart(Th/180*pi,Rd); % See Bug 149
	% [x,y] = pol2cart(Th,Rd);
elseif isstruct(locs_file)  % a locs struct
        [tmpeloc labels Th Rd ind] = readlocs(locs_file);
        if max(abs(Rd))>0.5
          fprintf('gradplot(): Shrinking electrode arc_lengths from (max) %4.3f to (max) 0.5\n',...
                          max(abs(Rd)));
          Rd = Rd/(2*max(abs(Rd))); % shrink to max radius = 0.5
        end
	[x,y] = pol2cart(Th,Rd);
        if length(x) ~= NCHANS
         fprintf(...
          'gradplot(): channels in locs file (%d) ~= channels in maps (%d)\n', ...
                                            length(x),              NCHANS);
         return
       end
else                % what is this complex number format ??? -sm
	x = real(locs_file);
	y = imag(locs_file);
end;	

% ------------------------------------------------
% Locate nearest position of electrode in the grid 
% ------------------------------------------------
xi = linspace(-0.5,0.5,GRID_SCALE);   % x-axis description (row vector)
yi = linspace(-0.5,0.5,GRID_SCALE);   % y-axis description (row vector)
delta = xi(2)-xi(1); % length of grid entry
horizidx=zeros(1, NCHANS); %preallocation
vertidx=zeros(1, NCHANS); % preallocation
    
for c=1:NCHANS
   [useless_var horizidx(c)] = min(abs(y(c) - xi)); % find pointers to electrode
   [useless_var  vertidx(c)] = min(abs(x(c) - yi));  % positions in Zi
end
   
% -------------------
% Compute gradient(s)
% -------------------
for m=1:size(map,2) 
   	[Xi,Yi,Zi] = griddata(y,x,map(:,m),yi',xi, 'v4'); % interpolate data
   	[FX,FY] = gradient(Zi);
	positions = horizidx + (vertidx-1)*GRID_SCALE;

	gradx(:,m) = FX(positions(:));
	grady(:,m) = FY(positions(:));

	% ----------------
	% Draw gradient(s)
	% ----------------
	if exist('draw')

        % Take data within head
        mask = (sqrt(Xi.^2+Yi.^2) <= MAX_RADIUS); 
        mask = find(mask==0);
        Zi(mask) = NaN;
        FX(mask) = NaN;
        FY(mask) = NaN;
        width = max(Xi)-min(Xi);

        subplot(ceil(sqrt(size(map,2))), ceil(sqrt(size(map,2))), m);

%        surface(1+width*(0.5+Xi-delta/2),...
%                1+width*(0.5+Yi-delta/2),...
%                    zeros(size(Zi)),Zi,'EdgeColor','none',...
%                    'FaceColor','flat'); hold on

         contour(imresize(Zi,0.5)); hold on
         quiver(imresize(FX, 0.5), imresize(FY, 0.5)); 
        title(['Map ' int2str(m)]);


% %%% Draw Head %%%%
ax = axis; 
width = ax(2)-ax(1);
expansion = 0.3;
axis([ax(1)-expansion ax(2)+expansion ax(3)-expansion ax(4)+expansion])
steps = 0:2*pi/100:2*pi;
basex = .18*MAX_RADIUS;  
tip = MAX_RADIUS*1.15; 
base = MAX_RADIUS-.004;
EarX = [.497 .510 .518 .5299 .5419 .54 .547 .532 .510 .489];
EarY = [.0555 .0775 .0783 .0746 .0555 -.0055 -.0932 -.1313 -.1384 -.1199];

HCOLOR = 'k';
HLINEWIDTH = 1.8;

% Plot Head, Ears, Nose
hold on
plot(1+width/2+cos(steps).*MAX_RADIUS*width,...
     1+width/2+sin(steps).*MAX_RADIUS*width,...
    'color',HCOLOR,'Linestyle','-','LineWidth',HLINEWIDTH); % head

plot(1+width/2+[.18*MAX_RADIUS*width;0;-.18*MAX_RADIUS*width],...
     1+width/2+[base;tip;base]*width,...
    'Color',HCOLOR,'LineWidth',HLINEWIDTH);                 % nose
   
plot(1+width/2+EarX*width,...
     1+width/2+EarY*width,...
           'color',HCOLOR,'LineWidth',HLINEWIDTH)           % l ear
plot(1+width/2-EarX*width,...
     1+width/2+EarY*width,...
           'color',HCOLOR,'LineWidth',HLINEWIDTH)           % r ear

hold off
axis off


	end
end

return;
