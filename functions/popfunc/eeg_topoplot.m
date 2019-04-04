% eeg_topoplot() - plot scalp map
%
% eeg_topoplot( vals, chanlocs, 'key', 'val');
%
% Input:
%   vals     - values, one per channel
%   chanlocs - channel structure, same size as vals
%
% Optional inputs:
%   'colormap'   - colormap. Possible colormaps are 'blueredyellow', ...
%                'yellowredblue', 'bluered' or any Matlab colormap ('cool',
%                'jet', 'hsv', ...). It can also be a text file 'xxx.txt'. 
%                The text file must contain 3 columns and idealy 64 rows 
%                defining the colors in RGB format.
%   'maplimits'  - can be [min max]. This help defines the color scale for
%                maps.
%   'electrodes' - can be 'on' to show electrode dots, 'off', or 
%               'labels' to show electrode labels. Default is 'on'
%   'dotsize'    - size of electrode dots. Default is 5.
%   'shading'    - 'flat','interp'  {default: 'interp'}
%   'exclude'    - labels or indices of electrodes not to be plotted. From the
%                compiled files, these must be entered using underscores
%                for separators (e.g., "cz_pz").
%   'sphspline'  - can be 'on' or 'off'. If 'on' spherical splines are used
%                for interpolation of the scalp map. If 'off' standard 
%                planar inverse distance interpolation is used.
%   'shrink'     - shrink electrode positions (default is 0.75 to be able to
%                plot electrode at the head limit if spherical interpolation
%                is set and 0.95 for planar 2-D interpolation).
%
% References for spline interpolation:
%   [1] Perrin, F., Pernier, J., Bertrand, O., & Echallier, J. F.
%       (1989). Spherical splines for scalp potential and current
%       density mapping. Electroencephalography and Clinical
%       Neurophysiology, 72, 184-187
%   [2] Ferree, T. C. (2000). Spline Interpolation of the Scalp EEG.
%       Retrieved March 26, 2006, from
%       www.egi.com/Technotes/SplineInterpolation.pdf
%
% limitation: does not plot anything below the upper part of the head

% Copyright (C) Arnaud Delorme, SCCN, INC, 2010
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

function eeg_topoplot(values, chanlocs, varargin);

g = [];
for index = 1:2:length(varargin)
    g = setfield(g, varargin{index}, varargin{index+1});
end
if ~isfield(g, 'electrodes'), g.electrodes = 'on'; end
if ~isfield(g, 'colormap'),   g.colormap   = jet;  end
if ~isfield(g, 'maplimits'),  g.maplimits  = [];   end
if ~isfield(g, 'headrad'),    g.headrad    = [];   end
if ~isfield(g, 'sphspline'),  g.sphspline  = 'on'; end
if ~isfield(g, 'shading'),    g.shading    = 'interp'; end
if ~isfield(g, 'contour'),    g.contour    = 'on'; end
if ~isfield(g, 'dotsize'),    g.dotsize    = 5;    end
if ~isfield(g, 'mark'),       g.mark       = [];   end
if ~isfield(g, 'exclude'),    g.exclude    = [];   end
if ~isfield(g, 'linewidth'),  g.linewidth  = 2;    end
if ~isfield(g, 'shrink'),     g.shrink     = 1;    end
if ischar(g.dotsize), g.dotsize = str2num(g.dotsize); end
if any(values == 0)
    inds = find(values == 0);
    if ~isempty( [ chanlocs(inds).theta ])
        g.contour = 'off';
        g.sphspline = 'off';
    end
end

% exclude electrodes
% ------------------
if ~isempty(g.exclude)
    chanlocs(g.exclude) = [];
    values(g.exclude)   = [];
end

% find channel coordinates
% ------------------------
emptyvals = cellfun('isempty', { chanlocs.theta }); 
th = [ chanlocs.theta ];
rd = [ chanlocs.radius ];
[y x] = pol2cart(th/180*pi, rd); x=-x;
x = x*g.shrink;
y = y*g.shrink;
newvalues            = values;
newvalues(emptyvals) = [];
labls = { chanlocs.labels }; 
labls(emptyvals) = [];

if strcmpi(g.sphspline, 'on')
    % default head radius
    % -------------------
    g.headrad = 0.5;
    
    % spherical plotting
    % ------------------
    xelec = [ chanlocs.X ];
	yelec = [ chanlocs.Y ];
	zelec = [ chanlocs.Z ];

    dist = sqrt(xelec.^2+yelec.^2+zelec.^2);
	xelec = xelec./dist;
	yelec = yelec./dist;
	zelec = zelec./dist;
    
    if g.shrink ~= 1
        [th phi rad] = cart2sph(xelec, yelec, zelec);
        phi = (phi-pi/2)*g.shrink+pi/2;
        [xelec, yelec, zelec] = sph2cart(th, phi, rad);
    end;        
    
	[xsph, ysph, zsph, valsph] = spheric_spline(xelec,yelec,zelec,newvalues); 
    surf(-ysph/2,xsph/2,zsph/2,double(valsph), 'edgecolor', 'none'); view([0 0 1]);hold on;
    shading(g.shading);
    top = max(abs(valsph(:)))*1000;
    
    if strcmpi(g.contour, 'on')
    	[c h] = contour3(-ysph/2, xsph/2, valsph+top/10, 5); view([0 0 1]);
        set(h, 'cdata', [], 'edgecolor', 'k')
    end
    
	% coordinates for electrodes
	% --------------------------
    xelec(find(zelec < 0)) = [];
    yelec(find(zelec < 0)) = [];
    x = yelec/2;
    y = xelec/2;
    
else
    % default head radius
    % -------------------
    if isempty(g.headrad);
        g.headrad = max(sqrt(x.^2+y.^2));
    end
    
    % data points for 2-D data plot
    % -----------------------------
    pnts = linspace(0,2*pi,200/0.25*(g.headrad.^2));
    xx = sin(pnts)*g.headrad;
    yy = cos(pnts)*g.headrad;

	% make grid and add circle
	% ------------------------
    gridres = 30;
	coords = linspace(-g.headrad, g.headrad, gridres);
	ay = repmat(coords,  [gridres 1]);
	ax = repmat(coords', [1 gridres]);
	for ind=1:length(xx)
        [tmp closex] = min(abs(xx(ind)-coords));
        [tmp closey] = min(abs(yy(ind)-coords));
        ax(closex,closey) = xx(ind);
        ay(closex,closey) = yy(ind);
	end
	xx2 = sin(pnts)*(g.headrad-0.01);
	yy2 = cos(pnts)*(g.headrad-0.01);
	for ind=1:length(xx)
        [tmp closex] = min(abs(xx2(ind)-coords));
        [tmp closey] = min(abs(yy2(ind)-coords));
        ax(closex,closey) = xx(ind);
        ay(closex,closey) = yy(ind);
	end
	
	% linear interpolation and removal of values outside circle
	% ---------------------------------------------------------
    a = griddata(x, y, newvalues, -ay, ax, 'v4');
	aradius = sqrt(ax.^2 + ay.^2);
	indoutcircle = find(aradius(:) > g.headrad+0.01);
	a(indoutcircle) = NaN;
	surf(ay, ax, a, 'edgecolor', 'none'); view([0 0 1]); hold on;
    shading(g.shading);
    top = max(values)*1.5;

	% plot level lines
	% ----------------
    if strcmpi(g.contour, 'on')
        [c h] = contour3(ay, ax, a, 5);
        set(h, 'cdata', [], 'edgecolor', 'k')
    end
end

% plot electrodes as dots
% -----------------------
if strcmpi(g.electrodes, 'on') || strcmpi(g.electrodes, 'labels')
    rad = sqrt(x.^2 + y.^2);
    x(find(rad > g.headrad)) = [];
    y(find(rad > g.headrad)) = [];
    plot3( -x, y, ones(size(x))*top, 'k.', 'markersize', g.dotsize);
    for i = g.mark,      plot3( -x(i), y(i), double(top), 'y.', 'markersize', 4*g.dotsize); plot3( -x(i), y(i), double(top), 'r.', 'markersize', 2*g.dotsize); end
    if strcmpi(g.electrodes, 'labels')
        for index = 1:length(x)
            text( -x(index)+0.02, y(index), double(top), labls{index});
        end
    end
else
    % invisible electrode that avoid plotting problem (no surface, only
    % contours)
    plot3( -x, y, -ones(size(x))*top, 'k.', 'markersize', 0.001); 
end

% plot dipoles if any
% -------------------
if ~isempty(g.dipole)  
    hold on;
    for index = 1:size(g.dipole,1)
        g.dipole(index,:)   = g.dipole(index,:)*0.5;
        g.dipole(index,3:5) = g.dipole(index,3:5)/norm(g.dipole(index,3:end))*0.2;
        if ~any(g.dipole(index,:))
            fprintf('Note: dipole contains 0 - not plotted\n')
        elseif sum(g.dipole(index,3:4).^2) <= 0.00001 
            fprintf('Note: dipole is length 0 - not plotted\n')
        elseif sum(g.dipole(index,1:2).^2) > g.headrad
            fprintf('Note: dipole is outside plotting area - not plotted\n')
        else
            hh = plot3( -g.dipole(index, 2), g.dipole(index, 1), top, '.');
            set(hh, 'color', 'k', 'markersize', 30);
            hh = line( -[g.dipole(index, 2) g.dipole(index, 2)+g.dipole(index, 4)]', ...
                [g.dipole(index, 1) g.dipole(index, 1)+g.dipole(index, 3)]',[top top]);
            set(hh, 'color', 'k', 'linewidth', 30/7);
        end
    end
end

% special colormaps
% -----------------
if ischar(g.colormap) 
    if ~isempty(strmatch(g.colormap, { 'hsv' 'jet' 'gray' 'hot' 'cool' 'bone' ...
            'copper', 'pink' 'flag' 'prism' }, 'exact'))
    else % read text file
        g.colormap = load('-ascii', g.colormap);
    end
end;    
colormap(g.colormap);

if ~isempty(g.maplimits)
    if ~ischar(g.maplimits) && ~isempty(g.maplimits) && ~isnan(g.maplimits(1))
        caxis(g.maplimits);
    end
end

% main circle
% -----------
radiuscircle = 0.5;
pnts   = linspace(0,2*pi,200);
xc     = sin(pnts)*radiuscircle;
yc     = cos(pnts)*radiuscircle;
sf     = 1; % scaling factor
plot3(xc*sf,yc*sf,ones(size(xc))*top, 'k', 'linewidth', g.linewidth); hold on;

% ears & nose
% -----------
rmax  = 0.5;
base  = rmax-.0046;
basex = 0.18*rmax;                   % nose width
tip   = 1.15*rmax; 
tiphw = .04*rmax;                    % nose tip half width
tipr  = .01*rmax;                    % nose tip rounding
q = .04; % ear lengthening
EarX  = [.497-.005  .510  .518  .5299 .5419  .54    .547   .532   .510   .489-.005]; % rmax = 0.5
EarY  = [q+.0555 q+.0775 q+.0783 q+.0746 q+.0555 -.0055 -.0932 -.1313 -.1384 -.1199];

plot3(EarX*sf,EarY*sf,ones(size(EarX))*top,'color','k','LineWidth',g.linewidth)    % plot left ear
plot3(-EarX*sf,EarY*sf,ones(size(EarY))*top,'color','k','LineWidth',g.linewidth)   % plot right ear
plot3([basex;tiphw;0;-tiphw;-basex]*sf,[base;tip-tipr;tip;tip-tipr;base]*sf,top*ones(size([basex;tiphw;0;-tiphw;-basex])),'color','k','LineWidth',g.linewidth);

% axis limits
% -----------
axis off;
set(gca, 'ydir', 'normal');
axis equal
ylimtmp = max(g.headrad, 0.58);
ylim([-ylimtmp ylimtmp]);

% ----------------
% spherical spline
% ----------------
function [x, y, z, Res] = spheric_spline( xelec, yelec, zelec, values);

SPHERERES = 40;
[x,y,z] = sphere(SPHERERES);
x(1:(length(x)-1)/2,:) = [];
y(1:(length(x)-1)/2,:) = [];
z(1:(length(x)-1)/2,:) = [];

Gelec = computeg(xelec,yelec,zelec,xelec,yelec,zelec);
Gsph  = computeg(x,y,z,xelec,yelec,zelec);

% equations are 
% Gelec*C + C0  = Potential (C unknow)
% Sum(c_i) = 0
% so 
%             [c_1]
%      *      [c_2]
%             [c_ ]
%    xelec    [c_n]
% [x x x x x]         [potential_1]
% [x x x x x]         [potential_ ]
% [x x x x x]       = [potential_ ]
% [x x x x x]         [potential_4]
% [1 1 1 1 1]         [0]

% compute solution for parameters C
% ---------------------------------
meanvalues = mean(values); 
values = values - meanvalues; % make mean zero
C = pinv([Gelec;ones(1,length(Gelec))]) * [values(:);0];

% apply results
% -------------
Res = zeros(1,size(Gsph,1));
for j = 1:size(Gsph,1)
    Res(j) = sum(C .* Gsph(j,:)');
end
Res = Res + meanvalues;
Res = reshape(Res, size(x));

% compute G function
% ------------------
function g = computeg(x,y,z,xelec,yelec,zelec)

unitmat = ones(length(x(:)),length(xelec));
EI = unitmat - ((repmat(x(:),1,length(xelec)) - repmat(xelec,length(x(:)),1)).^2 +... 
                (repmat(y(:),1,length(xelec)) - repmat(yelec,length(x(:)),1)).^2 +...
                (repmat(z(:),1,length(xelec)) - repmat(zelec,length(x(:)),1)).^2)/2;

g = zeros(length(x(:)),length(xelec));
m = 4; % 3 is linear, 4 is best according to Perrin's curve
for n = 1:7
    L = legendre(n,EI);
    g = g + ((2*n+1)/(n^m*(n+1)^m))*squeeze(L(1,:,:));
end
g = g/(4*pi);    

% find electrode indices
% ----------------------
function allinds = elecind( str, chanlocs, values );

    findmax = 0;
    findmin = 0;
    if ~iscell(str)
         if strmatch(str, 'max', 'exact'), findmax = 1; end
         if strmatch(str, 'min', 'exact'), findmin = 1; end;         
         indunderscore = [ 0 find( str == '_' ) length(str)+1 ];
    else indunderscore = [1:length(str)+1];
    end
     
    % find maximum or minimum
    % -----------------------
    if findmax, [tmp allinds] = max(values); return; end
    if findmin, [tmp allinds] = min(values); return; end
    
    % find indices for labels
    % -----------------------
    labels = lower({ chanlocs.labels });
    for i = 1:length(indunderscore)-1
        if ~iscell(str)
             tmpstr = str(indunderscore(i)+1:indunderscore(i+1)-1);
        else tmpstr = str{i};
        end
        tmpind = strmatch(lower(tmpstr), labels, 'exact');
        if isempty(tmpind)
            if str2num(tmpstr) > 0
                tmpind = str2num(tmpstr);
            else
                error(sprintf('Could not find channel "%s"', tmpstr));
            end
        end
        allinds(i) = tmpind;
    end
    
        
