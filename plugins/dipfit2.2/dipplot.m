% dipplot() - Visualize EEG equivalent-dipole locations and orientations 
%             in the BESA spherical head model or in a averaged MRI model.
% Usage:
%   >> dipplot( sources, 'key', 'val', ...);
%   >> [sources X Y Z XE YE ZE] = dipplot( sources, 'key', 'val', ...);
%
% Inputs:
%   sources   - structure array of dipole information:
%                besaexent: BESA eccentricity of the dipole
%                besathloc: BESA azimuth angle of the dipole
%                besaphloc: BESA horizontal angle of the dipole
%                besathori: BESA azimuth angle of the dipole orientation
%                besaphori: BESA horiz. angle of the dipole orientation
%                optional fields: component: component number
%                rv:              residual variance
%
% Optional input:
%  'color'    - [cell array of color strings or (1,3) color arrays]. For
%               exemple { 'b' 'g' [1 0 0] } gives blue, green and red. 
%               Dipole colors will rotate through the given colors if
%               the number given is less than the number of dipoles to plot.
%  'view'     - 3-D viewing angle (see >> help view). In cartesian coords.,
%               [0 0 1] gives a sagittal view, [0 -1 0] a view from the rear;
%               [1 0 0] gives a view from the side of the head.
%  'mesh'     - ['on'|'off'] Display spherical mesh. {Default is 'on'}
%  'gui'      - ['on'|'off'] Display controls. {Default is 'on'} If gui 'off', 
%               a new figure is not created. Useful for incomporating a dipplot 
%               into a complex figure.
%  'num'      - ['on'|'off'] Display component number. Take into account
%               dipole size. {Default: 'off'}
%  'summary'  - Build a summary plot with three views (top, back, side)
%  'image'    - ['besa'|'mri'] Background image. {Default: 'besa'} 'mri' uses
%               mean-MRI brain images from the Montreal Neurological Institute.
%  'rvrange'  - [min max] Only plot dipoles with residual variace within the
%               given range. Default: plot all dipoles.
%  'projimg'  - ['on'|'off'] Project dipole(s) onto the 2-D images, for use
%               in making 3-D plots {Default 'off'}
%  'dipolesize' - Size of the dipole sphere(s) {Default: 30}
%  'dipolelength' - Length of the dipole bar(s) {Default: 30}
%  'pointout' - ['on'|'off'] Point the dipoles outward. {Default: 'off'}
%
% Outputs:
%   sources   - EEG.source structure with updated 'X', 'Y' and 'Z' fields
%   X,Y,Z     - Locations of dipole heads (Cartesian coordinates)
%   XE,YE,ZE  - Locations of dipole ends (Cartesian coordinates)
%
% Author: Arnaud Delorme, CNL / Salk Institute, 1st July 2002
%
% Notes: Visualized locations are not exactly the same as in BESA (because of
% manual tuning of the size of the head textures).  Because of a bug in the 
% Matlab warp() function, the side-view texture cannot be displayed. 
% To diplay head textures, the files 'besarear.pcx' and 'besasagittal.pcx' 
% are required. The Matlab image processing toolbox is also required.
%
% Example:
%  % position and orientation of the first dipole
%  sources(1).besaexent= 69.036;
%  sources(1).besathloc= -26.71;
%  sources(1).besaphloc= 19.702;
%  sources(1).besathori= -80.02;
%  sources(1).besaphori= 87.575;
%  % position and orientation of the second dipole
%  sources(2).besaexent= 69.036;
%  sources(2).besathloc= 46;
%  sources(2).besaphloc= 39;
%  sources(2).besathori= 150;
%  sources(2).besaphori= -3;
%  % plot of the two dipoles (first in green, second in blue)
%  dipplot( sources, 'color', { 'g' 'b' }); 
%
%  % To make a stereographic plot
%  figure; 
%  subplot(1,2,1); dipplot( sources, 'view', [43 10], 'gui', 'off');
%  subplot(1,2,2); dipplot( sources, 'view', [37 10], 'gui', 'off');
%
%  % To make a summary plot
%  dipplot( sources, 'summary', 'on', 'dipolesize', 15, 'mesh', 'off');
%
% See also: eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2002 Arnaud Delorme
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

% README -- Plotting strategy:
% - All buttons have a tag 'tmp' so they can be removed
% - The component-number buttons have 'userdata' equal to 'editor' and 
%   can be found easily by other buttons find('userdata', 'editor')
% - All dipoles have a tag 'dipolex' (x=their number) and can be made 
%     visible/invisible
% - All dipoles have also 'userdata' equal to 'dipole' and can be
%     selected easily
% - The gcf object stores the handle of the dipole that is currently
%     being modified

function [sources, x, y, z, xe, ye, ze] = dipplot( sources, varargin )
    
    DEFAULTVIEW = [0 0 1];
        
    if nargin < 1
        help dipplot;
        return;
    end;
    
    if isfield(sources, 'besathloc')
        sources = convertbesaoldformat(sources);
    end;
    if isfield(sources, 'dip')
        for index = 1:length(sources)
            sources(index).posxyz = sources(index).dip.pos/4;
            tmp = sources(index).posxyz(:,1);
            sources(index).posxyz(:,1) = sources(index).posxyz(:,2);
            sources(index).posxyz(:,2) = -tmp;
            sources(index).momxyz = sources(index).dip.mom'/10000;
            tmp = sources(index).momxyz(:,1);
            sources(index).momxyz(:,1) = sources(index).momxyz(:,2);
            sources(index).momxyz(:,2) = -tmp;
        end;
    end;
    if ~isfield(sources, 'posxyz')
        sources = computexyzforbesa(sources);
    end;        
    
    % reading and testing arguments
    % -----------------------------
    if ~isstruct(sources)
        error('dipplot: ''sources'' must be a strcuture');
    end;
    
    %                             key        type       range             default
    g = finputcheck( varargin, { 'color'    ''         []                 [];
                                 'mesh'     'string'   { 'on' 'off' }     'off';
                                 'gui'      'string'   { 'on' 'off' }     'on';
                                 'summary'  'string'   { 'on' 'off' }     'off';
                                 'view'     'real'     []                 [];
                                 'rvrange'  'real'     [0 Inf]            [];
                                 'num'      'string'   { 'on' 'off' }     'off';
                                 'projimg'  'string'   { 'on' 'off' }     'off';
                                 'pointout'  'string'   { 'on' 'off' }     'off';
                                 'dipolesize' 'real'   [0 Inf]            30;
                                 'dipolelength' 'real'   [0 Inf]          30;
                                 'image'    'string'   { 'besa' 'mri' }   'besa' }, 'dipplot');
    if isstr(g), error(g); end;
    
    % remove sources with out of bound Residual variance
    % --------------------------------------------------
    if isfield(sources, 'rv') & ~isempty(g.rvrange)
        for index = length(sources):-1:1
            if sources(index).rv < g.rvrange(1) | sources(index).rv > g.rvrange(2)
                sources(index) = [];
            end;
        end;
    end;
    
    if strcmpi(g.image, 'besa')
        IMAGESLOC   = { 'besatop.pcx' 'besarear.pcx' 'besaside.pcx' };
        IMAGESOFFSET  = { [0.0  0.07 NaN ]  [0    NaN -0.01]  [ NaN -0.01   -0.025 ] };
        IMAGESMULT    = { [1.23 1.20 NaN ]  [1.28 NaN  1.13]  [ NaN 1.15   1.15 ] };
        %IMAGESOFFSET = { [0 -0.07]   [0 -0.03]  [ -0.12 -0.02 ] };
        %IMAGESMULT   = { [1.22 1.22] [1.28 1.3] [ 1.325 1.02 ] };
        COLORMESH = [.5 .5 .5];
        BACKCOLOR = 'w';
        AXISLIM   = [-1.2 1.2 -1.2 1.2 -1.2 1.2];
    else
        IMAGESLOC   = { 'mritop.pcx' 'mrirear.pcx' 'mriside.pcx' };
        IMAGESOFFSET = { [-0.01 0.005  NaN]   [-0.02 NaN 0.11]  [NaN 0.05 0.31] } ;% [ -0.12 0] };
        IMAGESMULT   = { [1.14  1.08   NaN]   [1.13  NaN 1.04]  [NaN 1.13  0.88]} ;%[ 1.4 1.15 ] };
        COLORMESH = 'w';
        BACKCOLOR = 'k';
        %AXISLIM   = [-1.2 1.2 -1.2 1.2 -1.2 1.2];
        AXISLIM   = [-1.4 1.4 -1.1 1.1 -1.2 1.2];
    end;

    % build summarized figure
    % -----------------------
    if strcmp(g.summary, 'on')
        figure;
        options = { 'gui', 'off', 'dipolesize', g.dipolesize,'dipolelength', g.dipolelength,'color', g.color, 'mesh', g.mesh, 'num', g.num, 'image', g.image };
        axes('position', [0 0 0.5 0.5]);  dipplot(sources, 'view', [1 0 0] , options{:}); axis off;
        axes('position', [0 0.5 0.5 .5]); dipplot(sources, 'view', [0 0 1] , options{:}); axis off;
        axes('position', [.5 .5 0.5 .5]); dipplot(sources, 'view', [0 -1 0], options{:}); axis off;
        axes('position', [0.5 0 0.5 0.5]); 
        %p = get(gcf, 'position');
        %p(2) = p(2)+p(4)-800;
        %p(4) = 800;
        %p(3) = 800;
        %set(gcf, 'position', p);
        colorcount = 1;
        if isfield(sources, 'component')
            for index = 1:length(sources)
                if index~=1 
                    if sources(index).component ~= sources(index-1).component
                        textgui(colorcount) = { sprintf('Component %d (R.V. %3.2f)', sources(index).component, 100*sources(index).rv) };
                        colorcount = colorcount+1;
                    end;
                else 
                    textgui(colorcount) = { sprintf('Component %d (R.V. %3.2f)', sources(index).component, 100*sources(index).rv) };
                    colorcount = colorcount+1;
                end;
            end;
            colorcount = colorcount-1;
            allstr = strvcat(textgui{:});
            h = text(0,0.5, allstr);
            if colorcount >= 15, set(h, 'fontsize', 8);end;
            if colorcount >= 20, set(h, 'fontsize', 6);end;
            if strcmp(BACKCOLOR, 'k'), set(h, 'color', 'w'); end;
        end;
        axis off;
        return;
    end;
    
    if isempty(g.color)
        g.color = { 'g' 'b' 'r' 'm' 'c' 'y' };
        if strcmp(BACKCOLOR, 'w'), g.color = { g.color{:} 'k' }; end;
    end;
    g.color = g.color(mod(0:length(sources)-1, length(g.color)) +1);
    if ~iscell(g.color)
        error('dipplot: ''color'' must be a cell array');
    end;
    
    % internal constants
    % ------------------
    scaling = 0.0105;
    ori = 90;

    % plot head graph in 3D
    % ---------------------
    if strcmp(g.gui, 'on')
        figure;
    end;
    try,
        warning off;  a = double(imread(IMAGESLOC{1}))/255; warning on;
        if ndims(a) == 2, a(:,:,2) = a; a(:,:,3) = a(:,:,1); end;
        wx = ([-1 1; -1 1]+IMAGESOFFSET{1}(1))*IMAGESMULT{1}(1);
        wy = ([-1 -1; 1 1]+IMAGESOFFSET{1}(2))*IMAGESMULT{1}(2);
        wz = [-1 -1; -1 -1]*1.07;
        % DISPLAY IF REAL RATIO OF IMAGE CORRESPOND TO RATIO SHOWN ON SCREEN
        %fprintf('Image ratio %3.2f\tCoord ratio:%3.2f\n', size(a,2)/size(a,1),  ...
        %       IMAGESMULT{1}(1)/ IMAGESMULT{1}(2) / (AXISLIM(2)-AXISLIM(1)) * (AXISLIM(4)-AXISLIM(3)) );
        surface(wx, wy, wz, a(end:-1:1,:,:), 'FaceColor','texturemap', 'EdgeColor','none', 'CDataMapping','direct');
        hold on; %%fill3([-2 -2 2 2], [-2 2 2 -2], wz(:)-1, BACKCOLOR);
    catch, error(lasterr); end;
    try,
        warning off; a = double(imread(IMAGESLOC{2}))/255;  warning on;
        if ndims(a) == 2, a(:,:,2) = a; a(:,:,3) = a(:,:,1); end;
        wx = ([-1 1; -1 1]+IMAGESOFFSET{2}(1))*IMAGESMULT{2}(1);
        wy = [1 1; 1 1]*1.07;
        wz = ([1 1; -1 -1]+IMAGESOFFSET{2}(3))*IMAGESMULT{2}(3);
        hold on; surface(wx, wy, wz, a, 'FaceColor','texturemap', 'EdgeColor','none', 'CDataMapping','direct');
        %fprintf('Image ratio %3.2f\tCoord ratio:%3.2f\n', size(a,2)/size(a,1),  ...
        %        IMAGESMULT{2}(1)/ IMAGESMULT{2}(3) / (AXISLIM(2)-AXISLIM(1)) * (AXISLIM(6)-AXISLIM(5)));
        %%fill3([-2 -2 2 2], wy(:)-1, [-2 2 2 -2], BACKCOLOR);
    catch, error(lasterr); end;
    try,
        warning off; a = double(imread(IMAGESLOC{3}))/255;  warning on;
        if ndims(a) == 2, a(:,:,2) = a; a(:,:,3) = a(:,:,1); end;
        wx = [-1 -1; -1 -1]*1.07;
        wy = ([-1 1; -1 1]+IMAGESOFFSET{3}(2))*IMAGESMULT{3}(2);
        wz = ([1 1; -1 -1]+IMAGESOFFSET{3}(3))*IMAGESMULT{3}(3);
        hold on; surface(wx, wy, wz, a, 'FaceColor','texturemap', 'EdgeColor','none', 'CDataMapping','direct');
        %fprintf('Image ratio %3.2f\tCoord ratio:%3.2f\n', size(a,2)/size(a,1), ...
        %        IMAGESMULT{3}(2)/ IMAGESMULT{3}(3) / (AXISLIM(4)-AXISLIM(3)) * (AXISLIM(6)-AXISLIM(5)));
        %%fill3([-2 -2 2 2], wy(:)-1, [-2 2 2 -2], BACKCOLOR);
    catch, error(lasterr); end;
    set(gca, 'color', BACKCOLOR);
    %warning off; a = imread('besaside.pcx'); warning on;
    % BECAUSE OF A BUG IN THE WARP FUNCTION, THIS DOES NOT WORK (11/02)
    %hold on; warp([], wy, wz, a);
    if ~isempty(g.view)
        view(g.view);
    else
        view(DEFAULTVIEW);
    end;
    axis square;
    axis(AXISLIM);
    axis off;
    
    % plot sphere mesh and nose
    % -------------------------
    [x y z] = sphere(10);
    x = x*100*scaling; y = y*100*scaling; z=z*100*scaling;
    hold on; 
    h = line(x,y,z); set(h, 'color', COLORMESH, 'linestyle', '--', 'tag', 'mesh');
    h = line(x,z,y); set(h, 'color', COLORMESH, 'linestyle', '--', 'tag', 'mesh');
    %h = line([0 0;0 0],[-1 -1.2; -1.2 -1], [-0.3 -0.7; -0.7 -0.7]);
    %set(h, 'color', COLORMESH, 'linewidth', 3, 'tag', 'noze');
    
    % determine max length if besatextori exist
    % -----------------------------------------
    sizedip = [];
    for index = 1:length(sources)
        sizedip = [ sizedip sources(index).momxyz(3) ]; 
    end;
    maxlength = max(sizedip);
    
    for index = 1:length(sources)
        nbdip = 1;
        if size(sources(index).posxyz, 1) > 1 & any(sources(index).posxyz(2,:)) nbdip = 2; end;
        for dip = 1:nbdip
        
            x = sources(index).posxyz(dip,1);
            y = sources(index).posxyz(dip,2);
            z = sources(index).posxyz(dip,3);
            xo = sources(index).momxyz(dip,1)*g.dipolelength;
            yo = sources(index).momxyz(dip,2)*g.dipolelength;
            zo = sources(index).momxyz(dip,3)*g.dipolelength;
            
            if abs([x+xo,y+yo,z+zo]) >= abs([x,y,z])
                xo = x+xo;
                yo = y+yo;
                zo = z+zo;
            elseif strcmpi(g.pointout,'on')
                xo = x-xo; % make dipole point outward from head center
                yo = y-yo;
                zo = z-zo;
            else
                xo = x+xo;
                yo = y+yo;
                zo = z+zo;
            end
            x = -x; xo=-xo;
            y = -y; yo=-yo;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% draw dipole bar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tag = [ 'dipole' num2str(index) ];
            h = line( [x xo]', [y yo]', [z zo]');
            set(h, 'userdata', 'dipole', 'tag', tag, 'color','k', 'linewidth', g.dipolesize/7.5);
            if strcmp(BACKCOLOR, 'k'), set(h, 'color', g.color{index}); end;
            
            % draw point
            hold on;
            h = plot3(x,  y,  z); 
            set(h, 'userdata', 'dipole', 'tag', tag, ...
                   'marker', '.', 'markersize', g.dipolesize, 'color', g.color{index});
            
            % project onto images
            if strcmpi(g.projimg, 'on')
                if isstr(g.color{index})
                    switch g.color{index}
                     case 'y', g.color{index} = [1 1 0]; % yellow
                     case 'm', g.color{index} = [1 0 1];
                     case 'c', g.color{index} = [0 1 1];
                     case 'r', g.color{index} = [1 0 0];
                     case 'g', g.color{index} = [0 1 0];
                     case 'b', g.color{index} = [0 0 1];
                     case 'w', g.color{index} = [1 1 1];
                     case 'k', g.color{index} = [0 0 0];
                    end;
                end;
                tmpcolor = g.color{index} / 2;
                
                % project onto z axis
                tag = [ 'dipole' num2str(index) ];
                h = line( [x xo]', [y yo]', [-1 -1]');
                set(h, 'userdata', 'dipole', 'tag', tag, 'color','k', 'linewidth', g.dipolesize/7.5);
                if strcmp(BACKCOLOR, 'k'), set(h, 'color', tmpcolor); end;
                h = plot3(x,  y,  -1); 
                set(h, 'userdata', 'dipole', 'tag', tag, ...
                       'marker', '.', 'markersize', g.dipolesize, 'color', tmpcolor);
                
                % project onto x axis
                tag = [ 'dipole' num2str(index) ];
                h = line( [x xo]', [-1 -1]', [z zo]');
                set(h, 'userdata', 'dipole', 'tag', tag, 'color','k', 'linewidth', g.dipolesize/7.5);
                if strcmp(BACKCOLOR, 'k'), set(h, 'color', tmpcolor); end;
                h = plot3(x,  -1,  z); 
                set(h, 'userdata', 'dipole', 'tag', tag, ...
                       'marker', '.', 'markersize', g.dipolesize, 'color', tmpcolor);
            end;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% draw circle  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if isfield(sources, 'stdX') & sources(index).stdX ~= 0
                stdX = num2str(sources(index).stdX); strX = num2str(x);
                stdY = num2str(sources(index).stdY); strY = num2str(y);
                stdZ = num2str(sources(index).stdZ); strZ = num2str(z);
                if isempty(g.view) | g.view(3) ~= 0
                    h = myezplot3([strX '+cos(t)*' stdX], [strY '+sin(t)*' stdY], strZ, [0,2*pi]);
                    set(h, 'color', g.color{index}, 'tag', tag, 'userdata', 'dipole' );
                end;
                if isempty(g.view) | g.view(2) ~= 0
                    h = myezplot3([strX '+cos(t)*' stdX], strY, [strZ '+sin(t)*' stdZ], [0,2*pi]);
                    set(h, 'color', g.color{index}, 'tag', tag, 'userdata', 'dipole' );
                end;
                if isempty(g.view) | g.view(1) ~= 0
                    h = myezplot3(strX, [strY '+cos(t)*' stdY], [strZ '+sin(t)*' stdZ], [0,2*pi]);
                    set(h, 'color', g.color{index}, 'tag', tag, 'userdata', 'dipole' );
                end;
            end;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% draw text  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if isfield(sources, 'component')
                if strcmp(g.num, 'on')
                    h = text(x,  y,  z, [ '  ' int2str(sources(index).component)]);
                    set(h, 'userdata', 'dipole', 'tag', tag, 'fontsize', g.dipolesize/2 );
                end;
                if index~=length(sources) 
                    if sources(index).component ~= sources(index+1).component
                        textgui(index) = { sprintf('C %d (%3.2f)', sources(index).component, sources(index).rv*100) };
                    end;
                else 
                    textgui(index) = { sprintf('C %d (%3.2f)', sources(index).component, sources(index).rv*100) };
                end;
            else
                textgui(index) = { '' };
            end;
        end;
    end;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% buttons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h = uicontrol( 'unit', 'normalized', 'position', [0 0 .15 .05], 'tag', 'tmp', ...
                  'style', 'pushbutton', 'string', 'Top view', 'callback', 'view([0 0 1]);');
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.05 .15 .05], 'tag', 'tmp', ...
                  'style', 'pushbutton', 'string', 'Coronal view', 'callback', 'view([0 -1 0]);');
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.1 .15 .05], 'tag', 'tmp', ...
                  'style', 'pushbutton', 'string', 'Sagital view', 'callback', 'view([1 0 0]);');
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.15 .15 .05], 'tag', 'tmp', ...
                  'style', 'text', 'string', '');
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.2  .15 .05], 'tag', 'tmp', ...
                  'style', 'pushbutton', 'string', 'Mesh on', 'callback', ...
                   'set(findobj(''parent'', gca, ''tag'', ''mesh''), ''visible'', ''on'');');
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.25 .15 .05], 'tag', 'tmp', ...
                  'style', 'pushbutton', 'string', 'Mesh off', 'callback', ...
                   'set(findobj(''parent'', gca, ''tag'', ''mesh''), ''visible'', ''off'');');
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.3 .15 .05], 'tag', 'tmp', ...
                  'style', 'text', 'string', 'Display:','fontweight', 'bold' );
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.35 .15 .05], 'tag', 'tmp', ...
                  'style', 'text', 'string', '');
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.4 .15 .05], 'tag', 'tmp', ...
                  'style', 'pushbutton', 'fontweight', 'bold', 'string', 'No controls', 'callback', ...
                   'set(findobj(''parent'', gcbf, ''tag'', ''tmp''), ''visible'', ''off'');');
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.45 .15 .05], 'tag', 'tmp', ...
                  'style', 'text', 'string', '');
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.50 .15 .05], 'tag', 'tmp', ...
                  'style', 'pushbutton', 'string', 'Keep|Prev', 'callback', ...
                [ 'editobj = findobj(''parent'', gcf, ''userdata'', ''editor'');' ...
                  'set(editobj, ''string'', num2str(str2num(get(editobj, ''string''))-1));' ...
                  'tmpobj = get(gcf, ''userdata'');' ...
                  'eval(get(editobj, ''callback''));' ...
                  'set(tmpobj, ''visible'', ''on'');' ...
                  'clear editobj tmpobj;' ]);
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.55 .15 .05], 'tag', 'tmp', ...
                  'style', 'pushbutton', 'string', 'Prev', 'callback',  ...
                [ 'editobj = findobj(''parent'', gcf, ''userdata'', ''editor'');' ...
                  'set(editobj, ''string'', num2str(str2num(get(editobj, ''string''))-1));' ...
                  'eval(get(editobj, ''callback''));' ...
                  'clear editobj;' ]);
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.6 .15 .05], 'tag', 'tmp', ...
                  'style', 'pushbutton', 'string', 'Next', 'callback', ...
                [ 'editobj = findobj(''parent'', gcf, ''userdata'', ''editor'');' ...
                  'set(editobj, ''string'', num2str(str2num(get(editobj, ''string''))+1));' ...
                  'eval(get(editobj, ''callback''));' ...
                  'clear editobj;' ]);
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.65 .15 .05], 'tag', 'tmp', ...
                  'style', 'pushbutton', 'string', 'Keep|Next', 'callback', ...
                [ 'editobj = findobj(''parent'', gcf, ''userdata'', ''editor'');' ...
                  'set(editobj, ''string'', num2str(str2num(get(editobj, ''string''))+1));' ...
                  'tmpobj = get(gcf, ''userdata'');' ...
                  'eval(get(editobj, ''callback''));' ...
                  'set(tmpobj, ''visible'', ''on'');' ...
                  'clear editobj tmpobj;' ]);
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.7 .15 .05], 'tag', 'tmp', 'userdata', 'rv', ...
                  'style', 'text', 'string', '');
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.75 .15 .05], 'tag', 'tmp', 'userdata', 'editor', ...
                  'style', 'edit', 'string', '1', 'callback', ...
                [ 'editobj = findobj(''parent'', gcf, ''userdata'', ''editor'');' ...
                  'tmpnum = str2num(get(editobj, ''string''));' ...
                  'if tmpnum < 1, tmpnum = 1; end;' ...
                  'if tmpnum >' num2str(length(sources)) ', tmpnum = ' num2str(length(sources)) '; end;' ...
                  'set(editobj, ''string'', num2str(tmpnum));' ...
                  'set(get(gcf, ''userdata''), ''visible'', ''off'');' ...
                  'newdip = findobj(''parent'', gca, ''tag'',' ...
                       '[''dipole'' get(editobj, ''string'')]);' ...
                  'set(newdip, ''visible'', ''on'');' ...
                  'set(gcf, ''userdata'', newdip);' ...
                  'tmprvobj = findobj(''parent'', gcf, ''userdata'', ''rv'');' ...
                  'tmpallrv = get( gca, ''userdata'');' ...
                  'set( tmprvobj, ''string'', tmpallrv{tmpnum});' ...
                  'clear newdip editobj tmpnum tmprvobj tmpallrv;' ] );
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.8 .15 .05], 'tag', 'tmp', ...
                  'style', 'pushbutton', 'string', 'Plot one', 'callback', ...
                   [ 'set(findobj(''parent'', gca, ''userdata'', ''dipole''), ''visible'', ''off'');' ...
                  'eval(get(findobj(''parent'', gcf, ''userdata'', ''editor''), ''callback''));' ]);
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.85 .15 .05], 'tag', 'tmp', ...
                  'style', 'pushbutton', 'string', 'Plot All', 'callback', ...
                   'set(findobj(''parent'', gca, ''userdata'', ''dipole''), ''visible'', ''on'');');
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.9 .15 .05], 'tag', 'tmp', ...
                  'style', 'text', 'string', [num2str(length(sources)) ' dipoles:'], 'fontweight', 'bold' );
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.95 .15 .05], 'tag', 'tmp', ...
                  'style', 'text', 'string', '   ' );
    set(gcf, 'userdata', findobj('parent', gca, 'tag', 'dipole1'));
    set(gca, 'userdata', textgui);
    set(gcf, 'color', BACKCOLOR);
    
    if strcmp(g.gui, 'off')
        set(findobj('parent', gcf, 'tag', 'tmp'), 'visible', 'off');
    end;
    if strcmp(g.mesh, 'off')
        set(findobj('parent', gca, 'tag', 'mesh'), 'visible', 'off');
    end;
        
    % HORIZONTAL LIMITS PLOTING
    % -------------------------
    %hold on;
    %for i=0:20:360
    %    [x y z]    = sph2cart(i/180*pi, 0, 100*scaling);
    %    h = plot3(x,  y,  z, 'g'); set(h, 'marker', '.', 'markersize', 30);
    %end;
    
    % VERTICAL LIMITS PLOTING
    % -----------------------
    %hold on;
    %for i=0:20:360
    %    [x y z]    = sph2cart(0, i/180*pi, 100*scaling);
    %    h = plot3(x,  y,  z, 'g'); set(h, 'marker', '.', 'markersize', 30);
    %end;

    try,
        x = cell2mat( {sources.X} );
        y = cell2mat( {sources.Y} );
        z = cell2mat( {sources.Z} );
        xe = cell2mat( {sources.XE} );
        ye = cell2mat( {sources.YE} );
        ze = cell2mat( {sources.ZE} );
    catch, end;
    rotate3d;
return;

    % 2D projection
    % -------------
    % internal constants
    scaling = 0.89;
    ori = 90;
    radiusoriinit = 35;

    figure; 
    warning off; a = imread(IMAGESLOC{1}); warning on;
    imagesc(a);
    colormap('gray');
    axis image;
    
    % plot projected coordinates
    tmpxlim = get(gca, 'xlim'); centerx = tmpxlim(2)/2+0.2;
    tmpylim = get(gca, 'ylim'); centery = tmpylim(2)/2+6;
    realradius = sin(source.sph_phi/180*pi)*source.sph_radius;
    absloc = centerx - cos((source.sph_theta+ori)/180*pi)*realradius*scaling;
    ordloc = centery - sin((source.sph_theta+ori)/180*pi)*realradius*scaling;
    absmax = centerx - cos((source.sph_theta+ori)/180*pi)*100*scaling;
    ordmax = centery - sin((source.sph_theta+ori)/180*pi)*100*scaling;
    circle(centerx, centery, 4, 'r', 'r');
    circle(absloc, ordloc, 4, 'k', 'k');
    circle(absmax, ordmax, 4, 'k', 'k');

    radiusori = sin(source.sph_phiori/180*pi)*radiusoriinit;
    absori = absloc - cos((source.sph_thetaori+ori)/180*pi)*radiusori*scaling;
    ordori = ordloc - sin((source.sph_thetaori+ori)/180*pi)*radiusori*scaling;
    h = line( [absloc absori]', [ordloc ordori]');
    set(h, 'color', 'k', 'linewidth', 4);

function h = myezplot3(strX, strY, strZ, range);
    figure; h = ezplot3(strX, strY, strZ, range);
    xdata = get(h, 'xdata');
    ydata = get(h, 'ydata');
    zdata = get(h, 'zdata');
    close;
    h = plot3(xdata, ydata, zdata);

function newsrc = convertbesaoldformat(src);
    newsrc = [];
    count = 1;
    if ~isfield(src, 'besaextori'), src.besaextori = []; end;
    for index = 1:length(src)
        countdip = 1;
        
        % convert format
        % --------------
        if isempty(src(index).besaextori), src(index).besaextori = 20; end; % 20 mm
        newsrc(index).possph(countdip,:) = [ src(index).besathloc src(index).besaphloc src(index).besaexent];
        newsrc(index).momsph(countdip,:) = [ src(index).besathori src(index).besaphori src(index).besaextori];
        
        % copy other fields
        % -----------------
        if isfield(src, 'rv')
            newsrc(index).rv = src(index).rv;
        end;
        if isfield(src, 'elecrv')
            newsrc(index).rvelec = src(index).elecrv;
        end;
        if isfield(src, 'component')
            newsrc(index).component = src(index).component;
            if index ~= 1 & src(index).component == src(index-1).component
                countdip = countdip + 1;
            else
                count = count + 1;
            end;
        else
            count = count + 1;
        end;
    end; 

function src = computexyzforbesa(src);
    scaling = 0.0105; % image scaling factor
    
    for index = 1:length( src )
        for index2 = size( src(index).possph, 1 )

            % compute coordinates
            % -------------------
            postmp = src(index).possph(index2,:);
            momtmp = src(index).momsph(index2,:);
            
            phi      = postmp(1)+90; %% %%%%%%%%%%%%%%% USE BESA COORDINATES %%%%%
            theta    = postmp(2);    %% %%%%%%%%%%%%%%% USE BESA COORDINATES %%%%%
            phiori   = momtmp(1)+90; %% %%%%%%%%%%%% USE BESA COORDINATES %%%%%
            thetaori = momtmp(2);    %% %%%%%%%%%%%% USE BESA COORDINATES %%%%%
            [x y z]    = sph2cart(theta/180*pi, phi/180*pi, postmp(3)*scaling);
            [xo yo zo] = sph2cart(thetaori/180*pi, phiori/180*pi, momtmp(3)*scaling);
            src(index).posxyz(index2,:) = [x y z];
            src(index).momxyz(index2,:) = [xo yo zo];
                    
        end;
    end;