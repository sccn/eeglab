% dipplot() - Visualize EEG equivalent-dipole locations and orientations 
%             in the BESA spherical head model or in a averaged MRI model.
% Usage:
%   >> dipplot( sources, 'key', 'val', ...);
%   >> [sources X Y Z XE YE ZE] = dipplot( sources, 'key', 'val', ...);
%
% Inputs:
%   sources   -  structure array of dipole information: can contain
%                either BESA or DIPFIT dipole information
%                besaexent: BESA eccentricity of the dipole
%                besathloc: BESA azimuth angle of the dipole
%                besaphloc: BESA horizontal angle of the dipole
%                besathori: BESA azimuth angle of the dipole orientation
%                besaphori: BESA horiz. angle of the dipole orientation
%                posxyz: DIPFIT dipole 3D carthesian position in mm
%                momxyz: DIPFIT dipole 3D carthesian orientation
%                optional fields for BESA and DIPFIT dipole info
%                     component: component number
%                     rv:        residual variance
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
%  'dipolelength' - Length of the dipole bar(s) {Default: 1}
%  'pointout' - ['on'|'off'] Point the dipoles outward. {Default: 'off'}
%  'sphere'   - [float] radius of sphere corresponding to the skin. Default is 1.
%  'normlen'  - ['on'|'off'] Normalize length of all dipoles. {Default: 'off'}
%  'std'      - [cell array] plot standard deviation of dipoles. i.e.
%               { [1:6] [7:12] } plot two elipsoids that best fit all the dipoles
%               from 1 to 6 and 7 to 12 with radius 1 standard deviation.
%               { { [1:6] 2 'linewidth' 2 } [7:12] } do the same but now the
%               first elipsoid is 2 standard-dev and the lines are thicker.
%
% Outputs:
%   sources   - EEG.source structure with updated 'X', 'Y' and 'Z' fields
%   X,Y,Z     - Locations of dipole heads (Cartesian coordinates). If there is
%               more than one dipole per components, the last dipole is returned.
%   XE,YE,ZE  - Locations of dipole ends (Cartesian coordinates). The same
%               remark as above applies.
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
% - All dipoles have a tag 'dipoleX' (X=their number) and can be made 
%     visible/invisible
% - The gcf object 'userdat' field stores the handle of the dipole that 
%     is currently being modified
% - Gca 'userdata' stores imqge names and position

%$Log: not supported by cvs2svn $
%Revision 1.32  2003/06/10 19:04:11  arno
%nothing
%
%Revision 1.31  2003/06/03 16:37:16  arno
%tag images
%
%Revision 1.30  2003/05/30 17:16:22  arno
%nothing
%
%Revision 1.29  2003/05/30 17:09:00  arno
%for index = 1:size(x, 2);
%    dipstruct(index).posxyz = [x(index) y(index) z(index)];
%    %dipstruct(index).posxyz = tmp(index).posxyz;
%    dipstruct(index).momxyz = [1 1 1];
%    dipstruct(index).component = index;
%    dipstruct(index).rv = 0.1;
%end;
%dipplot(dipstruct);
%making xyz output compatible
%
%Revision 1.28  2003/05/30 16:06:27  arno
%nothing
%
%Revision 1.27  2003/05/14 21:36:36  arno
%nothing
%
%Revision 1.26  2003/04/30 18:42:07  arno
%calibrating roughtly the slice selection
%
%Revision 1.25  2003/04/30 16:19:28  arno
%calibrating infants
%
%Revision 1.24  2003/04/30 02:05:24  arno
%changing axis properties for images
%
%Revision 1.23  2003/04/30 01:31:53  arno
%infant option
%
%Revision 1.22  2003/04/23 18:35:11  arno
%allow to plot elipses
%
%Revision 1.21  2003/04/22 21:18:44  arno
%standard dev
%
%Revision 1.20  2003/04/19 01:15:07  arno
%debugging 2 besa dipoles
%
%Revision 1.19  2003/04/19 00:55:53  arno
%correct normalized dipole length
%
%Revision 1.18  2003/04/19 00:46:41  arno
%correcting projection
%
%Revision 1.17  2003/04/19 00:37:43  arno
%changing dipole size for BESA
%
%Revision 1.16  2003/04/11 17:26:45  arno
%accurate plotting in fullMRI
%
%Revision 1.15  2003/04/10 17:37:18  arno
%multi layer MRI plot
%
%Revision 1.14  2003/03/14 17:06:42  arno
%adding error message if plotting non-BESA dipoles on the BESA head model
%
%Revision 1.13  2003/03/14 02:11:33  arno
%automatic scaling for dipfit
%
%Revision 1.12  2003/03/11 23:33:27  arno
%typo
%
%Revision 1.11  2003/03/11 23:27:09  arno
%adding normlen parameter
%
%Revision 1.10  2003/03/11 01:20:32  arno
%updating default besaextori, debuging summary
%
%Revision 1.9  2003/03/07 00:32:53  arno
%debugging textforgui
%
%Revision 1.8  2003/03/06 17:01:10  arno
%textgui -> textforgui
%
%Revision 1.7  2003/03/06 16:50:08  arno
%adding log message
%

function [outsources, XX, YY, ZZ, XO, YO, ZO] = dipplot( sourcesori, varargin )
    
    DEFAULTVIEW = [0 0 1];
        
    if nargin < 1
        help dipplot;
        return;
    end;
        
    % reading and testing arguments
    % -----------------------------
    sources = sourcesori;
    if ~isstruct(sources)
       updatedipplot(sources(1), sources(2)); 
       % sources(1) countain the figure handler, sources(2) the number of sources
       % (this is a callback function call to update the dipole plot)
       return
    end;
    
    %                             key        type       range             default
    g = finputcheck( varargin, { 'color'    ''         []                 [];
                                 'mesh'     'string'   { 'on' 'off' }     'off';
                                 'gui'      'string'   { 'on' 'off' }     'on';
                                 'summary'  'string'   { 'on' 'off' }     'off';
                                 'view'     'real'     []                 [];
                                 'rvrange'  'real'     [0 Inf]            [];
                                 'normlen'  'string'   { 'on' 'off' }     'off';
                                 'num'      'string'   { 'on' 'off' }     'off';
                                 'std'      'cell'     []                 {};
                                 'projimg'  'string'   { 'on' 'off' }     'off';
                                 'pointout'  'string'   { 'on' 'off' }     'off';
                                 'dipolesize' 'real'   [0 Inf]            30;
                                 'dipolelength' 'real' [0 Inf]            1;
                                 'sphere'   'real'     [0 Inf]              1;
                                 'image'    'string'   { 'besa' 'mri' 'mriinfant' 'fullmri'}   'besa' }, 'dipplot');
    if isstr(g), error(g); end;
    
    % axis image and limits
    % ---------------------
    if strcmpi(g.image, 'besa')
        IMAGESLOC   = { { 'besatop.pcx' } { 'besarear.pcx' } { 'besaside.pcx' } };
        IMAGESAXIS  = { -1  1 -1 };
        IMAGESOFFSET  = { [0.0  0.08 NaN ]  [0.01    NaN -0.01]  [ NaN -0.01   -0.025 ] };
        IMAGESMULT    = { 1.01 1.06 0.96 };
        COLORMESH = [.5 .5 .5];
        BACKCOLOR = 'w';
        AXISLIM   = [-1.2 1.2 -1.2 1.2 -1.2 1.2];
    elseif strcmpi(g.image, 'mri')
        IMAGESLOC   = { { 'mritop.pcx' } { 'mrirear.pcx' } { 'mriside.pcx' } };
        IMAGESAXIS  = { -1  1 -1 };
        IMAGESOFFSET = { [-0.01 0.005  NaN]   [-0.02 NaN 0.11]  [NaN 0.04 0.30] } ;% [ -0.12 0] };
        IMAGESMULT   = { 0.8 0.8 1 };
        COLORMESH = 'w';
        BACKCOLOR = 'k';
        %AXISLIM   = [-1.2 1.2 -1.2 1.2 -1.2 1.2];
        AXISLIM   = [-1.4 1.4 -1.1 1.1 -1.2 1.2];
    elseif strcmpi(g.image, 'mriinfant')
        IMAGESLOC   = { { 'transv_infant.pcx' } { 'cor_infant.pcx' } { 'sag_infant.pcx' } };
        IMAGESAXIS  = { -1  1 -1 };
        IMAGESOFFSET = { [-0.01 0.08  NaN]   [-0.02 NaN 0.05]  [NaN 0.05 0.1] } ;
        IMAGESMULT   = { 0.9 0.91 1.05 } ;
        COLORMESH = 'w';
        BACKCOLOR = 'k';
        AXISLIM   = [-1.25 1.25 -1.1 1.1 -1.2 1.2];
    else 
       [IMAGESLOC IMAGESAXIS] = getmriimgs;
       addpath('/data/common/matlab/MRIimages');
       IMAGESOFFSET = { [-0.01 0.005  NaN]   [-0.02 NaN 0.11]  [NaN 0 0.17] } ;% [ -0.12 0] };
       IMAGESMULT   = { 1.2 1.2 1.2 } ;%[ 1.4 1.15 ] };
       %IMAGESOFFSET = { [0 0  NaN]   [0 NaN 0]  [NaN 0 0] } ;% [ -0.12 0] };
       %IMAGESMULT   = { [1 1  NaN]   [1 NaN 1]  [NaN 1 1]} ;%[ 1.4 1.15 ] };
       COLORMESH = 'w';
       BACKCOLOR = 'k';
       %AXISLIM   = [-1.2 1.2 -1.2 1.2 -1.2 1.2];
       AXISLIM   = [-1.4 1.4 -1.1 1.1 -1.2 1.2];
    end;

    % conversion
    % ----------
    if strcmpi(g.normlen, 'on')
        try, sources = rmfield(sources, 'besaextori'); catch, end;
    end;
    if ~isfield(sources, 'besathloc') & strcmpi(g.image, 'besa') & ~is_sccn
        error(['For copyright reasons, it is not possible to use the BESA ' ...
               'head model to plot non-BESA dipoles']);
    end;
    
    if isfield(sources, 'besathloc')
        sources = convertbesaoldformat(sources);
    end;
    if ~isfield(sources, 'posxyz')
        sources = computexyzforbesa(sources);
    end;        

    % normalize position to unit sphere
    % ---------------------------------
    maxi = 0;
    for index = 1:length(sources)
        maxi = max(maxi,max(abs(sources(index).posxyz(:))));
    end;
    if maxi > 1.01
        disp('Non-normalized dipole positions, normalizing by standard head radius 84.747 mm'); 
        g.sphere = 84.747;
        fact = 0.1;
    else 
        fact = 1;
    end;
    
    % transform coordinates
    % ---------------------
    outsources = sources;
    for index = 1:length(sources)
        sources(index).posxyz = sources(index).posxyz/g.sphere;
        tmp = sources(index).posxyz(:,1);
        sources(index).posxyz(:,1) = sources(index).posxyz(:,2);
        sources(index).posxyz(:,2) = -tmp;
        sources(index).momxyz = sources(index).momxyz/g.sphere*0.05*fact;
        tmp = sources(index).momxyz(:,1);
        sources(index).momxyz(:,1) = sources(index).momxyz(:,2);
        sources(index).momxyz(:,2) = -tmp;
        if isfield(sources, 'stdX')
            tmp = sources(index).stdX;
            sources(index).stdX = sources(index).stdY;
            sources(index).stdY = -tmp;
        end;
        if strcmpi(g.normlen, 'on')
            warning off;
            sources(index).momxyz(1,:) = 0.2*sources(index).momxyz(1,:)/ norm(abs(sources(index).momxyz(1,:)));
            if size(sources(index).momxyz,1) > 1 & sources(index).momxyz(1) ~= 0
                sources(index).momxyz(2,:) = 0.2*sources(index).momxyz(2,:)/ norm(abs(sources(index).momxyz(2,:)));
            end;
            warning on;
        end;
    end;
    
    % remove sources with out of bound Residual variance
    % --------------------------------------------------
    if isfield(sources, 'rv') & ~isempty(g.rvrange)
        for index = length(sources):-1:1
            if sources(index).rv < g.rvrange(1) | sources(index).rv > g.rvrange(2)
                sources(index) = [];
            end;
        end;
    end;
    
    % build summarized figure
    % -----------------------
    if strcmp(g.summary, 'on')
        figure;
        options = { 'gui', 'off', 'dipolesize', g.dipolesize,'dipolelength', g.dipolelength, ...
                    'color', g.color, 'mesh', g.mesh, 'num', g.num, 'image', g.image };
        axes('position', [0 0 0.5 0.5]);  dipplot(sourcesori, 'view', [1 0 0] , options{:}); axis off;
        axes('position', [0 0.5 0.5 .5]); dipplot(sourcesori, 'view', [0 0 1] , options{:}); axis off;
        axes('position', [.5 .5 0.5 .5]); dipplot(sourcesori, 'view', [0 -1 0], options{:}); axis off;
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
                        textforgui(colorcount) = { sprintf('Component %d (R.V. %3.2f)', sources(index).component, 100*sources(index).rv) };
                        colorcount = colorcount+1;
                    end;
                else 
                    textforgui(colorcount) = { sprintf('Component %d (R.V. %3.2f)', sources(index).component, 100*sources(index).rv) };
                    colorcount = colorcount+1;
                end;
            end;
            colorcount = colorcount-1;
            allstr = strvcat(textforgui{:});
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
	 plotimgs(IMAGESLOC, IMAGESOFFSET, IMAGESMULT, IMAGESAXIS, AXISLIM, [1 1 1]);
    set(gca, 'color', BACKCOLOR);
    %warning off; a = imread('besaside.pcx'); warning on;
    % BECAUSE OF A BUG IN THE WARP FUNCTION, THIS DOES NOT WORK (11/02)
    %hold on; warp([], wy, wz, a);
    if ~isempty(g.view)
         view(g.view);
    else view(DEFAULTVIEW);
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

            % copy for output
            % ---------------
            XX(index) = -y;
            YY(index) = x;
            ZZ(index) = z;
            XO(index) = -yo;
            YO(index) = xo;
            ZO(index) = zo;
            
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
            dipstruct.pos3d = [x y z];
            dipstruct.rv    = sprintf('C %d (%3.2f)', sources(index).component, sources(index).rv*100);
            set(h, 'userdata', dipstruct, 'tag', tag, 'color','k', 'linewidth', g.dipolesize/7.5);
            if strcmp(BACKCOLOR, 'k'), set(h, 'color', g.color{index}); end;
            
            % draw point
            hold on;
            h = plot3(x,  y,  z); 
            set(h, 'userdata', dipstruct, 'tag', tag, ...
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
                if ~strcmpi(g.image, 'besa')
                    h = line( [x xo]', [y yo]', [-1 -1]');
                    set(h, 'userdata', 'proj', 'tag', tag, 'color','k', 'linewidth', g.dipolesize/7.5);
                end;
                if strcmp(BACKCOLOR, 'k'), set(h, 'color', tmpcolor); end;
                h = plot3(x,  y,  -1); 
                set(h, 'userdata', 'proj', 'tag', tag, ...
                       'marker', '.', 'markersize', g.dipolesize, 'color', tmpcolor);
                
                % project onto x axis
                tag = [ 'dipole' num2str(index) ];
                if ~strcmpi(g.image, 'besa')
                    h = line( [x xo]', [1 1]', [z zo]');
                    set(h, 'userdata', 'proj', 'tag', tag, 'color','k', 'linewidth', g.dipolesize/7.5);
                end;
                if strcmp(BACKCOLOR, 'k'), set(h, 'color', tmpcolor); end;
                h = plot3(x,  1,  z); 
                set(h, 'userdata', 'proj', 'tag', tag, ...
                       'marker', '.', 'markersize', g.dipolesize, 'color', tmpcolor);
                
                % project onto y axis
                tag = [ 'dipole' num2str(index) ];
                if ~strcmpi(g.image, 'besa')
                    h = line( [-1 -1]', [y yo]', [z zo]');
                    set(h, 'userdata', 'proj', 'tag', tag, 'color','k', 'linewidth', g.dipolesize/7.5);
                end;
                if strcmp(BACKCOLOR, 'k'), set(h, 'color', tmpcolor); end;
                h = plot3(-1,  y,  z); 
                set(h, 'userdata', 'proj', 'tag', tag, ...
                       'marker', '.', 'markersize', g.dipolesize, 'color', tmpcolor);
            end;
                        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% draw text  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if isfield(sources, 'component')
                if strcmp(g.num, 'on')
                    h = text(x,  y,  z, [ '  ' int2str(sources(index).component)]);
                    set(h, 'userdata', dipstruct, 'tag', tag, 'fontsize', g.dipolesize/2 );
                end;
            end;
        end;
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% draw elipse for group of dipoles  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~isempty(g.std)
        for index = 1:length(g.std)
            if ~iscell(g.std{index})
                plotellipse(sources, g.std{index}, 1);
            else
                sc = plotellipse(sources, g.std{index}{1}, g.std{index}{2});
                if length( g.std{index} ) > 2
                    set(sc, g.std{index}{3:end});
                end;
            end;
        end;
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% buttons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nbsrc = int2str(length(sources));
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
                  'dipplot([ gcbf ' nbsrc ' ]);' ...
                  'clear editobj;' ]);
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.65 .15 .05], 'tag', 'tmp', ...
                  'style', 'pushbutton', 'string', 'Keep|Next', 'callback', ...
                [ 'editobj = findobj(''parent'', gcf, ''userdata'', ''editor'');' ...
                  'set(editobj, ''string'', num2str(str2num(get(editobj, ''string''))+1));' ...
                  'tmpobj = get(gcf, ''userdata'');' ...
                  'dipplot([ gcbf ' nbsrc ' ]);' ...
                  'set(tmpobj, ''visible'', ''on'');' ...
                  'clear editobj tmpobj;' ]);
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.7 .15 .05], 'tag', 'tmp', 'userdata', 'rv', ...
                  'style', 'text', 'string', '');
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.75 .15 .05], 'tag', 'tmp', 'userdata', 'editor', ...
                  'style', 'edit', 'string', '1', 'callback', ...
                  [ 'dipplot([ gcbf ' nbsrc ' ]);' ] );
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.8 .15 .05], 'tag', 'tmp', ...
                   'style', 'pushbutton', 'string', 'Plot one', 'callback', ...
             	    [ 'for tmpi = 1:' nbsrc ',' ...
                   '   set(findobj(''parent'', gca, ''tag'', [ ''dipole'' int2str(tmpi) ]), ''visible'', ''off'');' ...
                   'end; clear tmpi;' ...
                   'dipplot([ gcbf ' nbsrc ' ]);' ]);
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.85 .15 .05], 'tag', 'tmp', ...
                   'style', 'pushbutton', 'string', 'Plot All', 'callback', ...
                   [ 'for tmpi = 1:' nbsrc ',' ...
                   '   set(findobj(''parent'', gca, ''tag'', [ ''dipole'' int2str(tmpi) ]), ''visible'', ''on'');' ...
                   'end; clear tmpi;' ]);
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.9 .15 .05], 'tag', 'tmp', ...
                  'style', 'text', 'string', [num2str(length(sources)) ' dipoles:'], 'fontweight', 'bold' );
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.95 .15 .05], 'tag', 'tmp', ...
                  'style', 'text', 'string', '   ' );
    set(gcf, 'userdata', findobj('parent', gca, 'tag', 'dipole1'));
    set(gca, 'userdata', { IMAGESLOC IMAGESOFFSET, IMAGESMULT IMAGESAXIS AXISLIM } );
    set(gcf, 'color', BACKCOLOR);
    
    if strcmp(g.gui, 'off')
        set(findobj('parent', gcf, 'tag', 'tmp'), 'visible', 'off');
    end;
    if strcmp(g.mesh, 'off')
        set(findobj('parent', gca, 'tag', 'mesh'), 'visible', 'off');
    end;
        
    rotate3d;
return;

    % 2D projection DEPRECATED
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

function sc = plotellipse(sources, ind, nstd);

    for i = 1:length(ind)
        tmpval(1,i) = -sources(ind(i)).posxyz(1);    
        tmpval(2,i) = -sources(ind(i)).posxyz(2);    
        tmpval(3,i) = sources(ind(i)).posxyz(3);    
    end;
    
    % mean and covariance
    C = cov(tmpval');
    M = mean(tmpval,2);
    [U,L] = eig(C);
    
    % For N standard deviations spread of data, the radii of the eliipsoid will
    % be given by N*SQRT(eigenvalues).
    radii = nstd*sqrt(diag(L));
    
    % generate data for "unrotated" ellipsoid
    [xc,yc,zc] = ellipsoid(0,0,0,radii(1),radii(2),radii(3), 10);
    
    % rotate data with orientation matrix U and center M
    a = kron(U(:,1),xc); b = kron(U(:,2),yc); c = kron(U(:,3),zc);
    data = a+b+c;  n = size(data,2);
    x = data(1:n,:)+M(1); y = data(n+1:2*n,:)+M(2); z = data(2*n+1:end,:)+M(3);
    
    % now plot the rotated ellipse
    c = ones(size(z));
    sc = mesh(x,y,z);
    alpha(0.5)
    
function newsrc = convertbesaoldformat(src);
    newsrc = [];
    count = 1;
    countdip = 1;
    if ~isfield(src, 'besaextori'), src(1).besaextori = []; end;
    for index = 1:length(src)
        
        % convert format
        % --------------
        if isempty(src(index).besaextori), src(index).besaextori = 300; end; % 20 mm
        newsrc(count).possph(countdip,:) = [ src(index).besathloc src(index).besaphloc src(index).besaexent];
        newsrc(count).momsph(countdip,:) = [ src(index).besathori src(index).besaphori src(index).besaextori/300];
        
        % copy other fields
        % -----------------
        if isfield(src, 'stdX')
            newsrc(count).stdX = -src(index).stdY;
            newsrc(count).stdY = src(index).stdX;
            newsrc(count).stdZ = src(index).stdZ;
        end;
        if isfield(src, 'rv')
            newsrc(count).rv = src(index).rv;
        end;
        if isfield(src, 'elecrv')
            newsrc(count).rvelec = src(index).elecrv;
        end;
        if isfield(src, 'component')
            newsrc(count).component = src(index).component;
            if index ~= length(src) & src(index).component == src(index+1).component
                countdip = countdip + 1;
            else
                count = count + 1; countdip = 1;
            end;
        else
            count = count + 1; countdip = 1;
        end;
    end; 

function src = computexyzforbesa(src);
    
    for index = 1:length( src )
        for index2 = 1:size( src(index).possph, 1 )

            % compute coordinates
            % -------------------
            postmp = src(index).possph(index2,:);
            momtmp = src(index).momsph(index2,:);
            
            phi      = postmp(1)+90; %% %%%%%%%%%%%%%%% USE BESA COORDINATES %%%%%
            theta    = postmp(2);    %% %%%%%%%%%%%%%%% USE BESA COORDINATES %%%%%
            phiori   = momtmp(1)+90; %% %%%%%%%%%%%% USE BESA COORDINATES %%%%%
            thetaori = momtmp(2);    %% %%%%%%%%%%%% USE BESA COORDINATES %%%%%
            % exentricities are in % of the radius of the head sphere
            [x y z]    = sph2cart(theta/180*pi, phi/180*pi, postmp(3)/100); 
            [xo yo zo] = sph2cart(thetaori/180*pi, phiori/180*pi, momtmp(3)*5); % exentricity scaled for compatibility with DIPFIT
            src(index).posxyz(index2,:) = [-y x z];
            src(index).momxyz(index2,:) = [-yo xo zo];
                    
        end;
     end;
     
% update dipplot (callback call)
% ------------------------------
function updatedipplot(fig, nbsources)
   
   % find current dipole index and test for authorized range
   editobj = findobj('parent', fig, 'userdata', 'editor');
   tmpnum = str2num(get(editobj, 'string'));
   if tmpnum < 1, tmpnum = 1; end;
   if tmpnum > nbsources, tmpnum = nbsources; end;
   set(editobj, 'string', num2str(tmpnum));
   
   % hide current dipole
   set(get(gcf, 'userdata'), 'visible', 'off');
   
   % find next dipole and show it
   newdip = findobj('parent', gca, 'tag', [ 'dipole' get(editobj, 'string')]);
   set(newdip, 'visible', 'on');
   set(gcf, 'userdata', newdip);
   
   % set the new RV
   tmprvobj = findobj('parent', gcf, 'userdata', 'rv');
   userdat = get(newdip(1), 'userdata');
   set( tmprvobj, 'string', userdat.rv);
   
   % adapt the MRI to the dipole depth
   imginfos  = get(gca, 'userdata');
   imglocs   = imginfos{1};
   imgoffset = imginfos{2};
   imgmult   = imginfos{3};
   imgaxis   = imginfos{4};
   axislim   = imginfos{5};
   if length(imginfos{1}) > 1 % several images
      delete(findobj('parent', gca, 'tag', 'img'));
      
      tmp = userdat.pos3d(1)
      indx = minpos(imgaxis{1} - userdat.pos3d(3) + 8/84.747);
      indy = minpos(imgaxis{2} - userdat.pos3d(2) - 4/84.747);
      indz = minpos(imgaxis{3} - userdat.pos3d(1) + 4/84.747);
      tmp = userdat.pos3d*84.747
      [ imgaxis{3}(indx) imgaxis{2}(indy) imgaxis{1}(indz)]*84.747
      plotimgs( imglocs, imgoffset, imgmult, imgaxis, axislim, [indx indy indz]);
   end;
   	
% plot images
% -----------
function plotimgs(IMAGESLOC, IMAGESOFFSET, IMAGESMULT, IMAGESAXIS, AXISLIM, index);
   
    %try,
        fprintf('Reading img: %s\n', IMAGESLOC{1}{index(1)} );
        warning off;  a = double(imread(IMAGESLOC{1}{index(1)}))/255; warning on;
        if ndims(a) == 2, a(:,:,2) = a; a(:,:,3) = a(:,:,1); end;
        aspect_ratio   = size(a,1)/size(a,2);
        wx = ([-1 1; -1 1])*IMAGESMULT{1}*AXISLIM(2)+IMAGESOFFSET{1}(1);
        wy = ([-1 -1; 1 1])*IMAGESMULT{1}*AXISLIM(4)*aspect_ratio+IMAGESOFFSET{1}(2);
        wz = [ 1 1; 1 1]*IMAGESAXIS{1}(index(1))*1.07;
        % DISPLAY IF REAL RATIO OF IMAGE CORRESPOND TO RATIO SHOWN ON SCREEN
        %fprintf('Image ratio %3.2f\tCoord ratio:%3.2f\n', size(a,2)/size(a,1),  ...
        %       IMAGESMULT{1}(1)/ IMAGESMULT{1}(2) / (AXISLIM(2)-AXISLIM(1)) * (AXISLIM(4)-AXISLIM(3)) );
        surface(wx, wy, wz, a(end:-1:1,:,:), 'FaceColor','texturemap', ...
           'EdgeColor','none', 'CDataMapping','direct','tag','img', 'facelighting', 'none');
        hold on; %%fill3([-2 -2 2 2], [-2 2 2 -2], wz(:)-1, BACKCOLOR);
    %catch, error(lasterr); end;
    try,
        fprintf('Reading img: %s\n', IMAGESLOC{2}{index(2)} );
        warning off; a = double(imread(IMAGESLOC{2}{index(2)}))/255;  warning on;
        if ndims(a) == 2, a(:,:,2) = a; a(:,:,3) = a(:,:,1); end;
        aspect_ratio   = size(a,1)/size(a,2);
        wx = ([-1 1; -1 1])*IMAGESMULT{2}*AXISLIM(2)+IMAGESOFFSET{2}(1);
        wz = ([1 1; -1 -1])*IMAGESMULT{2}*AXISLIM(6)*aspect_ratio+IMAGESOFFSET{2}(3);
        wy = [1 1; 1 1]*IMAGESAXIS{2}(index(2))*1.07;
        hold on; surface(wx, wy, wz, a, 'FaceColor','texturemap', ...
           'EdgeColor','none', 'CDataMapping','direct','tag','img', 'facelighting', 'none');
        %fprintf('Image ratio %3.2f\tCoord ratio:%3.2f\n', size(a,2)/size(a,1),  ...
        %        IMAGESMULT{2}(1)/ IMAGESMULT{2}(3) / (AXISLIM(2)-AXISLIM(1)) * (AXISLIM(6)-AXISLIM(5)));
        %%fill3([-2 -2 2 2], wy(:)-1, [-2 2 2 -2], BACKCOLOR);
    catch, error(lasterr); end;
    try,
        fprintf('Reading img: %s\n', IMAGESLOC{3}{index(3)} );
        warning off; a = double(imread(IMAGESLOC{3}{index(3)}))/255;  warning on;
        if ndims(a) == 2, a(:,:,2) = a; a(:,:,3) = a(:,:,1); end;
        aspect_ratio   = size(a,1)/size(a,2);
        wx = [ 1  1;  1 1]*IMAGESAXIS{3}(index(3))*1.07;
        wy = ([-1 1; -1 1])*IMAGESMULT{3}*AXISLIM(4)+IMAGESOFFSET{3}(2);
        wz = ([1 1; -1 -1])*IMAGESMULT{3}*AXISLIM(6)*aspect_ratio+IMAGESOFFSET{3}(3);
        hold on; surface(wx, wy, wz, a, 'FaceColor','texturemap', ...
           'EdgeColor','none', 'CDataMapping','direct','tag','img', 'facelighting', 'none');
        %fprintf('Image ratio %3.2f\tCoord ratio:%3.2f\n', size(a,2)/size(a,1), ...
        %        IMAGESMULT{3}(2)/ IMAGESMULT{3}(3) / (AXISLIM(4)-AXISLIM(3)) * (AXISLIM(6)-AXISLIM(5)));
        %%fill3([-2 -2 2 2], wy(:)-1, [-2 2 2 -2], BACKCOLOR);
    catch, error(lasterr); end;
    rotate3d
    
% return MRI images
function [LOCS, AXIS] = getmriimgs;
coronal = { ...		
'C13.JPG'	-100	;
'C15.JPG'	-96	;
'C17.JPG'	-92	;
'C19.JPG'	-88	;
'C21.JPG'	-84	;
'C23.JPG'	-80	;
'C25.JPG'	-76	;
'C27.JPG'	-72	;
'C29.JPG'	-68	;
'C31.JPG'	-64	;
'C33.JPG'	-60	;
'C35.JPG'	-56	;
'C37.JPG'	-52	;
'C39.JPG'	-48	;
'C41.JPG'	-44	;
'C43.JPG'	-40	;
'C45.JPG'	-36	;
'C47.JPG'	-32	;
'C49.JPG'	-28	;
'C51.JPG'	-24	;
'C53.JPG'	-20	;
'C55.JPG'	-16	;
'C57.JPG'	-12	;
'C59.JPG'	-8	;
'C61.JPG'	-4	;
'C63.JPG'	0	;
'C65.JPG'	4	;
'C67.JPG'	8	;
'C69.JPG'	12	;
'C71.JPG'	16	;
'C73.JPG'	20	;
'C75.JPG'	24	;
'C77.JPG'	28	;
'C79.JPG'	32	;
'C81.JPG'	36	;
'C83.JPG'	40	;
'C85.JPG'	44	;
'C87.JPG'	48	;
'C89.JPG'	52	;
'C91.JPG'	56	;
'C93.JPG'	60	;
'C95.JPG'	64	;
'C97.JPG'	68	;
'C99.JPG'	72	;
'C101.JPG'	76	;
'C103.JPG'	80	;
'C105.JPG'	84	;
'C107.JPG'	88	;
'C109.JPG'	92	;
'C111.JPG'	96	;
'C113.JPG'	100	;
'C115.JPG'	104	;
'C117.JPG'	108	};
sagital = { ...		
'S20.JPG'	-88	;
'S22.JPG'	-84	;
'S24.JPG'	-80	;
'S26.JPG'	-76	;
'S28.JPG'	-72	;
'S30.JPG'	-68	;
'S32.JPG'	-64	;
'S34.JPG'	-60	;
'S36.JPG'	-56	;
'S38.JPG'	-52	;
'S40.JPG'	-48	;
'S42.JPG'	-44	;
'S44.JPG'	-40	;
'S46.JPG'	-36	;
'S48.JPG'	-32	;
'S50.JPG'	-28	;
'S52.JPG'	-24	;
'S54.JPG'	-20	;
'S56.JPG'	-16	;
'S58.JPG'	-12	;
'S60.JPG'	-8	;
'S62.JPG'	-4	;
'S64.JPG'	0	;
'S66.JPG'	4	;
'S68.JPG'	8	;
'S70.JPG'	12	;
'S72.JPG'	16	;
'S74.JPG'	20	;
'S76.JPG'	24	;
'S78.JPG'	28	;
'S80.JPG'	32	;
'S82.JPG'	36	;
'S84.JPG'	40	;
'S86.JPG'	44	;
'S88.JPG'	48	;
'S90.JPG'	52	;
'S92.JPG'	56	;
'S94.JPG'	60	;
'S96.JPG'	64	;
'S98.JPG'	68	;
'S100.JPG'	84	;
'S102.JPG'	88	;
'S104.JPG'	92	;
'S106.JPG'	96	;
'S108.JPG'	100	};
transverse = { ...		
'T17.JPG'	-72	;
'T19.JPG'	-68	;
'T21.JPG'	-64	;
'T23.JPG'	-60	;
'T25.JPG'	-56	;
'T27.JPG'	-52	;
'T29.JPG'	-48	;
'T31.JPG'	-44	;
'T33.JPG'	-40	;
'T35.JPG'	-36	;
'T37.JPG'	-32	;
'T39.JPG'	-28	;
'T41.JPG'	-24	;
'T43.JPG'	-20	;
'T45.JPG'	-16	;
'T47.JPG'	-12	;
'T49.JPG'	-8	;
'T51.JPG'	-4	;
'T53.JPG'	0	;
'T55.JPG'	4	;
'T57.JPG'	8	;
'T59.JPG'	12	;
'T61.JPG'	16	;
'T63.JPG'	20	;
'T65.JPG'	24	;
'T67.JPG'	28	;
'T69.JPG'	32	;
'T71.JPG'	36	;
'T73.JPG'	40	;
'T75.JPG'	44	;
'T77.JPG'	48	;
'T79.JPG'	52	;
'T81.JPG'	56	;
'T83.JPG'	60	;
'T85.JPG'	64	;
'T87.JPG'	68	;
'T89.JPG'	72	;
'T91.JPG'	76	;
'T93.JPG'	80	;
'T95.JPG'	84	;
'T97.JPG'	88	;
'T99.JPG'	92	;
'T101.JPG'	72	;
'T103.JPG'	76	;
'T105.JPG'	80	};
LOCS = { transverse(:,1)' coronal(:,1)' sagital(:,1)' };
AXIS = { cell2mat(transverse(:,2)')/84.747 cell2mat(coronal(:,2)')/84.747 cell2mat(sagital(:,2)')/84.747 };

function index = minpos(vals);
	vals(find(vals < 0)) = inf;
	[tmp index] = min(vals);
