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
%  'view'     - 3-D viewing angle in cartesian coords.,
%               [0 0 1] gives a sagittal view, [0 -1 0] a view from the rear;
%               [1 0 0] gives a view from the side of the head.
%  'mesh'     - ['on'|'off'] Display spherical mesh. {Default is 'on'}
%  'axistight' - ['on'|'off'] For MRI only, display the closest MRI
%               slide. {Default is 'off'}
%  'gui'      - ['on'|'off'] Display controls. {Default is 'on'} If gui 'off', 
%               a new figure is not created. Useful for incomporating a dipplot 
%               into a complex figure.
%  'num'      - ['on'|'off'] Display component number. Take into account
%               dipole size. {Default: 'off'}
%  'summary'  - Build a summary plot with three views (top, back, side)
%  'image'    - ['besa'|'mri'] Background image. {Default: 'besa'} 'mri' uses
%               mean-MRI brain images from the Montreal Neurological Institute.
%  'cornermri' - ['on'|'off'] force MRI images to the corner of the MRI volume
%               (usefull when background is not black). Default: 'off'.
%  'drawedges' - ['on'|'off'] draw edges of the 3-D MRI (black in axistight,
%                white otherwise.) Default is 'off'.
%  'rvrange'  - [min max] Only plot dipoles with residual variace within the
%               given range. Default: plot all dipoles.
%  'projimg'  - ['on'|'off'] Project dipole(s) onto the 2-D images, for use
%               in making 3-D plots {Default 'off'}
%  'projlines' - ['on'|'off'] Plot lines connecting dipole with 2-D projection.
%                Color is dashed black for BESA head and dashed black for the
%                MNI brain {Default 'off'}
%  'projcol'   - [color] color for the projected line {Default is same as dipole}
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
%  sources(1).rv = 0.1;
%  % position and orientation of the second dipole
%  sources(2).besaexent= 69.036;
%  sources(2).besathloc= 46;
%  sources(2).besaphloc= 39;
%  sources(2).besathori= 150;
%  sources(2).besaphori= -3;
%  sources(2).rv = 0.05;
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
%Revision 1.67  2004/02/23 19:16:16  arno
%fixing 2 dipole problem
%
%Revision 1.66  2004/02/23 19:09:05  arno
%*** empty log message ***
%
%Revision 1.65  2004/01/28 15:26:14  arno
%detecting empty dipoles
%
%Revision 1.64  2003/11/05 20:30:51  arno
%nothing
%
%Revision 1.63  2003/11/05 18:50:28  arno
%nothing
%
%Revision 1.62  2003/10/31 23:29:29  arno
%dicreases dipole size for summary mode
%
%Revision 1.61  2003/10/31 19:50:38  arno
%reordering 2 dipoles, delete old edges
%
%Revision 1.60  2003/10/31 18:59:41  arno
%divide rvrange by 100
%
%Revision 1.59  2003/10/30 03:25:17  arno
%adding projection lines
%
%Revision 1.58  2003/10/29 23:44:38  arno
%don't know
%
%Revision 1.57  2003/10/27 18:50:06  arno
%adding edges
%
%Revision 1.56  2003/10/09 01:02:12  arno
%nothing
%
%Revision 1.55  2003/09/11 00:53:04  arno
%same
%
%Revision 1.54  2003/09/11 00:52:38  arno
%documenting new arg cornermri\
%
%Revision 1.53  2003/09/11 00:51:19  arno
%adding new param cornermri
%
%Revision 1.52  2003/08/13 00:46:07  arno
%*** empty log message ***
%
%Revision 1.51  2003/08/04 21:29:06  arno
%scale besa for summary plot
%
%Revision 1.50  2003/08/04 21:15:54  arno
%num color on MRI
%
%Revision 1.49  2003/08/04 19:08:09  arno
%[Afixing normlen in summarize mode
%
%Revision 1.48  2003/07/31 23:37:07  arno
%adding property to structure attatched to dipole
%
%Revision 1.47  2003/07/31 23:30:05  arno
%*** empty log message ***
%
%Revision 1.46  2003/07/24 17:08:03  arno
%making component number for dipfit
%
%Revision 1.45  2003/07/22 01:10:06  arno
%*** empty log message ***
%
%Revision 1.44  2003/07/21 22:04:09  arno
%debug for Matlab 5.3
%
%Revision 1.43  2003/07/21 21:55:14  arno
%fixing MNI brain for distribution
%
%Revision 1.42  2003/07/16 18:39:14  arno
%nothing
%
%Revision 1.41  2003/07/03 00:02:43  arno
%undo last change
%
%Revision 1.40  2003/07/02 23:48:56  arno
%changing outputs
%
%Revision 1.39  2003/07/02 23:38:58  arno
%debuging projections and summary
%
%Revision 1.38  2003/07/01 23:52:50  arno
%test for sphere 1 before renormalizing
%
%Revision 1.37  2003/07/01 23:49:21  arno
%try to implement contextual menu: does not work in 3-D
%
%Revision 1.36  2003/07/01 22:10:14  arno
%debuging for 2 dipoles/component
%
%Revision 1.35  2003/07/01 19:04:13  arno
%fixing view problem
%
%Revision 1.34  2003/07/01 00:21:18  arno
%implementing 3-D MNI volume
%
%Revision 1.33  2003/06/12 23:49:56  arno
%dipplot normalization
%
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
        updatedipplot(sources(1)); 
        % sources countain the figure handler
        return
   end;
    
    %                             key        type       range             default
    g = finputcheck( varargin, { 'color'     ''         []                 [];
                                 'axistight' 'string'  { 'on' 'off' }     'off';
                                 'drawedges' 'string'  { 'on' 'off' }     'off';
                                 'mesh'      'string'   { 'on' 'off' }     'off';
                                 'gui'       'string'   { 'on' 'off' }     'on';
                                 'summary'   'string'   { 'on' 'off' }     'off';
                                 'view'      'real'     []                 [0 0 1];
                                 'rvrange'   'real'     [0 Inf]            [];
                                 'normlen'   'string'   { 'on' 'off' }     'off';
                                 'num'       'string'   { 'on' 'off' }     'off';
                                 'cornermri' 'string'   { 'on' 'off' }     'off';
                                 'std'       'cell'     []                 {};
                                 'projimg'   'string'   { 'on' 'off' }     'off';
                                 'projcol'   { 'string' 'real' }   []      [];
                                 'projlines' 'string'   { 'on' 'off' }     'off';
                                 'pointout'  'string'   { 'on' 'off' }     'off';
                                 'dipolesize' 'real'   [0 Inf]            30;
                                 'dipolelength' 'real' [0 Inf]            1;
                                 'sphere'    'real'     [0 Inf]              1;
                                 'image'     'string'   { 'besa' 'mri' 'mriinfant' 'fullmri'}   'besa' }, 'dipplot');
    if isstr(g), error(g); end;
    g.zoom = 1500;
    
    % axis image and limits
    % ---------------------
    dat.mode       = g.image;
    dat.maxcoord   = [ 90 100 100 ]; % location of images in 3-D, Z then Y then X
    dat.axistight  = strcmpi(g.axistight, 'on');
    dat.drawedges  = g.drawedges;
    radius = 85;
    if strcmpi(g.image, 'besa')
        scaling = 1.05;
        
        % read besa images
        % ----------------
        warning off; imgt = double(imread('besatop.pcx' ))/255; warning on;
        warning off; imgc = double(imread('besarear.pcx'))/255; warning on;
        warning off; imgs = double(imread('besaside.pcx'))/255; warning on;
        dat.imgs        = { imgt imgc imgs };
        
        allcoords1 = ([-1.12 1.12]*size(imgt,2)/200+0.01)*radius; 
        allcoords2 = ([-1.12 1.12]*size(imgt,1)/200+0.08)*radius;
        allcoords3 = [-1.12 1.12]*size(imgs,1)/200*radius; 
        dat.imgcoords = { allcoords3        allcoords2        allcoords1 };

        valinion  = [ 1  0  0 ]*radius;
        valnasion = [-1  0  0 ]*radius;
        vallear   = [ 0 -1  0 ]*radius;
        valrear   = [ 0  1  0 ]*radius;
        valvertex = [ 0  0  1 ]*radius;
        dat.tcparams = { valinion valnasion vallear valrear valvertex 0 };

        %dat.imageaxis   = { -1  1 -1 };
        %dat.imageoffset = { [0.0  0.08 NaN ]  [0.01    NaN -0.01]  [ NaN -0.01   -0.025 ] };
        %dat.imagemult   = { 1.01 1.06 0.96 };
        %dat.axislim     = [-1.2 1.2 -1.2 1.2 -1.2 1.2];
        COLORMESH       = [.5 .5 .5];
        BACKCOLOR       = 'w';
    else 
        fid = fopen('VolumeMNI.bin', 'rb', 'ieee-le');
        if fid == -1
            error('Cannot find MRI data file');
        end;
        V = double(fread(fid, [108 129*129], 'uint8'))/255;
        V = reshape(V, 108, 129, 129);
        fclose(fid);
        
        %V = floatread('VolumeMNI.fdt');
        %load('/home/arno/matlab/MNI/VolumeMNI.mat');
        
        dat.imgs   = V; %smooth3(V,'gaussian', [3 3 3]);
        coordinc   = 2; % 2 mm
        allcoords1 = [0.5:coordinc:size(V,1)*coordinc]-size(V,1)/2*coordinc; 
        allcoords2 = [0.5:coordinc:size(V,2)*coordinc]-size(V,2)/2*coordinc;
        allcoords3 = [0.5:coordinc:size(V,3)*coordinc]-size(V,3)/2*coordinc;
        
        dat.imgcoords = { allcoords3        allcoords2        allcoords1 };
        if strcmpi(g.cornermri, 'on') % make the corner of the MRI volume match
            dat.maxcoord = [max(dat.imgcoords{1}) max(dat.imgcoords{2}) max(dat.imgcoords{3})];
        end;

        COLORMESH = 'w';
        BACKCOLOR = 'k';
        %valinion  = [ 58.5413  -10.5000  -30.8419 ]*2;
        %valnasion = [-56.8767  -10.5000  -30.9566 ]*2;
        %vallear   = [ 0.1040   -59.0000  -30.9000 ]*2;
        %valrear   = [ 0.1040    38.0000  -30.9000 ]*2;
        %valvertex = [ 0.0238   -10.5000   49.8341 ]*2;
        valinion  = [ 52.5413  -10.5000  -30.8419 ]*2;
        valnasion = [-50.8767  -10.5000  -30.9566 ]*2;
        vallear   = [ 0.1040   -51.0000  -30.9000 ]*2;
        valrear   = [ 0.1040    31.0000  -30.9000 ]*2;
        valvertex = [ 0.0238   -10.5000   40.8341 ]*2;
        zoffset   = 27.1190/(27.1190+radius) * (valvertex(3)-vallear(3));
        dat.tcparams = { valinion valnasion vallear valrear valvertex zoffset };
       
        %plotimgs(IMAGESLOC, IMAGESOFFSET, IMAGESMULT, IMAGESAXIS, AXISLIM, [57 85 65]);
        %view(30, 45); axis equal; return;
    end;

    % point 0
    % -------
    [xx yy zz] = transformcoords(0,0,0, dat.tcparams);
    dat.zeroloc = [ xx yy zz ];
    
    % conversion
    % ----------
    if strcmpi(g.normlen, 'on')
        if isfield(sources, 'besaextori')
            sources = rmfield(sources, 'besaextori'); 
        end;
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
    if ~isfield(sources, 'component')
        disp('No component indices, making incremental ones...');
        for index = 1:length(sources)
            sources(index).component = index;
        end;
    end;        

    % normalize position to unit sphere
    % ---------------------------------
    maxi = 0;
    for index = 1:length(sources)
        maxi = max(maxi,max(abs(sources(index).posxyz(:))));
    end;
    if maxi > 1.01 & g.sphere == 1
        disp('Non-normalized dipole positions, normalizing by standard head radius 85 mm'); 
        g.sphere = 85;
        fact = 0.1;
    else 
        fact = 1;
    end;
    
    % find non-empty sources
    % ----------------------
    noempt = cellfun('isempty', { sources.posxyz } );
    sources = sources( find(~noempt) );
    
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
            if sources(index).rv < g.rvrange(1)/100 | sources(index).rv > g.rvrange(2)/100
                sources(index) = [];
            end;
        end;
    end;
    
    % color array
    % -----------
    if isempty(g.color)
        g.color = { 'g' 'b' 'r' 'm' 'c' 'y' };
        if strcmp(BACKCOLOR, 'w'), g.color = { g.color{:} 'k' }; end;
    end;
    g.color = g.color(mod(0:length(sources)-1, length(g.color)) +1);
    if ~iscell(g.color)
        error('dipplot: ''color'' must be a cell array');
    end;
    
    % build summarized figure
    % -----------------------
    if strcmp(g.summary, 'on')
        figure;
        options = { 'gui', 'off', 'dipolesize', g.dipolesize/1.5,'dipolelength', g.dipolelength, 'sphere', g.sphere ...
                    'color', g.color, 'mesh', g.mesh, 'num', g.num, 'image', g.image 'normlen' g.normlen };
        axes('position', [0 0 0.5 0.5]);  dipplot(sourcesori, 'view', [1 0 0] , options{:}); axis off; if strcmpi(g.image, 'besa'), scalegca(0.1); end;
        axes('position', [0 0.5 0.5 .5]); dipplot(sourcesori, 'view', [0 0 1] , options{:}); axis off; if strcmpi(g.image, 'besa'), scalegca(0.1); end;
        axes('position', [.5 .5 0.5 .5]); dipplot(sourcesori, 'view', [0 -1 0], options{:}); axis off; if strcmpi(g.image, 'besa'), scalegca(0.1); end;
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
        
    % plot head graph in 3D
    % ---------------------
    if strcmp(g.gui, 'on')
        figure; 
        pos = get(gca, 'position');
        set(gca, 'position', [pos(1)+0.05 pos(2:end)]);
    end;
    plotimgs( dat, [1 1 1]);
    
    set(gca, 'color', BACKCOLOR);
    %warning off; a = imread('besaside.pcx'); warning on;
    % BECAUSE OF A BUG IN THE WARP FUNCTION, THIS DOES NOT WORK (11/02)
    %hold on; warp([], wy, wz, a);
    % set camera target 
    % -----------------
    
    % format axis (BESA or MRI)
    axis equal;
    set(gca, 'cameraviewanglemode', 'manual'); % disable change size
    camzoom(1.2^2);
    view(g.view);
    %set(gca, 'cameratarget',   dat.zeroloc); % disable change size
    %set(gca, 'cameraposition', dat.zeroloc+g.view*g.zoom); % disable change size
    axis off;
        
    % plot sphere mesh and nose
    % -------------------------
    [x y z] = sphere(20);
    hold on; 
    [xx yy zz] = transformcoords(x, y, z, dat.tcparams);
    if strcmpi(COLORMESH, 'w')
        hh = mesh(xx, yy, zz, 'cdata', ones(21,21,3), 'tag', 'mesh'); hidden off;
    else
        hh = mesh(xx, yy, zz, 'cdata', zeros(21,21,3), 'tag', 'mesh'); hidden off;
    end;
    %x = x*100*scaling; y = y*100*scaling; z=z*100*scaling;
    %h = line(xx,yy,zz); set(h, 'color', COLORMESH, 'linestyle', '--', 'tag', 'mesh');
    %h = line(xx,zz,yy); set(h, 'color', COLORMESH, 'linestyle', '--', 'tag', 'mesh');
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
        
        % reorder dipoles for plotting
        if nbdip == 2
            if sources(index).posxyz(1,1) > sources(index).posxyz(2,1)
                tmp = sources(index).posxyz(2,:);
                sources(index).posxyz(2,:) = sources(index).posxyz(1,:);
                sources(index).posxyz(1,:) = tmp;
                tmp = sources(index).momxyz(2,:);
                sources(index).momxyz(2,:) = sources(index).momxyz(1,:);
                sources(index).momxyz(1,:) = tmp;
            end;
        end;
        
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
                xo1 = x+xo;
                yo1 = y+yo;
                zo1 = z+zo;
            elseif strcmpi(g.pointout,'on')
                xo1 = x-xo; % make dipole point outward from head center
                yo1 = y-yo;
                zo1 = z-zo; 
            else
                xo1 = x+xo;
                yo1 = y+yo;
                zo1 = z+zo;
            end
            x = -x; xo1 = -xo1; 
            y = -y; yo1 = -yo1; 
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% draw dipole bar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tag = [ 'dipole' num2str(index) ];
            [xx   yy   zz]   = transformcoords(x,   y,   z,   dat.tcparams); 
            [xxo1 yyo1 zzo1] = transformcoords(xo1, yo1, zo1,  dat.tcparams); 
            h1 = line( [xx xxo1]', [yy yyo1]', [zz zzo1]' );
            dipstruct.pos3d  = [xx yy zz];            % value used for fitting MRI
            dipstruct.posxyz = sources(index).posxyz; % value used for other purposes
            dipstruct.rv     = sprintf('C %d (%3.2f)', sources(index).component, sources(index).rv*100);
            set(h1, 'userdata', dipstruct, 'tag', tag, 'color','k', 'linewidth', g.dipolesize/7.5);
            if strcmp(BACKCOLOR, 'k'), set(h1, 'color', g.color{index}); end;
            
            % draw point
            hold on;
            h = plot3(xx,  yy,  zz); 
            set(h, 'userdata', dipstruct, 'tag', tag, ...
                   'marker', '.', 'markersize', g.dipolesize, 'color', g.color{index});
            
            % project onto images
            % -------------------
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
                    h = line( [xx xxo1]', [yy yyo1]', [-1 -1]'*dat.maxcoord(1));
                    set(h, 'userdata', 'proj', 'tag', tag, 'color','k', 'linewidth', g.dipolesize/7.5);
                end;
                if strcmp(BACKCOLOR, 'k'), set(h, 'color', tmpcolor); end;
                h = plot3(xx,  yy,  -dat.maxcoord(1)); 
                set(h, 'userdata', 'proj', 'tag', tag, ...
                       'marker', '.', 'markersize', g.dipolesize, 'color', tmpcolor);
                
                % project onto x axis
                tag = [ 'dipole' num2str(index) ];
                if ~strcmpi(g.image, 'besa')
                    h = line( [xx xxo1]', [1 1]'*dat.maxcoord(2), [zz zzo1]');
                    set(h, 'userdata', 'proj', 'tag', tag, 'color','k', 'linewidth', g.dipolesize/7.5);
                end;
                if strcmp(BACKCOLOR, 'k'), set(h, 'color', tmpcolor); end;
                h = plot3(xx,  dat.maxcoord(2),  zz); 
                set(h, 'userdata', 'proj', 'tag', tag, ...
                       'marker', '.', 'markersize', g.dipolesize, 'color', tmpcolor);
                
                % project onto y axis
                tag = [ 'dipole' num2str(index) ];
                if ~strcmpi(g.image, 'besa')
                    h = line( [-1 -1]'*dat.maxcoord(3), [yy yyo1]', [zz zzo1]');
                    set(h, 'userdata', 'proj', 'tag', tag, 'color','k', 'linewidth', g.dipolesize/7.5);
                end;
                if strcmp(BACKCOLOR, 'k'), set(h, 'color', tmpcolor); end;
                h = plot3(-dat.maxcoord(3),  yy,  zz); 
                set(h, 'userdata', 'proj', 'tag', tag, ...
                       'marker', '.', 'markersize', g.dipolesize, 'color', tmpcolor);
            end;

            % line projections
            % -------------------
            if strcmpi(g.projlines, 'on')
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
                tmpcolor = g.color{index};
                
                % project onto z axis
                tag = [ 'dipole' num2str(index) ];
                h(1) = line( [xx xx]', [yy yy]', [zz -dat.maxcoord(1)]');
                set(h(1), 'userdata', 'proj', 'linestyle', '--', 'tag', tag, 'color', tmpcolor, 'linewidth', g.dipolesize/7.5/5);
                
                % project onto x axis
                tag = [ 'dipole' num2str(index) ];
                h(2) = line( [xx xx]', [yy dat.maxcoord(2)]', [zz zz]');
                set(h(2), 'userdata', 'proj', 'linestyle', '--', 'tag', tag, 'color', tmpcolor, 'linewidth', g.dipolesize/7.5/5);
                
                % project onto y axis
                tag = [ 'dipole' num2str(index) ];
                h(3) = line( [xx -dat.maxcoord(3)]', [yy yy]', [zz zz]');
                set(h(3), 'userdata', 'proj', 'linestyle', '--', 'tag', tag, 'color', tmpcolor, 'linewidth', g.dipolesize/7.5/5);
                if ~isempty(g.projcol)
                    set(h, 'color', g.projcol);
                end;
            end;
                        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% draw text  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if isfield(sources, 'component')
                if strcmp(g.num, 'on')
                    h = text(xx,  yy,  zz, [ '  ' int2str(sources(index).component)]);
                    set(h, 'userdata', dipstruct, 'tag', tag, 'fontsize', g.dipolesize/2 );
                    if ~strcmpi(g.image, 'besa'), set(h, 'color', 'w'); end;
                end;
            end;
        end;
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% draw elipse for group of dipoles  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~isempty(g.std)
        for index = 1:length(g.std)
            if ~iscell(g.std{index})
                plotellipse(sources, g.std{index}, 1, dat.tcparams);
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
    cbmesh = [ 'if get(gcbo, ''userdata''), ' ...
               '    set(findobj(''parent'', gca, ''tag'', ''mesh''), ''visible'', ''off'');' ...
               '    set(gcbo, ''string'', ''Mesh on'');' ...
               '    set(gcbo, ''userdata'', 0);' ...
               'else,' ...
               '    set(findobj(''parent'', gca, ''tag'', ''mesh''), ''visible'', ''on'');' ...
               '    set(gcbo, ''string'', ''Mesh off'');' ...
               '    set(gcbo, ''userdata'', 1);' ...
               'end;' ];
    cbview = [ 'tmpuserdat = get(gca, ''userdata'');' ...
             'if tmpuserdat.axistight, ' ...
             '    set(gcbo, ''string'', ''Tight view'');' ...
             'else,' ...
             '    set(gcbo, ''string'', ''Loose view'');' ...
             'end;' ...
             'tmpuserdat.axistight = ~tmpuserdat.axistight;' ...
             'set(gca, ''userdata'', tmpuserdat);' ...
             'clear tmpuserdat;' ...
             'dipplot(gcbf);' ];
    viewstyle  = fastif(strcmpi(dat.mode, 'besa'), 'text', 'pushbutton');
    viewstring = fastif(dat.axistight, 'Loose view', 'Tight view');
    
    h = uicontrol( 'unit', 'normalized', 'position', [0 0 .15 .05], 'tag', 'tmp', ...
                  'style', 'pushbutton', 'string', 'Top view', 'callback', 'view([0 0 1]);');
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.05 .15 .05], 'tag', 'tmp', ...
                  'style', 'pushbutton', 'string', 'Coronal view', 'callback', 'view([0 -1 0]);');
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.1 .15 .05], 'tag', 'tmp', ...
                  'style', 'pushbutton', 'string', 'Sagital view', 'callback', 'view([1 0 0]);');
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.15 .15 .05], 'tag', 'tmp', ...
                  'style', 'text', 'string', '');
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.2  .15 .05], 'tag', 'tmp', ...
                  'style',  viewstyle   , 'string', viewstring, 'callback', cbview);
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.25 .15 .05], 'tag', 'tmp', ...
                  'style', 'pushbutton', 'string', 'Mesh on', 'userdata', 0, 'callback', cbmesh);
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
                  'dipplot(gcbf);' ...
                  'clear editobj;' ]);
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.65 .15 .05], 'tag', 'tmp', ...
                  'style', 'pushbutton', 'string', 'Keep|Next', 'callback', ...
                [ 'editobj = findobj(''parent'', gcf, ''userdata'', ''editor'');' ...
                  'set(editobj, ''string'', num2str(str2num(get(editobj, ''string''))+1));' ...
                  'tmpobj = get(gcf, ''userdata'');' ...
                  'dipplot(gcbf);' ...
                  'set(tmpobj, ''visible'', ''on'');' ...
                  'clear editobj tmpobj;' ]);
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.7 .15 .05], 'tag', 'tmp', 'userdata', 'rv', ...
                  'style', 'text', 'string', '');
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.75 .15 .05], 'tag', 'tmp', 'userdata', 'editor', ...
                  'style', 'edit', 'string', '1', 'callback', ...
                  [ 'dipplot(gcbf);' ] );
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.8 .15 .05], 'tag', 'tmp', ...
                   'style', 'pushbutton', 'string', 'Plot one', 'callback', ...
             	    [ 'for tmpi = 1:' nbsrc ',' ...
                   '   set(findobj(''parent'', gca, ''tag'', [ ''dipole'' int2str(tmpi) ]), ''visible'', ''off'');' ...
                   'end; clear tmpi;' ...
                   'dipplot(gcbf);' ]);
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

    dat.nbsources  = length(sources);
    set(gca, 'userdata', dat ); % last param=1 for MRI view tight/loose
    set(gcf, 'color', BACKCOLOR);
        
    if strcmp(g.gui, 'off')
        set(findobj('parent', gcf, 'tag', 'tmp'), 'visible', 'off');
    end;
    if strcmp(g.mesh, 'off')
        set(findobj('parent', gca, 'tag', 'mesh'), 'visible', 'off');
    end;
    updatedipplot(gcf);
    
    rotate3d on;
return;

function [x,y,z] = transformcoords(x, y, z, TCPARAMS);
    if iscell(TCPARAMS)
        
        valinion  = TCPARAMS{1};
        valnasion = TCPARAMS{2};
        vallear   = TCPARAMS{3};
        valrear   = TCPARAMS{4};
        valvertex = TCPARAMS{5};
        zoffset   = TCPARAMS{6};
        
        % scale axis
        % ----------
        y = y/2*( valinion(1)  - valnasion(1) );
        x = x/2*( valrear(2)   - vallear(2) );
        z = z*  ( valvertex(3) - vallear(3) - zoffset );
        
        % recenter
        % --------
        y = y + (vallear (1) + valrear  (1))/2;
        x = x + (valinion(2) + valnasion(2))/2;
        z = z + (vallear (3) + valrear  (3))/2 + zoffset;
        return;

        % rotate
        % ------
        %% pitch (x-axis); roll = y axis rotation; yaw = z axis
        %% see http://bishopw.loni.ucla.edu/AIR5/homogenous.html
        %pitch = atan2(vallear(1)  -valrear(1),  vallear(3)-valrear(3));    fprintf('Pitch: %1.1g degree\n', pitch/pi*180);
        %roll  = atan2(valnasion(1)-valinion(1), valnasion(3)-valinion(3)); fprintf('Roll:  %1.1g degree\n', roll/pi*180);
        %yaw   = 0;
        %cp = cos(pitch); sp = sin(pitch);
        %cr = cos(roll);  sr = sin(roll);
        %cy = cos(yaw);   sy = sin(yaw);
        
        %rot3d = [ cy*cr+sy*sp*sr    sy*cr-cy*sp*sr     cp*sr  ;
        %          -sy*cp            cy*cp               sp     ;
        %          sy*sp*cr-cy*sr    -cy*sp*cr-sy*sr     cp*cr  ];
        %elec = rot3d*elec; % ROTATE NOT NECESSARY FOR MNI HEAD
    elseif isnumeric(TCPARAMS)
        y = y*TCPARAMS;
        x = x*TCPARAMS;
        z = z*TCPARAMS;        
    end;
    
function sc = plotellipse(sources, ind, nstd, TCPARAMS);

    for i = 1:length(ind)
        tmpval(1,i) = -sources(ind(i)).posxyz(1);    
        tmpval(2,i) = -sources(ind(i)).posxyz(2);    
        tmpval(3,i) = sources(ind(i)).posxyz(3);
        [tmpval(1,i) tmpval(2,i) tmpval(3,i)] = transformcoords(tmpval(1,i), tmpval(2,i), tmpval(3,i), TCPARAMS);
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
function updatedipplot(fig)
   
   % find current dipole index and test for authorized range
   % -------------------------------------------------------
   dat     = get(gca, 'userdata');
   editobj = findobj('parent', fig, 'userdata', 'editor');
   tmpnum  = str2num(get(editobj(end), 'string'));
   if tmpnum < 1,             tmpnum = 1;             end;
   if tmpnum > dat.nbsources, tmpnum = dat.nbsources; end;
   set(editobj(end), 'string', num2str(tmpnum));
   
   % hide current dipole, find next dipole and show it
   % -------------------------------------------------
   set(get(gcf, 'userdata'), 'visible', 'off');
   newdip = findobj('parent', gca, 'tag', [ 'dipole' get(editobj(end), 'string')]);
   set(newdip, 'visible', 'on');
   set(gcf, 'userdata', newdip);
   
   % set the new RV
   % --------------
   tmprvobj = findobj('parent', fig, 'userdata', 'rv');
   index = 1;
   while ~isstruct( get(newdip(index), 'userdata') )
       index = index+1;
   end;
   userdat  = get(newdip(index), 'userdata');
   
   set( tmprvobj(end), 'string', userdat.rv);
   
   % adapt the MRI to the dipole depth
   % ---------------------------------
   %if ~strcmpi(dat.mode, 'besa') % not besa mode
   delete(findobj('parent', gca, 'tag', 'img'));
   
   tmpdiv = 1;
   if ~dat.axistight
       [xx yy zz] = transformcoords(0,0,0, dat.tcparams);
       indx = minpos(dat.imgcoords{1}-zz);
       indy = minpos(dat.imgcoords{2}-yy);
       indz = minpos(dat.imgcoords{3}-xx);
   else
       indx = minpos(dat.imgcoords{1} - userdat.pos3d(3) + 4/tmpdiv);
       indy = minpos(dat.imgcoords{2} - userdat.pos3d(2) - 4/tmpdiv);
       indz = minpos(dat.imgcoords{3} - userdat.pos3d(1) + 4/tmpdiv);
   end;
   plotimgs( dat, [indx indy indz]);
   %end;
   	
% plot images
% -----------
function plotimgs(dat, index);
   
    % loading images
    % --------------
    if strcmpi(dat.mode, 'besa')
        imgt = dat.imgs{1};
        imgc = dat.imgs{2};
        imgs = dat.imgs{3};
    else
        imgt = rot90(squeeze(dat.imgs(:,:,index(1)))); 
        imgc = rot90(squeeze(dat.imgs(:,index(2),:))); 
        imgs = rot90(squeeze(dat.imgs(index(3),:,:))); 
    end;
    if ndims(imgt) == 2, imgt(:,:,3) = imgt; imgt(:,:,2) = imgt(:,:,1); end;
    if ndims(imgc) == 2, imgc(:,:,3) = imgc; imgc(:,:,2) = imgc(:,:,1); end;
    if ndims(imgs) == 2, imgs(:,:,3) = imgs; imgs(:,:,2) = imgs(:,:,1); end;
    
    % computing coordinates for planes
    % --------------------------------    
    wxt = [min(dat.imgcoords{3}) max(dat.imgcoords{3}); min(dat.imgcoords{3}) max(dat.imgcoords{3})];
    wyt = [min(dat.imgcoords{2}) min(dat.imgcoords{2}); max(dat.imgcoords{2}) max(dat.imgcoords{2})];
    wxc = [min(dat.imgcoords{3}) max(dat.imgcoords{3}); min(dat.imgcoords{3}) max(dat.imgcoords{3})];
    wzc = [max(dat.imgcoords{1}) max(dat.imgcoords{1}); min(dat.imgcoords{1}) min(dat.imgcoords{1})];
    wys = [min(dat.imgcoords{2}) max(dat.imgcoords{2}); min(dat.imgcoords{2}) max(dat.imgcoords{2})];
    wzs = [max(dat.imgcoords{1}) max(dat.imgcoords{1}); min(dat.imgcoords{1}) min(dat.imgcoords{1})];
    if dat.axistight
        wzt = [ 1 1; 1 1]*dat.imgcoords{1}(index(1));
        wyc = [ 1 1; 1 1]*dat.imgcoords{2}(index(2));
        wxs = [ 1 1; 1 1]*dat.imgcoords{3}(index(3));
    else
        wzt = -[ 1 1; 1 1]*dat.maxcoord(1);
        wyc =  [ 1 1; 1 1]*dat.maxcoord(2);
        wxs = -[ 1 1; 1 1]*dat.maxcoord(3);
    end;
    
    % ploting surfaces
    % ----------------
    options = { 'FaceColor','texturemap', 'EdgeColor','none', 'CDataMapping', ...
                'direct','tag','img', 'facelighting', 'none' };
    hold on;
    surface(wxt, wyt, wzt, imgt(end:-1:1,:,:), options{:});
    surface(wxc, wyc, wzc, imgc              , options{:});
    surface(wxs, wys, wzs, imgs              , options{:});
    if strcmpi(dat.drawedges, 'on')
        % removing old edges if any
        delete(findobj( gcf, 'tag', 'edges'));
        if dat.axistight, col = 'k'; else col = [0.5 0.5 0.5]; end;
        h(1) = line([wxs(1) wxs(2)]', [wys(1) wyc(1)]', [wzt(1) wzt(2)]'); % sagital-transverse
        h(2) = line([wxs(1) wxc(3)]', [wyc(1) wyc(2)]', [wzt(1) wzt(2)]'); % coronal-tranverse
        h(3) = line([wxs(1) wxs(2)]', [wyc(1) wyc(2)]', [wzs(1) wzt(1)]'); % sagital-coronal
        set(h, 'color', col, 'linewidth', 2, 'tag', 'edges');
    end;
        
    %%fill3([-2 -2 2 2], [-2 2 2 -2], wz(:)-1, BACKCOLOR);
    %%fill3([-2 -2 2 2], wy(:)-1, [-2 2 2 -2], BACKCOLOR);
    rotate3d on 

function index = minpos(vals);
	vals(find(vals < 0)) = inf;
	[tmp index] = min(vals);

function scalegca(factor)
    xl = xlim; xf = ( xl(2) - xl(1) ) * factor;
    yl = ylim; yf = ( yl(2) - yl(1) ) * factor;
    zl = zlim; zf = ( zl(2) - zl(1) ) * factor;
    xlim( [ xl(1)-xf xl(2)+xf ]);
    ylim( [ yl(1)-yf yl(2)+yf ]);
    zlim( [ zl(1)-zf zl(2)+zf ]);
    