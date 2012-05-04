% dipplot() - Visualize EEG equivalent-dipole locations and orientations 
%             in the MNI average MRI head or in the BESA spherical head model. 
% Usage:
%   >> dipplot( sources, 'key', 'val', ...);
%   >> [sources X Y Z XE YE ZE] = dipplot( sources, 'key', 'val', ...);
%
% Inputs:
%   sources   -  structure array of dipole information: can contain
%                either BESA or DIPFIT dipole information. BESA dipole
%                information are still supported but may disapear in the 
%                future. For DIPFIT
%                   sources.posxyz: contains 3-D location of dipole in each 
%                                   column. 2 rows indicate 2 dipoles.
%                   sources.momxyz: contains 3-D moments for dipoles above.
%                   sources.rv    : residual variance from 0 to 1.
%                   other fields  : used for graphic interface.
%
% Optional input:
%  'rvrange'  - [min max] or [max] Only plot dipoles with residual variace
%               within the given range. Default: plot all dipoles.
%  'summary'  - ['on'|'off'|'3d'] Build a summary plot with three views (top, 
%               back, side). {default: 'off'}
%  'mri'      - Matlab file containing an MRI volume and a 4-D transformation
%               matrix to go from voxel space to electrode space:
%               mri.anatomy   contains a 3-D anatomical data array
%               mri.transfrom contains a 4-D homogenous transformation matrix.
%  'coordformat' - ['MNI'|'spherical'] Consider that dipole coordinates are in
%               MNI or spherical coordinates (for spherical, the radius of the 
%               head is assumed to be 85 (mm)). See also function sph2spm().
%  'transform' - [real array] traditional transformation matrix to convert
%               dipole coordinates to MNI space. Default is assumed from 
%               'coordformat' input above. Type help traditional for more 
%               information.
%  'image'    - ['besa'|'mri'] Background image. 
%               'mri' (or 'fullmri') uses mean-MRI brain images from the Montreal 
%               Neurological Institute. This option can also contain a 3-D MRI
%               volume (dim 1: left to right; dim 2: anterior-posterior; dim 3: 
%               superior-inferior). Use 'coregist' to coregister electrodes
%               with the MRI. {default: 'mri'} 
%  'verbose' - ['on'|'off'] comment on operations on command line {default:
%  'on'}.
%  'plot'    - ['on'|'off'] only return outputs                  {default: 'off'}.
%
% Plotting options:
%  'color'    - [cell array of color strings or (1,3) color arrays]. For
%               exemple { 'b' 'g' [1 0 0] } gives blue, green and red. 
%               Dipole colors will rotate through the given colors if
%               the number given is less than the number of dipoles to plot.
%               A single number will be used as color index in the jet colormap.
%  'view'     - 3-D viewing angle in cartesian coords.,
%               [0 0 1] gives a sagittal view, [0 -1 0] a view from the rear;
%               [1 0 0] gives a view from the side of the head.
%  'mesh'     - ['on'|'off'] Display spherical mesh. {Default is 'on'}
%  'meshdata' - [cell array|'file_name'] Mesh data in a cell array { 'vertices'
%               data 'faces' data } or a boundary element model filename (the 
%               function will plot the 3rd mesh in the 'bnd' sub-structure).
%  'axistight' - ['on'|'off'] For MRI only, display the closest MRI
%               slide. {Default is 'off'}
%  'gui'      - ['on'|'off'] Display controls. {Default is 'on'} If gui 'off', 
%               a new figure is not created. Useful for incomporating a dipplot 
%               into a complex figure.
%  'num'      - ['on'|'off'] Display component number. Take into account
%               dipole size. {Default: 'off'}
%  'cornermri' - ['on'|'off'] force MRI images to the corner of the MRI volume
%               (usefull when background is not black). Default: 'off'.
%  'drawedges' - ['on'|'off'] draw edges of the 3-D MRI (black in axistight,
%                white otherwise.) Default is 'off'.
%  'projimg'  - ['on'|'off'] Project dipole(s) onto the 2-D images, for use
%               in making 3-D plots {Default 'off'}
%  'projlines' - ['on'|'off'] Plot lines connecting dipole with 2-D projection.
%                Color is dashed black for BESA head and dashed black for the
%                MNI brain {Default 'off'}
%  'projcol'   - [color] color for the projected line {Default is same as dipole}
%  'dipolesize' - Size of the dipole sphere(s). This option may also contain one
%               value per dipole {Default: 30}
%  'dipolelength' - Length of the dipole bar(s) {Default: 1}
%  'pointout' - ['on'|'off'] Point the dipoles outward. {Default: 'off'}
%  'sphere'   - [float] radius of sphere corresponding to the skin. Default is 1.
%  'spheres'  - ['on'|'off'] {default: 'off'} plot dipole markers as 3-D spheres. 
%               Does not yet interact with gui buttons, produces non-gui mode.
%  'spheresize' - [real>0] size of spheres (if 'on'). {default: 5}
%  'normlen'  - ['on'|'off'] Normalize length of all dipoles. {Default: 'off'}
%  'dipnames' - [cell array] cell array of string with a name for each dipole (or
%               pair of dipole).
%  'holdon'   - ['on'|'off'] create a new dipplot figure or plot dipoles within an
%               an existing figure. Default is 'off'.
%  'camera'   - ['auto'|'set'] camera position. 'auto' is the default and 
%               an option using camera zoom. 'set' is a fixed view that
%               does not depend on the content being plotted.
%
% Outputs:
%   sources   - EEG.source structure with two extra fiels 'mnicoord' and 'talcoord'
%               containing the MNI and talairach coordinates of the dipoles. Note
%               that for the BEM model, dipoles are already in MNI coordinates.
%   X,Y,Z     - Locations of dipole heads (Cartesian coordinates in MNI space). 
%               If there is more than one dipole per components, the last dipole 
%               is returned.
%   XE,YE,ZE  - Locations of dipole ends (Cartesian coordinates). The same
%               remark as above applies.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 1st July 2002
%
% Notes: See DIPFIT web tutorial at sccn.ucsd.edu/eeglab/dipfittut/dipfit.html
%        for more details about MRI co-registration etc...
%
% Example:
%  % define dipoles
%  sources(1).posxyz = [-59 48 -28];   % position for the first dipole
%  sources(1).momxyz = [  0 58 -69];   % orientation for the first dipole
%  sources(1).rv     = 0.036;          % residual variance for the first dipole
%  sources(2).posxyz = [74 -4 -38];    % position for the second dipole
%  sources(2).momxyz = [43 -38 -16]; % orientation for the second dipole
%  sources(2).rv     = 0.027;          % residual variance for the second dipole
%
%  % plot of the two dipoles (first in green, second in blue)
%  dipplot( sources, 'color', { 'g' 'b' }); 
%
%  % To make a stereographic plot
%  figure( 'position', [153 553 1067 421]; 
%  subplot(1,3,1); dipplot( sources, 'view', [43 10], 'gui', 'off');
%  subplot(1,3,3); dipplot( sources, 'view', [37 10], 'gui', 'off');
%
%  % To make a summary plot
%  dipplot( sources, 'summary', 'on', 'num', 'on');
%
% See also: eeglab(), dipfit()

% old options
% -----------
%  'std'      - [cell array] plot standard deviation of dipoles. i.e.
%               { [1:6] [7:12] } plot two elipsoids that best fit all the dipoles
%               from 1 to 6 and 7 to 12 with radius 1 standard deviation.
%               { { [1:6] 2 'linewidth' 2 } [7:12] } do the same but now the
%               first elipsoid is 2 standard-dev and the lines are thicker.

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
    g = finputcheck( varargin, { 'color'     ''         []                  [];
                                 'axistight' 'string'   { 'on' 'off' }     'off';
                                 'camera'    'string'   { 'auto' 'set' }   'auto';
                                 'coordformat' 'string' { 'MNI' 'spherical' 'CTF' 'auto' } 'auto';
                                 'drawedges' 'string'   { 'on' 'off' }     'off';
                                 'mesh'      'string'   { 'on' 'off' }     'off';
                                 'gui'       'string'   { 'on' 'off' }     'on';
                                 'summary'   'string'   { 'on2' 'on' 'off' '3d' }     'off';
                                 'verbose'   'string'   { 'on' 'off' }     'on';
                                 'view'      'real'     []                 [0 0 1];
                                 'rvrange'   'real'     [0 Inf]             [];
                                 'transform' 'real'     [0 Inf]             [];
                                 'normlen'   'string'   { 'on' 'off' }     'off';
                                 'num'       'string'   { 'on' 'off' }     'off';
                                 'cornermri' 'string'   { 'on' 'off' }     'off';
                                 'mri'       { 'string' 'struct' } []      '';
                                 'dipnames'   'cell'     []                  {};
                                 'projimg'   'string'   { 'on' 'off' }     'off';
                                 'projcol'   ''         []       [];
                                 'projlines' 'string'   { 'on' 'off' }     'off';
                                 'pointout'  'string'   { 'on' 'off' }     'off';
                                 'holdon'    'string'   { 'on' 'off' }     'off';
                                 'dipolesize' 'real'    [0 Inf]             30;
                                 'dipolelength' 'real'  [0 Inf]             1;
                                 'sphere'    'real'   [0 Inf]               1;
                                 'spheres'   'string'  {'on' 'off'}         'off';
                                 'links'     'real'   []                    [];
                                 'image'     { 'string' 'real' } []         'mri';
                                 'plot'      'string'   { 'on' 'off' }      'on';
                                 'meshdata'  { 'string' 'cell' } []         '' }, 'dipplot');
    %                             'std'       'cell'     []                  {}; 
    %                             'coreg'     'real'     []                  [];

    if isstr(g), error(g); end;
    if strcmpi(g.holdon, 'on'), g.gui = 'off'; end;
    if length(g.dipolesize) == 1, g.dipolesize = repmat(g.dipolesize, [1 length(sourcesori)]); end;
    
    g.zoom = 1500;

    if strcmpi(g.image, 'besa')
        error('BESA image not supported any more. Use EEGLAB version 4.512 or earlier. (BESA dipoles can still be plotted in MNI brain.)');
    end;
    
    % trying to determine coordformat
    % -------------------------------
    if ~isfield(sources, 'momxyz')
        g.coordformat = 'spherical';
    end;
    if strcmpi(g.coordformat, 'auto')
        if ~isempty(g.meshdata)
            g.coordformat = 'MNI';
            if strcmpi(g.verbose, 'on'), 
                disp('Coordinate format unknown: using ''MNI'' since mesh data was provided as input');
            end
        else
            maxdiplen = 0;
            for ind = 1:length(sourcesori)
                maxdiplen = max(maxdiplen, max(abs(sourcesori(ind).momxyz(:))));
            end;
            if maxdiplen>2000
                if strcmpi(g.verbose, 'on'),
                    disp('Coordinate format unknown: using ''MNI'' because of large dipole moments');
                end
            else
                g.coordformat = 'spherical';
                if strcmpi(g.verbose, 'on'),
                    disp('Coordinate format unknown: using ''spherical'' since no mesh data was provided as input');
                end
            end;
        end;
    end;
    
    % axis image and limits
    % ---------------------
    dat.axistight  = strcmpi(g.axistight, 'on');
    dat.drawedges  = g.drawedges;
    dat.cornermri  = strcmpi(g.cornermri, 'on');
    radius = 85;

    % look up an MRI file if necessary
    % --------------------------------
    if isempty(g.mri)
        if strcmpi(g.verbose, 'on'),
            disp('No MRI file given as input. Looking up one.');
        end
        dipfitdefs;
        g.mri = template_models(1).mrifile;
    end;
        
    % read anatomical MRI using Fieldtrip and SPM2 functons
    % -----------------------------------------------------
    if isstr(g.mri);
        try, 
            g.mri = load('-mat', g.mri);
            g.mri = g.mri.mri;
        catch,
            disp('Failed to read Matlab file. Attempt to read MRI file using function read_fcdc_mri');
            try,
                warning off;
                g.mri = read_fcdc_mri(g.mri);
                %g.mri.anatomy(find(g.mri.anatomy > 255)) = 255;
                %g.mri.anatomy = uint8(g.mri.anatomy);
                g.mri.anatomy = round(gammacorrection( g.mri.anatomy, 0.8));
                g.mri.anatomy = uint8(round(g.mri.anatomy/max(reshape(g.mri.anatomy, prod(g.mri.dim),1))*255));
                % WARNING: if using double instead of int8, the scaling is different 
                % [-128 to 128 and 0 is not good]
                % WARNING: the transform matrix is not 1, 1, 1 on the diagonal, some slices may be 
                % misplaced
                warning on;
            catch,
                error('Cannot load file using read_fcdc_mri');
            end;
        end;
    end;
    
    if strcmpi(g.coordformat, 'spherical')
         dat.sph2spm    = sph2spm;
    elseif strcmpi(g.coordformat, 'CTF')
        dat.sph2spm    = traditionaldipfit([0 0 0 0 0 0 10 -10 10]);
    else
        dat.sph2spm    = []; %traditional([0 0 0 0 0 pi 1 1 1]);
    end;
    
    if ~isempty(g.transform), dat.sph2spm = traditionaldipfit(g.transform);
    end;
    if isfield(g.mri, 'anatomycol')
        dat.imgs       = g.mri.anatomycol;
    else
        dat.imgs       = g.mri.anatomy;
    end;
    dat.transform  = g.mri.transform;    
    
    % MRI coordinates for slices
    % --------------------------
    if ~isfield(g.mri, 'xgrid')
        g.mri.xgrid = [1:size(dat.imgs,1)]; 
        g.mri.ygrid = [1:size(dat.imgs,2)];
        g.mri.zgrid = [1:size(dat.imgs,3)];
    end;
    if strcmpi(g.coordformat, 'CTF')
        g.mri.zgrid = g.mri.zgrid(end:-1:1);
    end;

    dat.imgcoords = { g.mri.xgrid g.mri.ygrid g.mri.zgrid };            
    dat.maxcoord  = [max(dat.imgcoords{1}) max(dat.imgcoords{2}) max(dat.imgcoords{3})];
    COLORMESH = 'w';
    BACKCOLOR = 'k';
    
    % point 0
    % -------
    [xx yy zz] = transform(0, 0, 0, dat.sph2spm); % nothing happens for BEM since dat.sph2spm is empty
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
        if strcmpi(g.verbose, 'on'),
            disp('No component indices, making incremental ones...');
        end        
        for index = 1:length(sources)
            sources(index).component = index;
        end;
    end;        
    
    % find non-empty sources
    % ----------------------
    noempt = cellfun('isempty', { sources.posxyz } );
    sources = sources( find(~noempt) );
    
    % transform coordinates
    % ---------------------
    outsources = sources;
    for index = 1:length(sources)
        sources(index).momxyz = sources(index).momxyz/1000;
    end;
    
    % remove 0 second dipoles if any
    % ------------------------------
    for index = 1:length(sources)
        if size(sources(index).momxyz,1) == 2
            if all(sources(index).momxyz(2,:) == 0)
                sources(index).momxyz = sources(index).momxyz(1,:);
                sources(index).posxyz = sources(index).posxyz(1,:);
            end;
        end;
    end;

    % remove sources with out of bound Residual variance
    % --------------------------------------------------
    if isfield(sources, 'rv') & ~isempty(g.rvrange)
        if length(g.rvrange) == 1, g.rvrange = [ 0 g.rvrange ]; end;
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
    if ~isempty(g.color)
        g.color = strcol2real( g.color, jet(64) );
    end;
    if ~isempty(g.projcol)
        g.projcol = strcol2real( g.projcol, jet(64) );
        g.projcol = g.projcol(mod(0:length(sources)-1, length(g.projcol)) +1);        
    else
        g.projcol = g.color;
        for index = 1:length(g.color)
            g.projcol{index} =  g.projcol{index}/2;
        end;
    end;
    
    % build summarized figure
    % -----------------------
    if strcmpi(g.summary, 'on') | strcmpi(g.summary, 'on2')
        figure;
        options = { 'gui', 'off', 'dipolesize', g.dipolesize/1.5,'dipolelength', g.dipolelength, 'sphere', g.sphere ...
                    'color', g.color, 'mesh', g.mesh, 'num', g.num, 'image', g.image 'normlen' g.normlen ...
                    'coordformat' g.coordformat 'mri' g.mri 'meshdata' g.meshdata 'axistight' g.axistight };
        pos1 = [0 0 0.5 0.5];
        pos2 = [0 0.5 0.5 .5];
        pos3 = [.5 .5 0.5 .5]; if strcmp(g.summary, 'on2'), tmp = pos1; pos1 =pos3; pos3 = tmp; end;
        axes('position', pos1);  newsources = dipplot(sourcesori, 'view', [1 0 0] , options{:}); axis off; 
        axes('position', pos2); newsources = dipplot(sourcesori, 'view', [0 0 1] , options{:}); axis off; 
        axes('position', pos3); newsources = dipplot(sourcesori, 'view', [0 -1 0], options{:}); axis off; 
        axes('position', [0.5 0 0.5 0.5]); 
        colorcount = 1;
        if isfield(newsources, 'component')
            for index = 1:length(newsources)
                if isempty(g.dipnames), tmpname = sprintf( 'Comp. %d', newsources(index).component);
                else                    tmpname = char(g.dipnames{index});
                end;
                talpos = newsources(index).talcoord;
                if strcmpi(g.coordformat, 'CTF')
                    textforgui(colorcount) = { sprintf( [ tmpname ' (RV:%3.2f%%)' ], 100*newsources(index).rv) };
                elseif size(talpos,1) == 1
                    textforgui(colorcount) = { sprintf( [ tmpname ' (RV:%3.2f%%; Tal:%d,%d,%d)' ], ...
                                                        100*newsources(index).rv, ...
                                                        round(talpos(1,1)), round(talpos(1,2)), round(talpos(1,3))) };
                else
                    textforgui(colorcount) = { sprintf( [ tmpname ' (RV:%3.2f%%; Tal:%d,%d,%d & %d,%d,%d)' ], ...
                                                        100*newsources(index).rv, ...
                                                        round(talpos(1,1)), round(talpos(1,2)), round(talpos(1,3)), ...
                                                        round(talpos(2,1)), round(talpos(2,2)), round(talpos(2,3))) };
                end;
                colorcount = colorcount+1;
            end;
            colorcount = colorcount-1;
            allstr = strvcat(textforgui{:});
            h = text(0,0.45, allstr);
            if colorcount >= 15, set(h, 'fontsize', 8);end;
            if colorcount >= 20, set(h, 'fontsize', 6);end;
            if strcmp(BACKCOLOR, 'k'), set(h, 'color', 'w'); end;
        end;
        axis off;
        return;
    elseif strcmpi(g.summary, '3d')
        options = { 'gui', 'off', 'dipolesize', g.dipolesize/1.5,'dipolelength', g.dipolelength, 'sphere', g.sphere, 'spheres', g.spheres ...
                    'color', g.color, 'mesh', g.mesh, 'num', g.num, 'image', g.image 'normlen' g.normlen ...
                    'coordformat' g.coordformat 'mri' g.mri 'meshdata' g.meshdata 'axistight' g.axistight };
        figure('position', [ 100 600 600 200 ]); 
        axes('position', [-0.1 -0.1 1.2 1.2], 'color', 'k'); axis off; blackimg = zeros(10,10,3); image(blackimg);
        axes('position', [0   0 1/3 1], 'tag', 'rear'); dipplot(sourcesori, options{:}, 'holdon', 'on'); view([0 -1 0]);
        axes('position', [1/3 0 1/3 1], 'tag', 'top' ); dipplot(sourcesori, options{:}, 'holdon', 'on'); view([0  0 1]);
        axes('position', [2/3 0 1/3 1], 'tag', 'side'); dipplot(sourcesori, options{:}, 'holdon', 'on'); view([1 -0.01 0]);
        set(gcf, 'paperpositionmode', 'auto');
        return;
    end;
        
    % plot head graph in 3D
    % ---------------------
    if strcmp(g.gui, 'on')
        fig = figure('visible', g.plot); 
        pos = get(gca, 'position');
        set(gca, 'position', [pos(1)+0.05 pos(2:end)]);
    end;
    indx = ceil(dat.imgcoords{1}(end)/2);
    indy = ceil(dat.imgcoords{2}(end)/2);
    indz = ceil(dat.imgcoords{3}(end)/2);
    if strcmpi(g.holdon, 'off')
        plotimgs( dat, [indx indy indz], dat.transform);
    
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
        if strcmpi(g.coordformat, 'CTF'), g.view(2:3) = -g.view(2:3); end;
        view(g.view);
        %set(gca, 'cameratarget',   dat.zeroloc); % disable change size
        %set(gca, 'cameraposition', dat.zeroloc+g.view*g.zoom); % disable change size
        axis off;
    end;
        
    % plot sphere mesh and nose
    % -------------------------
    if strcmpi(g.holdon, 'off')
        if isempty(g.meshdata)
            SPHEREGRAIN = 20; % 20 is also Matlab default
            [x y z] = sphere(SPHEREGRAIN);
            hold on; 
            [xx yy zz] = transform(x*0.085, y*0.085, z*0.085, dat.sph2spm);
            [xx yy zz] = transform(x*85   , y*85   , z*85   , dat.sph2spm);
            %xx = x*100;
            %yy = y*100;
            %zz = z*100;
            if strcmpi(COLORMESH, 'w')
                hh = mesh(xx, yy, zz, 'cdata', ones(21,21,3), 'tag', 'mesh'); hidden off;
        else
            hh = mesh(xx, yy, zz, 'cdata', zeros(21,21,3), 'tag', 'mesh'); hidden off;
            end;
        else
            try, 
            if isstr(g.meshdata)
                tmp = load('-mat', g.meshdata);
                g.meshdata = { 'vertices' tmp.vol.bnd(1).pnt 'faces' tmp.vol.bnd(1).tri };
            end;
            hh = patch(g.meshdata{:}, 'facecolor', 'none', 'edgecolor', COLORMESH, 'tag', 'mesh');
            catch, disp('Unrecognize model file (probably CTF)'); end;
        end;
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
    
    % diph = gca; % DEBUG
    % colormap('jet');
    % cbar
    % axes(diph);

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
            if isfield(sources, 'active'),
                nbdip = length(sources(index).active);
            end;
        end;
        
        % dipole length
        % -------------
        multfactor = 1;
        if strcmpi(g.normlen, 'on')
            if nbdip == 1
                len    = sqrt(sum(sources(index).momxyz(1,:).^2));
            else
                len1   = sqrt(sum(sources(index).momxyz(1,:).^2));
                len2   = sqrt(sum(sources(index).momxyz(2,:).^2));
                len    = mean([len1 len2]);
            end;
            if strcmpi(g.coordformat, 'CTF'), len = len*10; end;
            if len ~= 0, multfactor = 15/len; end;
        else
            if strcmpi(g.coordformat, 'spherical')
                 multfactor = 100;
            else multfactor = 1.5;
            end;            
        end;
        
        for dip = 1:nbdip
        
            x = sources(index).posxyz(dip,1);
            y = sources(index).posxyz(dip,2);
            z = sources(index).posxyz(dip,3);
            
            xo = sources(index).momxyz(dip,1)*g.dipolelength*multfactor;
            yo = sources(index).momxyz(dip,2)*g.dipolelength*multfactor;
            zo = sources(index).momxyz(dip,3)*g.dipolelength*multfactor;
            
            xc = 0;
            yc = 0;
            zc = 0;
            
            centvec = [xo-xc yo-yc zo-zc]; % vector pointing into center
            dipole_orient = [x+xo y+yo z+zo]/norm([x+xo y+yo z+zo]);
            c = dot(centvec, dipole_orient);
            
            if strcmpi(g.pointout,'on')
                if (c < 0) | (abs([x+xo,y+yo,z+zo]) < abs([x,y,z]))
                    xo1 = x-xo; % make dipole point outward from head center
                    yo1 = y-yo;
                    zo1 = z-zo; 
                    %fprintf('invert because: %e  \n', c);
                else
                    xo1 = x+xo;
                    yo1 = y+yo;
                    zo1 = z+zo;
                    %fprintf('NO invert because: %e  \n', c);
                end
            else
                xo1 = x+xo;
                yo1 = y+yo;
                zo1 = z+zo;
                %fprintf('NO invert because: %e  \n', c);
            end
            
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% draw dipole bar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            tag = [ 'dipole' num2str(index) ];
            
            % from spherical to electrode space
            % ---------------------------------
            [xx   yy   zz]   = transform(x,   y,   z,   dat.sph2spm); % nothing happens for BEM
            [xxo1 yyo1 zzo1] = transform(xo1, yo1, zo1, dat.sph2spm); % because dat.sph2spm = []

            if ~strcmpi(g.spheres,'on') % plot dipole direction lines
               h1 = line( [xx xxo1]', [yy yyo1]', [zz zzo1]');
               
            elseif g.dipolelength>0 % plot dipole direction cylinders with end cap patch
            
               [xc yc zc] = cylinder( 2, 10);
               [xs ys zs] = sphere(10);
               xc = [ xc; -xs(7:11,:)*2  ];
               yc = [ yc; -ys(7:11,:)*2  ];
               zc = [ zc; zs(7:11,:)/5+1 ];
               
               colorarray = repmat(reshape(g.color{index}, 1,1,3), [size(zc,1) size(zc,2) 1]);
               handles = surf(xc, yc, zc, colorarray, 'tag', tag, 'edgecolor', 'none', ...
                              'backfacelighting', 'lit', 'facecolor', 'interp', 'facelighting', ...
                              'phong', 'ambientstrength', 0.3);
               [xc yc zc] = adjustcylinder2( handles, [xx yy zz], [xxo1 yyo1 zzo1] );

               cx = mean(xc,2); %cx = [(3*cx(1)+cx(2))/4; (cx(1)+3*cx(2))/4];
               cy = mean(yc,2); %cy = [(3*cy(1)+cy(2))/4; (cy(1)+3*cy(2))/4];
               cz = mean(zc,2); %cz = [(3*cz(1)+cz(2))/4; (cz(1)+3*cz(2))/4];
               tmpx = xc - repmat(cx, [1 size(xc, 2)]);
               tmpy = yc - repmat(cy, [1 size(xc, 2)]);
               tmpz = zc - repmat(cz, [1 size(xc, 2)]);
               l=sqrt(tmpx.^2+tmpy.^2+tmpz.^2);
               warning('off', 'MATLAB:divideByZero');                       % this is due to a Matlab 2008b (or later) 
               normals = reshape([tmpx./l tmpy./l tmpz./l],[size(tmpx) 3]); % in the rotate function in adjustcylinder2
               warning('off', 'MATLAB:divideByZero');                       % one of the z (the last row is not rotated)
               set( handles, 'vertexnormals', normals);
               
            end

            [xxmri   yymri   zzmri  ] = transform(xx,   yy,   zz,   pinv(dat.transform));
            [xxmrio1 yymrio1 zzmrio1] = transform(xxo1, yyo1, zzo1, pinv(dat.transform));
            dipstruct.mricoord  = [xxmri yymri zzmri];   % Coordinates in MRI space
            dipstruct.eleccoord = [ xx yy zz ];          % Coordinates in elec space
            dipstruct.posxyz    = sources(index).posxyz; % Coordinates in spherical space
            outsources(index).eleccoord(dip,:) = [xx yy zz];
            outsources(index).mnicoord(dip,:) = [xx    yy    zz];
            outsources(index).mricoord(dip,:) = [xxmri yymri zzmri];
            outsources(index).talcoord(dip,:) = mni2tal([xx yy zz]')';
            dipstruct.talcoord                = mni2tal([xx yy zz]')';
            
            % copy for output
            % ---------------
            XX(index) = xxmri;
            YY(index) = yymri;
            ZZ(index) = zzmri;
            XO(index) = xxmrio1;
            YO(index) = yymrio1;
            ZO(index) = zzmrio1;

            if isempty(g.dipnames)
                dipstruct.rv   = sprintf('%3.2f', sources(index).rv*100);
                dipstruct.name = sources(index).component;
            else
                dipstruct.rv   = sprintf('%3.2f', sources(index).rv*100);
                dipstruct.name = g.dipnames{index};
            end;
            if ~strcmpi(g.spheres,'on') % plot disk markers
               set(h1,'userdata',dipstruct,'tag',tag,'color','k','linewidth',g.dipolesize(index)/7.5);
               if strcmp(BACKCOLOR, 'k'), set(h1, 'color', g.color{index}); end;
            end
            
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%% draw sphere or disk marker %%%%%%%%%%%%%%%%%%%%%%%%%
            %
            hold on;
            if strcmpi(g.spheres,'on') % plot spheres
                if strcmpi(g.projimg, 'on')
                    if strcmpi(g.verbose, 'on'),
                        disp('Warning: projections cannot be plotted for 3-D sphere');
                    end
                    %tmpcolor = g.color{index} / 2;
                    %h = plotsphere([xx yy zz], g.dipolesize/6, 'color', g.color{index}, 'proj', ...
                    %               [dat.imgcoords{1}(1) dat.imgcoords{2}(end) dat.imgcoords{3}(1)]*97/100, 'projcol', tmpcolor);
  
                    %set(h(2:end), 'userdata', 'proj', 'tag', tag);
                else
                    %h = plotsphere([xx yy zz], g.dipolesize/6, 'color', g.color{index});
                end;                    
                h = plotsphere([xx yy zz], g.dipolesize(index)/6, 'color', g.color{index});
                set(h(1), 'userdata', dipstruct, 'tag', tag);
            else % plot dipole markers
               h = plot3(xx,  yy,  zz); 
               set(h, 'userdata', dipstruct, 'tag', tag, ...
                   'marker', '.', 'markersize', g.dipolesize(index), 'color', g.color{index});
            end
            
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% project onto images %%%%%%%%%%%%%%%%%%%%%%%%%
            %
            [tmp1xx   tmp1yy   tmp1zz  ] = transform( xxmri  , yymri  , dat.imgcoords{3}(1), dat.transform);
            [tmp1xxo1 tmp1yyo1 tmp1zzo1] = transform( xxmrio1, yymrio1, dat.imgcoords{3}(1), dat.transform);
            [tmp2xx   tmp2yy   tmp2zz  ] = transform( xxmri  , dat.imgcoords{2}(end), zzmri  , dat.transform);
            [tmp2xxo1 tmp2yyo1 tmp2zzo1] = transform( xxmrio1, dat.imgcoords{2}(end), zzmrio1, dat.transform);
            [tmp3xx   tmp3yy   tmp3zz  ] = transform( dat.imgcoords{1}(1), yymri  , zzmri  , dat.transform);
            [tmp3xxo1 tmp3yyo1 tmp3zzo1] = transform( dat.imgcoords{1}(1), yymrio1, zzmrio1, dat.transform);
            
            if strcmpi(g.projimg, 'on') & strcmpi(g.spheres, 'off')
                tmpcolor = g.projcol{index};
                
                % project onto z axis
                tag = [ 'dipole' num2str(index) ];
                if ~strcmpi(g.image, 'besa')
                    h = line( [tmp1xx tmp1xxo1]', [tmp1yy tmp1yyo1]', [tmp1zz tmp1zzo1]');
                    set(h, 'userdata', 'proj', 'tag', tag, 'color','k', 'linewidth', g.dipolesize(index)/7.5);
                end;
                if strcmp(BACKCOLOR, 'k'), set(h, 'color', tmpcolor); end;
                h = plot3(tmp1xx,  tmp1yy,  tmp1zz); 
                set(h, 'userdata', 'proj', 'tag', tag, ...
                       'marker', '.', 'markersize', g.dipolesize(index), 'color', tmpcolor);
                
                % project onto y axis
                tag = [ 'dipole' num2str(index) ];
                if ~strcmpi(g.image, 'besa')
                    h = line( [tmp2xx tmp2xxo1]', [tmp2yy tmp2yyo1]', [tmp2zz tmp2zzo1]');
                    set(h, 'userdata', 'proj', 'tag', tag, 'color','k', 'linewidth', g.dipolesize(index)/7.5);
                end;
                if strcmp(BACKCOLOR, 'k'), set(h, 'color', tmpcolor); end;
                h = plot3(tmp2xx,  tmp2yy,  tmp2zz); 
                set(h, 'userdata', 'proj', 'tag', tag, ...
                       'marker', '.', 'markersize', g.dipolesize(index), 'color', tmpcolor);
                   
                % project onto x axis
                tag = [ 'dipole' num2str(index) ];
                if ~strcmpi(g.image, 'besa')
                    h = line( [tmp3xx tmp3xxo1]', [tmp3yy tmp3yyo1]', [tmp3zz tmp3zzo1]');
                    set(h, 'userdata', 'proj', 'tag', tag, 'color','k', 'linewidth', g.dipolesize(index)/7.5);
                end;
                if strcmp(BACKCOLOR, 'k'), set(h, 'color', tmpcolor); end;
                h = plot3(tmp3xx,  tmp3yy,  tmp3zz); 
                set(h, 'userdata', 'proj', 'tag', tag, ...
                       'marker', '.', 'markersize', g.dipolesize(index), 'color', tmpcolor);
            end;

            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% project onto axes %%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if strcmpi(g.projlines, 'on')                
                clear h;
                % project onto z axis
                tag = [ 'dipole' num2str(index) ];
                h(1) = line( [xx tmp1xx]', [yy tmp1yy]', [zz tmp1zz]);
                set(h(1), 'userdata', 'proj', 'linestyle', '--', ...
                             'tag', tag, 'color', g.color{index}, 'linewidth', g.dipolesize(index)/7.5/5);
                
                % project onto x axis
                tag = [ 'dipole' num2str(index) ];
                h(2) = line( [xx tmp2xx]', [yy tmp2yy]', [zz tmp2zz]);
                set(h(2), 'userdata', 'proj', 'linestyle', '--', ...
                             'tag', tag, 'color', g.color{index}, 'linewidth', g.dipolesize(index)/7.5/5);
                
                % project onto y axis
                tag = [ 'dipole' num2str(index) ];
                h(3) = line( [xx tmp3xx]', [yy tmp3yy]', [zz tmp3zz]);
                set(h(3), 'userdata', 'proj', 'linestyle', '--', ...
                             'tag', tag, 'color', g.color{index}, 'linewidth', g.dipolesize(index)/7.5/5);
                if ~isempty(g.projcol)
                    set(h, 'color', g.projcol{index});
                end;
            end;
            %           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% draw text  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if isfield(sources, 'component')
                if strcmp(g.num, 'on')
                    h = text(xx,  yy,  zz, [ '  ' int2str(sources(index).component)]);
                    set(h, 'userdata', dipstruct, 'tag', tag, 'fontsize', g.dipolesize(index)/2 );
                    if ~strcmpi(g.image, 'besa'), set(h, 'color', 'w'); end;
                end;
            end;
        end;
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 3-D settings
    if strcmpi(g.spheres, 'on')
        lighting phong;
        material shiny;
        camlight left;
        camlight right;
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% draw elipse for group of dipoles  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % does not work because of new scheme, have to be reprogrammed
    
    %if ~isempty(g.std)
    %    for index = 1:length(g.std)
    %        if ~iscell(g.std{index})
    %            plotellipse(sources, g.std{index}, 1, dat.tcparams, dat.coreg);
    %        else
    %            sc = plotellipse(sources, g.std{index}{1}, g.std{index}{2}, dat.tcparams, dat.coreg);
    %            if length( g.std{index} ) > 2
    %                set(sc, g.std{index}{3:end});
    %            end;
    %         end;
    %     end;
    % end;

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
    cbplot = [ 'if strcmpi(get(gcbo, ''string''), ''plot one''),' ...
               '    for tmpi = 1:' nbsrc ',' ...
               '        set(findobj(''parent'', gca, ''tag'', [ ''dipole'' int2str(tmpi) ]), ''visible'', ''off'');' ...
               '    end; clear tmpi;' ...
               '    dipplot(gcbf);' ...               
               '    set(gcbo, ''string'', ''Plot all'');' ...
                   'else,' ...
               '    for tmpi = 1:' nbsrc ',' ...
               '        set(findobj(''parent'', gca, ''tag'', [ ''dipole'' int2str(tmpi) ]), ''visible'', ''on'');' ...
               '    end; clear tmpi;' ...
               '    set(gcbo, ''string'', ''Plot one'');' ...
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
    viewstring = fastif(dat.axistight, 'Loose view', 'Tight view');
    enmesh     = fastif(isempty(g.meshdata) & strcmpi(g.coordformat, 'MNI'), 'off', 'on');
    if strcmpi(g.coordformat, 'CTF'), viewcor = 'view([0 1 0]);';  viewtop = 'view([0 0 -1]);'; vis = 'off';
    else                              viewcor = 'view([0 -1 0]);'; viewtop = 'view([0 0 1]);';  vis = 'on';
    end;
    
    h = uicontrol( 'unit', 'normalized', 'position', [0 0  .15 1], 'tag', 'tmp', ...
                   'style', 'text', 'string',' ');
    h = uicontrol( 'unit', 'normalized', 'position', [0 0 .15 .05], 'tag', 'tmp', ...
                       'style', 'pushbutton', 'fontweight', 'bold', 'string', 'No controls', 'callback', ...
                   'set(findobj(''parent'', gcbf, ''tag'', ''tmp''), ''visible'', ''off'');');
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.05 .15 .05], 'tag', 'tmp', ...
                   'style', 'pushbutton', 'string', 'Top view', 'callback', viewtop);
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.1 .15 .05], 'tag', 'tmp', ...
                   'style', 'pushbutton', 'string', 'Coronal view', 'callback', viewcor);
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.15 .15 .05], 'tag', 'tmp', ...
                   'style', 'pushbutton', 'string', 'Sagittal view', 'callback', 'view([1 0 0]);');
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.2  .15 .05], 'tag', 'tmp', ...
                   'style', 'pushbutton', 'string', viewstring, 'callback', cbview);
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.25 .15 .05], 'tag', 'tmp', ...
                   'style', 'pushbutton', 'string', 'Mesh on', 'userdata', 0, 'callback', ...
                   cbmesh, 'enable', enmesh, 'visible', vis );
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.3 .15 .05], 'tag', 'tmp', ...
                   'style', 'text', 'string', 'Display:','fontweight', 'bold' );
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.35 .15 .02], 'tag', 'tmp',...
                   'style', 'text', 'string', '');
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.37 .15 .05], 'tag', 'tmp','userdata', 'z',...
                   'style', 'text', 'string', 'Z:', 'visible', vis );
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.42 .15 .05], 'tag', 'tmp','userdata', 'y', ...
                   'style', 'text', 'string', 'Y:', 'visible', vis );
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.47 .15 .05], 'tag', 'tmp', 'userdata', 'x',...
                   'style', 'text', 'string', 'X:', 'visible', vis );
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.52 .15 .05], 'tag', 'tmp', 'userdata', 'rv',...
                   'style', 'text', 'string', 'RV:' );
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.57 .15 .05], 'tag', 'tmp', 'userdata', 'comp', ...
                   'style', 'text', 'string', '');
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.62 .15 .05], 'tag', 'tmp', 'userdata', 'editor', ...
                   'style', 'edit', 'string', '1', 'callback', ...
                   [ 'dipplot(gcbf);' ] );
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.67 .15 .05], 'tag', 'tmp', ...
                       'style', 'pushbutton', 'string', 'Keep|Prev', 'callback', ...
                   [ 'editobj = findobj(''parent'', gcf, ''userdata'', ''editor'');' ...
                     'set(editobj, ''string'', num2str(str2num(get(editobj, ''string''))-1));' ...
                     'tmpobj = get(gcf, ''userdata'');' ...
                     'eval(get(editobj, ''callback''));' ...
                     'set(tmpobj, ''visible'', ''on'');' ...
                     'clear editobj tmpobj;' ]);
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.72 .15 .05], 'tag', 'tmp', ...
                   'style', 'pushbutton', 'string', 'Prev', 'callback',  ...
                   [ 'editobj = findobj(''parent'', gcf, ''userdata'', ''editor'');' ...
                     'set(editobj, ''string'', num2str(str2num(get(editobj, ''string''))-1));' ...
                     'eval(get(editobj, ''callback''));' ...
                     'clear editobj;' ]);
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.77 .15 .05], 'tag', 'tmp', ...
                   'style', 'pushbutton', 'string', 'Next', 'callback', ...
                   [ 'editobj = findobj(''parent'', gcf, ''userdata'', ''editor'');' ...
                     'set(editobj, ''string'', num2str(str2num(get(editobj, ''string''))+1));' ...
                     'dipplot(gcbf);' ...
                     'clear editobj;' ]);
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.82 .15 .05], 'tag', 'tmp', ...
                   'style', 'pushbutton', 'string', 'Keep|Next', 'callback', ...
                   [ 'editobj = findobj(''parent'', gcf, ''userdata'', ''editor'');' ...
                     'set(editobj, ''string'', num2str(str2num(get(editobj, ''string''))+1));' ...
                     'tmpobj = get(gcf, ''userdata'');' ...
                     'dipplot(gcbf);' ...
                     'set(tmpobj, ''visible'', ''on'');' ...
                     'clear editobj tmpobj;' ]);
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.87 .15 .05], 'tag', 'tmp', ...
                   'style', 'pushbutton', 'string', 'Plot one', 'callback', cbplot);
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.92 .15 .05], 'tag', 'tmp', ...
                   'style', 'text', 'string', [num2str(length(sources)) ' dipoles:'], 'fontweight', 'bold' );
    h = uicontrol( 'unit', 'normalized', 'position', [0 0.97 .15 .05], 'tag', 'tmp', ...
                   'style', 'text', 'string', '');
    set(gcf, 'userdata', findobj('parent', gca, 'tag', 'dipole1'));
    
    dat.nbsources  = length(sources);
    set(gca, 'userdata', dat ); % last param=1 for MRI view tight/loose
    set(gcf, 'color', BACKCOLOR);
    
    if strcmp(g.gui, 'off') | strcmpi(g.holdon, 'on')
        set(findobj('parent', gcf, 'tag', 'tmp'), 'visible', 'off');
    end;
    if strcmp(g.mesh, 'off')
        set(findobj('parent', gca, 'tag', 'mesh'), 'visible', 'off');
    end;
    
    updatedipplot(gcf);
    
    rotate3d on;
    
    % close figure if necessary
    if strcmpi(g.plot, 'off')
        try, close(fig); catch, end;
    end;
    
    if strcmpi(g.holdon, 'on')
        box off;
        axis equal;
        axis off;
    end;

    % set camera positon
    if strcmpi(g.camera, 'set')
        set(gca, 'CameraPosition', [2546.94 -894.981 689.613], ...
                 'CameraPositionMode', 'manual', ...
                 'CameraTarget', [0 -18 18], ...
                 'CameraTargetMode', 'manual', ...
                 'CameraUpVector', [0 0 1], ...
                 'CameraUpVectorMode', 'manual', ...
                 'CameraViewAngle', [3.8815], ...
                 'CameraViewAngleMode', 'manual');
    end;
    
return;

% electrode space to MRI space
% ============================
function [x,y,z] = transform(x, y, z, transmat);
    
    if isempty(transmat), return; end;
    for i = 1:size(x,1)
        for j = 1:size(x,2)
            tmparray = transmat * [ x(i,j) y(i,j) z(i,j) 1 ]';
            x(i,j) = tmparray(1);
            y(i,j) = tmparray(2);
            z(i,j) = tmparray(3);
        end;
    end;
    
    
% does not work any more
% ----------------------
function sc = plotellipse(sources, ind, nstd, TCPARAMS, coreg);

    for i = 1:length(ind)
        tmpval(1,i) = -sources(ind(i)).posxyz(1);    
        tmpval(2,i) = -sources(ind(i)).posxyz(2);    
        tmpval(3,i) = sources(ind(i)).posxyz(3);
        [tmpval(1,i) tmpval(2,i) tmpval(3,i)] = transform(tmpval(1,i), tmpval(2,i), tmpval(3,i), TCPARAMS);
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
            [x y z]    = sph2cart(theta/180*pi, phi/180*pi, postmp(3)/1.2); 
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
   
   % find all dipolar structures
   % ---------------------------
   index   = 1;
   count   = 1;
   for index = 1:length(newdip)
       if isstruct( get(newdip(index), 'userdata') )
           dip_mricoord(count,:) = getfield(get(newdip(index), 'userdata'), 'mricoord');
           count = count+1;
           foundind = index;
       end;
   end;
 
   % get residual variance
   % ---------------------
   if exist('foundind')
      tmp = get(newdip(foundind), 'userdata');
      tal = tmp.talcoord;
      if ~isstr( tmp.name )
           tmprvobj = findobj('parent', fig, 'userdata', 'comp'); set( tmprvobj(end), 'string', [ 'Comp: ' int2str(tmp.name) ] );
      else tmprvobj = findobj('parent', fig, 'userdata', 'comp'); set( tmprvobj(end), 'string', tmp.name );
      end;
      tmprvobj = findobj('parent', fig, 'userdata', 'rv');   set( tmprvobj(end), 'string', [ 'RV: ' tmp.rv '%' ] );
      tmprvobj = findobj('parent', fig, 'userdata', 'x');    set( tmprvobj(end), 'string', [ 'X tal: ' int2str(round(tal(1))) ]);
      tmprvobj = findobj('parent', fig, 'userdata', 'y');    set( tmprvobj(end), 'string', [ 'Y tal: ' int2str(round(tal(2))) ]);
      tmprvobj = findobj('parent', fig, 'userdata', 'z');    set( tmprvobj(end), 'string', [ 'Z tal: ' int2str(round(tal(3))) ]);
   end
   
   % adapt the MRI to the dipole depth
   % ---------------------------------
   delete(findobj('parent', gca, 'tag', 'img'));
   
   tmpdiv1 = dat.imgcoords{1}(2)-dat.imgcoords{1}(1);
   tmpdiv2 = dat.imgcoords{2}(2)-dat.imgcoords{2}(1);
   tmpdiv3 = dat.imgcoords{3}(2)-dat.imgcoords{3}(1);
   if ~dat.axistight
       [xx yy zz] = transform(0,0,0, pinv(dat.transform)); % elec -> MRI space
       indx = minpos(dat.imgcoords{1}-zz);
       indy = minpos(dat.imgcoords{2}-yy);
       indz = minpos(dat.imgcoords{3}-xx);
   else
       if ~dat.cornermri
           indx = minpos(dat.imgcoords{1} - mean(dip_mricoord(:,1))) - 3*tmpdiv1;
           indy = minpos(dat.imgcoords{2} - mean(dip_mricoord(:,2))) + 3*tmpdiv2;
           indz = minpos(dat.imgcoords{3} - mean(dip_mricoord(:,3))) - 3*tmpdiv3;
       else % no need to shift slice if not ploted close to the dipole
           indx = minpos(dat.imgcoords{1} - mean(dip_mricoord(:,1)));
           indy = minpos(dat.imgcoords{2} - mean(dip_mricoord(:,2)));
           indz = minpos(dat.imgcoords{3} - mean(dip_mricoord(:,3)));
       end;
   end;
   
   % middle of the brain
   % -------------------   
   plotimgs( dat, [indx indy indz], dat.transform);
   %end;
   	
% plot images (transmat is the uniform matrix MRI coords -> elec coords)
% ----------------------------------------------------------------------
function plotimgs(dat, mricoord, transmat);
   
    % loading images
    % --------------
    if ndims(dat.imgs) == 4 % true color data
        img1(:,:,3) = rot90(squeeze(dat.imgs(mricoord(1),:,:,3))); 
        img2(:,:,3) = rot90(squeeze(dat.imgs(:,mricoord(2),:,3))); 
        img3(:,:,3) = rot90(squeeze(dat.imgs(:,:,mricoord(3),3))); 
        img1(:,:,2) = rot90(squeeze(dat.imgs(mricoord(1),:,:,2))); 
        img2(:,:,2) = rot90(squeeze(dat.imgs(:,mricoord(2),:,2))); 
        img3(:,:,2) = rot90(squeeze(dat.imgs(:,:,mricoord(3),2))); 
        img1(:,:,1) = rot90(squeeze(dat.imgs(mricoord(1),:,:,1))); 
        img2(:,:,1) = rot90(squeeze(dat.imgs(:,mricoord(2),:,1))); 
        img3(:,:,1) = rot90(squeeze(dat.imgs(:,:,mricoord(3),1)));     
    else
        img1 = rot90(squeeze(dat.imgs(mricoord(1),:,:))); 
        img2 = rot90(squeeze(dat.imgs(:,mricoord(2),:))); 
        img3 = rot90(squeeze(dat.imgs(:,:,mricoord(3)))); 

        if ndims(img1) == 2, img1(:,:,3) = img1; img1(:,:,2) = img1(:,:,1); end;
        if ndims(img2) == 2, img2(:,:,3) = img2; img2(:,:,2) = img2(:,:,1); end;
        if ndims(img3) == 2, img3(:,:,3) = img3; img3(:,:,2) = img3(:,:,1); end;
    end;
    
    % computing coordinates for planes
    % --------------------------------    
    wy1 = [min(dat.imgcoords{2}) max(dat.imgcoords{2}); min(dat.imgcoords{2}) max(dat.imgcoords{2})];
    wz1 = [min(dat.imgcoords{3}) min(dat.imgcoords{3}); max(dat.imgcoords{3}) max(dat.imgcoords{3})];
    wx2 = [min(dat.imgcoords{1}) max(dat.imgcoords{1}); min(dat.imgcoords{1}) max(dat.imgcoords{1})];
    wz2 = [min(dat.imgcoords{3}) min(dat.imgcoords{3}); max(dat.imgcoords{3}) max(dat.imgcoords{3})];
    wx3 = [min(dat.imgcoords{1}) max(dat.imgcoords{1}); min(dat.imgcoords{1}) max(dat.imgcoords{1})];
    wy3 = [min(dat.imgcoords{2}) min(dat.imgcoords{2}); max(dat.imgcoords{2}) max(dat.imgcoords{2})];
    if dat.axistight & ~dat.cornermri
        wx1 = [ 1 1; 1 1]*dat.imgcoords{1}(mricoord(1));
        wy2 = [ 1 1; 1 1]*dat.imgcoords{2}(mricoord(2));
        wz3 = [ 1 1; 1 1]*dat.imgcoords{3}(mricoord(3));
    else
        wx1 =  [ 1 1; 1 1]*dat.imgcoords{1}(1);
        wy2 =  [ 1 1; 1 1]*dat.imgcoords{2}(end);
        wz3 =  [ 1 1; 1 1]*dat.imgcoords{3}(1);
    end;
    
    % transform MRI coordinates to electrode space
    % --------------------------------------------
    [ elecwx1 elecwy1 elecwz1 ] = transform( wx1, wy1, wz1, transmat);
    [ elecwx2 elecwy2 elecwz2 ] = transform( wx2, wy2, wz2, transmat);
    [ elecwx3 elecwy3 elecwz3 ] = transform( wx3, wy3, wz3, transmat);
    
    % ploting surfaces
    % ----------------
    options = { 'FaceColor','texturemap', 'EdgeColor','none', 'CDataMapping', ...
                'direct','tag','img', 'facelighting', 'none' };
    hold on;
    surface(elecwx1, elecwy1, elecwz1, img1(end:-1:1,:,:), options{:}); 
    surface(elecwx2, elecwy2, elecwz2, img2(end:-1:1,:,:), options{:});
    surface(elecwx3, elecwy3, elecwz3, img3(end:-1:1,:,:), options{:}); 
    %xlabel('x'); ylabel('y'); zlabel('z'); axis equal; dsaffd
    
    if strcmpi(dat.drawedges, 'on')
        % removing old edges if any
        delete(findobj( gcf, 'tag', 'edges'));
        if dat.axistight & ~dat.cornermri, col = 'k'; else col = [0.5 0.5 0.5]; end;
        h(1) = line([elecwx3(1) elecwx3(2)]', [elecwy3(1) elecwy2(1)]', [elecwz1(1) elecwz1(2)]'); % sagittal-transverse
        h(2) = line([elecwx3(1) elecwx2(3)]', [elecwy2(1) elecwy2(2)]', [elecwz1(1) elecwz1(2)]'); % coronal-tranverse
        h(3) = line([elecwx3(1) elecwx3(2)]', [elecwy2(1) elecwy2(2)]', [elecwz3(1) elecwz1(1)]'); % sagittal-coronal
        set(h, 'color', col, 'linewidth', 2, 'tag', 'edges');
    end;
        
    %%fill3([-2 -2 2 2], [-2 2 2 -2], wz(:)-1, BACKCOLOR);
    %%fill3([-2 -2 2 2], wy(:)-1, [-2 2 2 -2], BACKCOLOR);
    rotate3d on 

function index = minpos(vals);
	vals(find(vals < 0)) = inf;
	[tmp index] = min(vals);

function scalegca(multfactor)
    xl = xlim; xf = ( xl(2) - xl(1) ) * multfactor;
    yl = ylim; yf = ( yl(2) - yl(1) ) * multfactor;
    zl = zlim; zf = ( zl(2) - zl(1) ) * multfactor;
    xlim( [ xl(1)-xf xl(2)+xf ]);
    ylim( [ yl(1)-yf yl(2)+yf ]);
    zlim( [ zl(1)-zf zl(2)+zf ]);
    
function color = strcol2real(colorin, colmap)
    if ~iscell(colorin)
        for index = 1:length(colorin)
            color{index} = colmap(colorin(index),:);
        end;
    else
        color = colorin;
        for index = 1:length(colorin)
            if isstr(colorin{index})
                switch colorin{index}
                 case 'r', color{index} = [1 0 0];
                 case 'g', color{index} = [0 1 0];
                 case 'b', color{index} = [0 0 1];
                 case 'c', color{index} = [0 1 1];
                 case 'm', color{index} = [1 0 1];
                 case 'y', color{index} = [1 1 0];
                 case 'k', color{index} = [0 0 0];
                 case 'w', color{index} = [1 1 1];
                 otherwise, error('Unknown color');
                end;
            end;
        end;
    end;

function x = gammacorrection(x, gammaval);
    x = 255 * (double(x)/255).^ gammaval;
    % image is supposed to be scaled from 0 to 255
    % gammaval = 1 is identity of course
    
