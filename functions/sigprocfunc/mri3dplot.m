% mri3dplot() - plot 3-D density image translucently on top of the mean MR 
%               brain image used in dipplot(). Plot brain slices in directions
%               'top' (axial), or 'side' (sagittal), or 'rear' (coronal).
%               Creates a new figure(). Smoothing uses Matlab smooth3()
% Usage:
%      >> [smoothed_3ddens, mriplanes] = mri3dplot(array3d, mri, 'key', 'val', ...);
%
% Input: 
%   array3d     - 3-D array to plot translucently on top of MRI image planes
%                  (e.g., as returned by dipoledensity(), unit: dipoles/cc).
%   mri         - [string or struct] base MR image structure (as returned by 
%                 dipoledensity.m or mri file (matlab format or file format read 
%                 by fcdc_read_mri. See dipplot.m help for more information.
%
% Optional inputs:
%   'mriview'   - ['top'|'side'|rear'] MR image slices to plot. 'axial',
%                    'coronal', and 'saggital' are also recognized keywords
%                    {default|[]: 'top'|'axial'}. May also be a cell array
%                    of such keyword (one per slice).
%   'mrislices' - [real vector] MR slice plotting coordinates in 'mriview' 
%                    direction {default|[]: -50 to +50 in steps of 11 mm}
%   'kernel'    - 3-D smoothing ht, width, & depth in voxels {default|0: none}
%                    > 0 -> uses gaussian kernel of std. dev. 0.65 voxels)
%   'geom'      - [rows, cols] geometry of slice array in figure output
%                    {default: smallest enclosing near-square array}
%   'rotate'    - [0|90|180|270 deg] 2-D slice image rotation {default: 90}
%   'cmap'      - [float array] colormap for plotting the 3-D array 
%                    {default: 'hot'}
%   'cmax'      - [float] color palette max value {default: array3d max}
%   'cmin'      - [float] color palette min value {default: 0}
%   'cbar'      - ['on'|'off'] plot colorbar. Default is 'on'.  
%   'subplot'   - ['on'|'off'] for single slice only, plot within a sub-plot
%                 panel. If 'on', this automatically sets 'cbar' to 'off'.
%                 Default is 'off'.  
%   'mixfact'   - [float] factor for mixing the background image with the
%                 array3d information. Default is 0.5.
%   'mixmode'   - ['add'|'overwrite'] 'add' will allow for trasnparency
%                 while 'overwrite' will preserve the orginal MRI image 
%                 and overwrite the pixel showind density changes.
%
% Outputs:
%   smoothed_3ddens - the plotted (optionally smoothed) 3-D density matrix
%   mriplanes       - cell array of the plotted MR image slices
%
% Author: Arnaud Delorme & Scott Makeig, SCCN, 10 June 2003
%
% See also: plotmri()

% Copyright (C) Arnaud Delorme, sccn, INC, UCSD, 2003-
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
% Revision 1.5  2009/05/08 02:13:21  arno
% better aspect ratio for figure
%
% Revision 1.4  2009/05/08 02:04:43  arno
% better figure geometry
%
% Revision 1.3  2009/05/08 02:02:52  arno
% cbar etc...
%
% Revision 1.2  2009/05/07 23:35:36  arno
% custom view for each slice
%
% Revision 1.1  2009/05/07 23:32:47  arno
% adding mri3dplot
%
% Revision 1.18  2009/03/30 23:14:31  arno
% Implement the overwrite method
%
% Revision 1.16  2009/03/30 22:44:16  arno
% Same
%
% Revision 1.15  2009/03/30 22:43:42  arno
% mixmode implementation
%
% Revision 1.14  2009/03/28 00:22:02  nick
% allow for weighting of background image with foreground image
%
% Revision 1.13  2008/06/03 00:22:49  nima
% NaNs are now accepted and assumed to be points outside the brain.
%
% Revision 1.12  2008/05/23 22:07:01  nima
% cmin option added.
%
% Revision 1.11  2007/01/02 17:49:49  scott
% removed test for min color value - allow 0
%
% Revision 1.10  2007/01/02 16:51:39  scott
% added documentation to code!
%
% Revision 1.9  2006/12/22 05:53:37  scott
% added smoothing details to help msg
%
% Revision 1.8  2006/12/22 05:28:06  arno
% scott header edit
%
% Revision 1.7  2005/09/30 16:31:19  arno
% programming cmax; speed up by 10
%
% Revision 1.6  2005/09/30 15:41:12  arno
% fix background
%
% Revision 1.5  2005/09/28 22:07:42  arno
% adjusting colors
%
% Revision 1.4  2005/09/28 21:32:30  arno
% change default smoothing, remake axes
%
% Revision 1.3  2005/09/28 16:43:32  arno
% header
%
% Revision 1.2  2005/07/08 22:51:30  arno
% new version (completely remodeveled)
%

function [smoothprob3d, mriplanes] = mri3dplot(prob3d, mri, varargin)

    % REDUCEPATCH  Reduce number of patch faces.

    if nargin < 1
        help mri3dplot;
        return;
    end;

    DEFAULT_SPACING = 11;     % default MR slice interval (in mm)
    translucency    = 0.5;    % default alpha value
    mri_lim         = 85;     % +/- axis limits of MNI head image
    
    g = finputcheck( varargin, { 'mriview'   { 'string' 'cell' }  { { 'sagital' 'axial' 'coronal' ...
                                                          'top' 'side' 'rear' } {} }   'top';
                        'mixmode'   'string'   { 'add' 'overwrite' }     'add';
                        'mrislices' 'float'    []                        [];
                        'view'      'float'    []                        [];
                        'geom'      'float'    []                        [];
                        'cmap'      'float'    []                        jet;
                        'cmax'      'float'    []                        [];
                        'mixfact'   'float'    []                        0.5;
                        'cmin'      'float'    []                        0;
                        'cbar'      'string'   { 'on' 'off' }            'on';
                        'subplot'   'string'   { 'on' 'off' }            'off';
                        'rotate'    'integer'  { 0 90 180 270 }          90;
                        'kernel'    'float'    []                        0 });
    if isstr(g), error(g); end;
    
    if strcmpi(g.mriview,'sagittal'),    g.mriview = 'side'; 
    elseif strcmpi(g.mriview,'axial'),   g.mriview = 'top'; 
    elseif strcmpi(g.mriview,'coronal'), g.mriview = 'rear';
    end;
    if strcmpi(g.subplot, 'on') % plot colorbar
        g.cbar = 'off';
    end;
    if isstr(mri)
         try, 
            mri = load('-mat', mri);
            mri = mri.mri;
        catch,
            disp('Failed to read Matlab file. Attempt to read MRI file using function read_fcdc_mri');
            try,
                warning off;
                mri = read_fcdc_mri(mri);
                mri.anatomy = round(gammacorrection( mri.anatomy, 0.8));
                mri.anatomy = uint8(round(mri.anatomy/max(reshape(mri.anatomy, prod(mri.dim),1))*255));
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
    
    % normalize prob3d for 1 to ncolors and create 3-D dim
    % ----------------------------------------------------
    if g.kernel ~= 0
        disp('Smoothing...');
        smoothprob3d    = smooth3(prob3d, 'gaussian', g.kernel);
        prob3d          = smoothprob3d;
    end;

    ncolors = size(g.cmap,1);
    maxdens = max(prob3d(:));

    if isempty(g.cmax), g.cmax = maxdens; end;
        
    fprintf('Brightest color denotes a density of: %1.6f (presumed unit: dipoles/cc)\n', maxdens);

    prob3d    = round((prob3d-g.cmin)/(g.cmax - g.cmin)*(ncolors-1))+1; % project desnity image into the color space: [1:ncolors]
    prob3d( find(prob3d > ncolors) ) = ncolors;
    newprob3d = zeros(size(prob3d,1), size(prob3d,2), size(prob3d,3), 3);

    outOfBrainMask = find(isnan(prob3d)); % place NaNs in a mask, NaNs are assumed for points outside the brain
    prob3d(outOfBrainMask) = 1;
    
    tmp = g.cmap(prob3d,1); newprob3d(:,:,:,1) = reshape(tmp, size(prob3d));
    tmp = g.cmap(prob3d,2); newprob3d(:,:,:,2) = reshape(tmp, size(prob3d));
    tmp = g.cmap(prob3d,3); newprob3d(:,:,:,3) = reshape(tmp, size(prob3d));
    
    % plot MRI slices
    % ---------------
    if isempty(g.mrislices), 
        g.mrislices = linspace(-50, 50, DEFAULT_SPACING); 
    end;
    
    if strcmpi(g.cbar, 'on'), add1 = 1; else add1 = 0; end;
    if isempty(g.geom), 
        g.geom = ceil(sqrt(length(g.mrislices)+add1)); 
        g.geom(2) = ceil((length(g.mrislices)+add1)/g.geom);
    end;

    if strcmpi(g.subplot, 'off')
        fig = figure;
        
        pos = get(fig, 'position');
        set(fig, 'position', [ pos(1)+15 pos(2)+15 pos(3)/4*g.geom(1) pos(4)/3*g.geom(2) ]);
    end;
    
    disp('Plotting...');
    
    % custom view for each slice
    if ~iscell( g.mriview )
        g.mriview = { g.mriview };
        g.mriview(2:length( g.mrislices )) = g.mriview(1);
    end;
    
    for index = 1:length( g.mrislices ) %%%%%%% for each plotted MR image slice %%%%%%%%

        mysubplot(g.geom(1), g.geom(2), index); % get an image slice axis
        switch g.mriview{index}
         case 'side', coord = [  g.mrislices(index) 0 0 1 ]; 
         case 'top' , coord = [  0 0 g.mrislices(index) 1 ]; 
         case 'rear', coord = [  0 g.mrislices(index) 0 1 ]; 
        end;
        
        coord = round( pinv(mri.transform)*coord' )';
        
        switch g.mriview{index}
         case 'side', mriplot  = squeeze( mri.anatomy(coord(1), :, :) );
         case 'top' , mriplot  = squeeze( mri.anatomy(:, :, coord(3)) );
         case 'rear', mriplot  = squeeze( mri.anatomy(:, coord(2), :) );
        end;

        mriplot(:,:,2) = mriplot(:,:,1);
        mriplot(:,:,3) = mriplot(:,:,1);
        mriplot = rotatemat( mriplot, g.rotate );
        
        switch g.mriview{index}
         case 'side', densplot = squeeze( newprob3d  (coord(1), :, :, :) );
         case 'top' , densplot = squeeze( newprob3d  (:, :, coord(3), :) );
         case 'rear', densplot = squeeze( newprob3d  (:, coord(2), :, :) );
        end;

        densplot = rotatemat( densplot, g.rotate );

        if strcmpi(g.mixmode, 'add')
            densplot(isnan(densplot)) = 0; % do not plot colors outside the brain, as indicated by NaNs
            mriplot  = mriplot*g.mixfact + densplot*(1-g.mixfact); % Mix 1/2 MR image + 1/2 density image
        else
            indsnon0 = sum(densplot(:,:,:),3) > 0;
            tmpmri = mriplot(:,:,1); tmpdens = densplot(:,:,1); tmpmri(indsnon0) = tmpdens(indsnon0); mriplot(:,:,1) = tmpmri;
            tmpmri = mriplot(:,:,2); tmpdens = densplot(:,:,2); tmpmri(indsnon0) = tmpdens(indsnon0); mriplot(:,:,2) = tmpmri;
            tmpmri = mriplot(:,:,3); tmpdens = densplot(:,:,3); tmpmri(indsnon0) = tmpdens(indsnon0); mriplot(:,:,3) = tmpmri;
        end;
        
        mriplanes{index} = mriplot;
        
        imagesc(mriplot); % plot [background MR image + density image]

        axis off; hold on;

        xl = xlim;
        yl = ylim;
        zl = zlim;
        
        % options = { 'FaceColor','texturemap', 'EdgeColor','none', 'CDataMapping', ...
        %             'scaled','facelighting', 'none', 'facealpha', translucency};
        % h = surface( [xl(1) xl(2); xl(1) xl(2)], [yl(1) yl(1); yl(2) yl(2)], ...
        %              [1 1; 1 1], densplot, options{:});

        axis equal;
        if ~isempty(g.view), view(g.view); end;
        title( [ int2str(g.mrislices(index)) ' mm' ], 'color', 'w');
    end;
    
    % plot colorbar
    % -------------
    if strcmpi(g.cbar, 'on') % plot colorbar

        h = mysubplot(g.geom(1), g.geom(2), length(g.mrislices)+1);
        pos = get(h, 'position');
        pos(1) = pos(1)+pos(3)/3;
        pos(3) = pos(3)/6;
        pos(2) = pos(2)+pos(4)/5;
        pos(4) = pos(4)*3/5;
        axis off;
        h = axes('unit', 'normalized', 'position', pos); % position the colorbar

        if strcmpi(g.mixmode, 'add')
            tmpmap = g.cmap/2 + ones(size(g.cmap))/4;  % restrict range to [1/4, 3/4] of cmap
        else
            tmpmap = g.cmap; %g.cmap/2 + ones(size(g.cmap))/4;  % restrict range to [1/4, 3/4] of cmap
        end;
        colormap(tmpmap);
        cbar(h, [1:length(g.cmap)], [g.cmin g.cmax]);
        box off;
        set(h, 'ycolor', [0.7 0.7 0.7]);
    end;
    
    fprintf('\n');
    if exist('fig') == 1
        set(fig,'color', g.cmap(1,:)/2);
    end;
return;

function mat = rotatemat(mat, angle);
    
    if angle == 0, return; end;
    if angle >= 90,
        newmat(:,:,1) = rot90(mat(:,:,1));
        newmat(:,:,2) = rot90(mat(:,:,2));
        newmat(:,:,3) = rot90(mat(:,:,3));
        mat = newmat;
    end;
    if angle >= 180,
        newmat(:,:,1) = rot90(mat(:,:,1));
        newmat(:,:,2) = rot90(mat(:,:,2));
        newmat(:,:,3) = rot90(mat(:,:,3));
        mat = newmat;
    end;
    if angle >= 270,
        newmat(:,:,1) = rot90(mat(:,:,1));
        newmat(:,:,2) = rot90(mat(:,:,2));
        newmat(:,:,3) = rot90(mat(:,:,3));
        mat = newmat;
    end;

function h = mysubplot(geom1, geom2, coord);
    
    coord = coord-1;
    horiz_border = 0;
    vert_border  = 0.1;
    
    coordy = floor(coord/geom1);
    coordx = coord - coordy*geom1;
    
    posx   = coordx/geom1+horiz_border*1/geom1/2;
    posy   = 1-(coordy/geom2+vert_border*1/geom2/2)-1/geom2;
    width  = 1/geom1*(1-horiz_border);
    height = 1/geom2*(1- vert_border);
    
    h = axes('unit', 'normalized', 'position', [ posx posy width height ]);
    %h = axes('unit', 'normalized', 'position', [ coordx/geom1 1-coordy/geom2-1/geom2 1/geom1 1/geom2 ]);
    
