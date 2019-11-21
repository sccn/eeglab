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
%   'plotintersect' - ['on'|'off'] plot intersection between plotted slices.
%                 Default is 'on'.
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
% Example:
%   dipfitdefs;
%   load('-mat', template_models(1).mrifile); % load mri variable
%   array = gauss3d(91,109,91);
%   mri3dplot(array, mri);
%
% See also: plotmri()

% Copyright (C) Arnaud Delorme, sccn, INC, UCSD, 2003-
% 03/29/2013 Makoto. Line 370 added to avoid negative matrix indices.
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

function [newprob3dori, mriplanes] = mri3dplot(prob3d, mri, varargin)

    % REDUCEPATCH  Reduce number of patch faces.

    if nargin < 1
        help mri3dplot;
        return;
    end

    DEFAULT_SPACING = 11;     % default MR slice interval (in mm)
    translucency    = 0.5;    % default alpha value
    mri_lim         = 85;     % +/- axis limits of MNI head image
    
    g = finputcheck( varargin, { 'mriview'   { 'string','cell' 'real' }  { { 'sagital','axial','coronal','top','side','rear' } {} [] }   'top';
                        'mixmode'   'string'   { 'add','overwrite','min' }     'add';
                        'mrislices' 'float'    []                        [];
                        'view'      'float'    []                        [];
                        'geom'      'float'    []                        [];
                        'cmap'      'float'    []                        jet(64);
                        'cmax'      'float'    []                        [];
                        'mixfact'   'float'    []                        0.5;
                        'cmin'      'float'    []                        0;
                        'plotintersect' 'string'   { 'on','off' }            'on';
                        'cbar'      'string'   { 'on','off' }            'on';
                        'subplot'   'string'   { 'on','off' }            'off';
                        'rotate'    'integer'  { 0,90,180,270 }          90;
                        'kernel'    'float'    []                        0; 
                        'addrow'    'integer'  []                        0;
                        'fighandle' ''         []                        []});
    if ischar(g), error(g); end
    if ischar(g.mriview) == 1, g.plotintersect = 'off'; end
    if strcmpi(g.mriview,'sagittal'),    g.mriview = 'side'; 
    elseif strcmpi(g.mriview,'axial'),   g.mriview = 'top'; 
    elseif strcmpi(g.mriview,'coronal'), g.mriview = 'rear';
    end
    if strcmpi(g.subplot, 'on') % plot colorbar
        g.cbar = 'off';
    end
    if ischar(mri)
         try
            mri = load('-mat', mri);
            mri = mri.mri;
         catch
            disp('Failed to read Matlab file. Attempt to read MRI file using function read_fcdc_mri');
            try
                warning off;
                mri = read_fcdc_mri(mri);
                mri.anatomy = round(gammacorrection( mri.anatomy, 0.8));
                mri.anatomy = uint8(round(mri.anatomy/max(reshape(mri.anatomy, prod(mri.dim),1))*255));
                % WARNING: if using double instead of int8, the scaling is different 
                % [-128 to 128 and 0 is not good]
                % WARNING: the transform matrix is not 1, 1, 1 on the diagonal, some slices may be 
                % misplaced
                warning on;
            catch
                error('Cannot load file using read_fcdc_mri');
            end
         end
    end
    
    % normalize prob3d for 1 to ncolors and create 3-D dim
    % ----------------------------------------------------
    if ~iscell(prob3d), prob3d = { prob3d }; end
    if length(prob3d) > 1
        if isempty(g.cmax), g.cmax = max(max(prob3d{1}(:)),max(prob3d{2}(:))); end
        [newprob3d{1}] = prepare_dens(prob3d{1}, g, 'abscolor');
        [newprob3d{2}] = prepare_dens(prob3d{2}, g, 'abscolor');
    else
        if isempty(g.cmax), g.cmax = max(prob3d{1}(:)); end
        [newprob3d{1}, maxdens1] = prepare_dens(prob3d{1}, g, 'usecmap');
    end
    fprintf('Brightest color denotes a density of: %1.6f (presumed unit: dipoles/cc)\n', g.cmax);
    
    % plot MRI slices
    % ---------------
    if isempty(g.mrislices), 
        g.mrislices = linspace(-50, 50, DEFAULT_SPACING); 
    end
    
    if strcmpi(g.cbar, 'on'), add1 = 1; else add1 = 0; end
    if isempty(g.geom)
        g.geom = ceil(sqrt(length(g.mrislices)+add1)); 
        g.geom(2) = ceil((length(g.mrislices)+add1)/g.geom)+g.addrow;
    end

    if strcmpi(g.subplot, 'off')
        if isempty(g.fighandle)
             fig = figure;
        else
             fig = g.fighandle;
             clf(fig);
        end
        
        pos = get(fig, 'position');
        set(fig, 'position', [ pos(1)+15 pos(2)+15 pos(3)/4*g.geom(1) pos(4)/3*g.geom(2) ]);
    end
    
    disp('Plotting...');
    
    % custom view for each slice
    if ~iscell( g.mriview )
        g.mriview = { g.mriview };
        g.mriview(2:length( g.mrislices )) = g.mriview(1);
    end
    
    newprob3dori = newprob3d;
    for index = 1:length( g.mrislices ) %%%%%%% for each plotted MR image slice %%%%%%%%

        % plot intersection between plotted slices
        % ----------------------------------------
        newprob3d = newprob3dori;
        if strcmpi(g.plotintersect, 'on')
            for index2 = setdiff_bc(1:length( g.mrislices ), index)
                switch g.mriview{index2}
                    case 'side', coord = [  g.mrislices(index2) 0 0 1 ]; 
                    case 'top' , coord = [  0 0 g.mrislices(index2) 1 ]; 
                    case 'rear', coord = [  0 g.mrislices(index2) 0 1 ]; 
                end
                coord = round( pinv(mri.transform)*coord' )';
                for i = 1:length(newprob3d)
                    switch g.mriview{index2}
                     case 'side', newprob3d{i}(  coord(1), :, :, :) = 0;
                     case 'top' , newprob3d{i}(  :, :, coord(3), :) = 0;
                     case 'rear', newprob3d{i}(  :, coord(2), :, :) = 0;
                    end
                end
            end
        end

        % create axis if necessary
        % ------------------------
        if strcmpi(g.subplot, 'off')
            mysubplot(g.geom(1), g.geom(2), index); % get an image slice axis
        end
        
        % find coordinate
        % ---------------
        switch g.mriview{index}
         case 'side', coord = [  g.mrislices(index) 0 0 1 ]; 
         case 'top' , coord = [  0 0 g.mrislices(index) 1 ]; 
         case 'rear', coord = [  0 g.mrislices(index) 0 1 ]; 
        end
        
        coord = round( pinv(mri.transform)*coord' )';
        
        % get MRI slice
        % -------------
        switch g.mriview{index}
         case 'side', mriplot  = squeeze( mri.anatomy(coord(1), :, :) );
         case 'top' , mriplot  = squeeze( mri.anatomy(:, :, coord(3)) );
         case 'rear', mriplot  = squeeze( mri.anatomy(:, coord(2), :) );
        end

        mriplot(:,:,2) = mriplot(:,:,1);
        mriplot(:,:,3) = mriplot(:,:,1);
        mriplot = rotatemat( mriplot, g.rotate );
        
        % get dipole density slice
        % ------------------------
        for i = 1:length(newprob3d)
            switch g.mriview{index}
             case 'side', densplot{i} = squeeze( newprob3d{i}(coord(1), :, :, :) );
             case 'top' , densplot{i} = squeeze( newprob3d{i}(:, :, coord(3), :) );
             case 'rear', densplot{i} = squeeze( newprob3d{i}(:, coord(2), :, :) );
            end
            densplot{i} = rotatemat( densplot{i}, g.rotate );
        end
  
        if isa(mriplot, 'uint8')
            % check if densplot is in uint8
            % if not - transform to uint8
            for dlen = 1:length(densplot)
                if ~isa(densplot{dlen}, 'uint8')
                    % check if multiply by 255 (when double
                    % with values ranging from 0 - 1) and
                    % then transform to uint8
                    checkval = densplot{dlen} >= 0 & densplot{dlen} <= 1;
                    checkval = sum(checkval(:)) == numel(densplot{dlen});
                    if checkval
                        densplot{dlen} = uint8(densplot{dlen} * 255); %#ok<AGROW>
                        continue
                    end
                    
                    % check if it's ok to transform
                    % straight to int8
                    checkval = densplot{dlen} >= 0 & densplot{dlen} <= 255;
                    checkval = sum(checkval(:)) == numel(densplot);
                    testint = isequal(densplot{dlen}, round(densplot{dlen}));
                    
                    if checkval && testint
                        densplot{dlen} = uint8(densplot{dlen}); %#ok<AGROW>
                    end
                    
                end
            end
            clear dlen checkval testint
        end
        
        if length(densplot) == 1
            densplot = densplot{1};
            if strcmpi(g.mixmode, 'add')
                densplot(isnan(densplot)) = 0; % do not plot colors outside the brain, as indicated by NaNs
                mriplot  = mriplot*g.mixfact + densplot*(1-g.mixfact); % Mix 1/2 MR image + 1/2 density image
            elseif strcmpi(g.mixmode, 'overwrite')
                indsnon0 = sum(densplot(:,:,:),3) > 0;
                tmpmri = mriplot(:,:,1); tmpdens = densplot(:,:,1); tmpmri(indsnon0) = tmpdens(indsnon0); mriplot(:,:,1) = tmpmri;
                tmpmri = mriplot(:,:,2); tmpdens = densplot(:,:,2); tmpmri(indsnon0) = tmpdens(indsnon0); mriplot(:,:,2) = tmpmri;
                tmpmri = mriplot(:,:,3); tmpdens = densplot(:,:,3); tmpmri(indsnon0) = tmpdens(indsnon0); mriplot(:,:,3) = tmpmri;
            elseif strcmpi(g.mixmode, 'min')
                densplot(isnan(densplot)) = 0; % do not plot colors outside the brain, as indicated by NaNs
                mriplot  = min(mriplot, densplot); % min
            end
            clear densplot;
        else
            densplot{1}(isnan(densplot{1})) = 0; % do not plot colors outside the brain, as indicated by NaNs
            densplot{2}(isnan(densplot{2})) = 0; % do not plot colors outside the brain, as indicated by NaNs
            if strcmpi(g.mixmode, 'add')
                mriplot(:,:,1) = mriplot(:,:,1)*g.mixfact + densplot{1}(:,:)*(1-g.mixfact); % min
                mriplot(:,:,3) = mriplot(:,:,3)*g.mixfact + densplot{2}(:,:)*(1-g.mixfact); % min
                mriplot(:,:,2) = mriplot(:,:,2)*g.mixfact; % min
            elseif strcmpi(g.mixmode, 'overwrite')
                indsnon01 = densplot{1}(:,:,:) > 0;
                indsnon02 = densplot{2}(:,:,:) > 0;
                tmpmri = mriplot(:,:,1); tmpdens = densplot{1}(:,:); tmpmri(indsnon01) = tmpdens(indsnon01); mriplot(:,:,1) = tmpmri;
                tmpmri = mriplot(:,:,2); tmpdens = densplot{2}(:,:); tmpmri(indsnon02) = tmpdens(indsnon02); mriplot(:,:,2) = tmpmri;
            else
                mriplot(:,:,1) = max(mriplot(:,:,1), densplot{1}); % min
                mriplot(:,:,2) = max(mriplot(:,:,2), densplot{2}); % min
            end
        end
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
        if ~isempty(g.view), view(g.view); end
        tit = title( [ int2str(g.mrislices(index)) ' mm' ]);
        if strcmpi(g.subplot, 'off'), set(tit, 'color', 'w'); end
    end
    
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
        end
        colormap(tmpmap);
        cbar(h, [1:length(g.cmap)], [g.cmin g.cmax]);
        box off;
        set(h, 'ycolor', [0.7 0.7 0.7]);
    end
    
    fprintf('\n');
    if exist('fig') == 1 
        if length(prob3d) == 1
            set(fig,'color', g.cmap(1,:)/2);
        else
            set(fig,'color', [0.0471 0.0471 0.0471]/1.3);
        end
    end
return;

function mat = rotatemat(mat, angle);

    if angle == 0, return; end
    if ndims(mat) == 2
        if angle >= 90,  mat = rot90(mat); end
        if angle >= 180, mat = rot90(mat); end
        if angle >= 270, mat = rot90(mat); end
    else
        if angle >= 90,
            newmat(:,:,1) = rot90(mat(:,:,1));
            newmat(:,:,2) = rot90(mat(:,:,2));
            newmat(:,:,3) = rot90(mat(:,:,3));
            mat = newmat;
        end
        if angle >= 180,
            newmat(:,:,1) = rot90(mat(:,:,1));
            newmat(:,:,2) = rot90(mat(:,:,2));
            newmat(:,:,3) = rot90(mat(:,:,3));
            mat = newmat;
        end
        if angle >= 270,
            newmat(:,:,1) = rot90(mat(:,:,1));
            newmat(:,:,2) = rot90(mat(:,:,2));
            newmat(:,:,3) = rot90(mat(:,:,3));
            mat = newmat;
        end
    end
    
function [newprob3d maxdens] = prepare_dens(prob3d, g, col);

    if g.kernel ~= 0
        disp('Smoothing...');
        smoothprob3d    = smooth3(prob3d, 'gaussian', g.kernel);
        prob3d          = smoothprob3d;
    end

    maxdens = max(prob3d(:));
    ncolors = size(g.cmap,1);
    
    prob3d    = round((prob3d-g.cmin)/(g.cmax - g.cmin)*(ncolors-1))+1; % project desnity image into the color space: [1:ncolors]
    prob3d( find(prob3d > ncolors) ) = ncolors;
    prob3d( find(prob3d < 1))        = 1; % added by Makoto
    newprob3d = zeros(size(prob3d,1), size(prob3d,2), size(prob3d,3), 3);

    outOfBrainMask = find(isnan(prob3d)); % place NaNs in a mask, NaNs are assumed for points outside the brain
    prob3d(outOfBrainMask) = 1;
    
    if strcmpi(col, 'abscolor')
        newprob3d = prob3d/ncolors;
    else
        tmp = g.cmap(prob3d,1); newprob3d(:,:,:,1) = reshape(tmp, size(prob3d));
        tmp = g.cmap(prob3d,2); newprob3d(:,:,:,2) = reshape(tmp, size(prob3d));
        tmp = g.cmap(prob3d,3); newprob3d(:,:,:,3) = reshape(tmp, size(prob3d));
    end
    
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
    
