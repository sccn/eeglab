% plotmri() - Overlay activation on MNI normalized brain and vizualize it.
%
% Usage:
%  >> mri = plotmri( mri, act, 'key', val, ... );
%
% Inputs:
%      mri   - mri structure normalized to MNI brain. Supposed to contain 
%              a field anatomy and a field transform. 4-D array with color
%              from 0 to 1.
%      act   - activation value, one per voxel.
%
% optional inputs:
%    'mrigamma' - [float] mri gamma factor (alter contrast). 1 does not
%                 change contrast. < 1 increases contrast; > 1 decreases
%                 contrast.
%    'actgamma' - [float] activity gamma factor (alter contrast). 1 does not
%                 change contrast. < 1 increases contrast; > 1 decreases
%                 contrast.
%    'cmap'     - [hot] colormap for activity
%    'plot'     - [X Y Z] plot one MRI slices at MNI coordinate [X Y Z].
%    'scroll'   - ['on'|'off'] GUI for selecting slices. Default is 'on'.
%
% 
% Example: see loreta_importcomp()
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2005

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, 2005 arno@salk.edu
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

function mri = plotmri( mri, activations, varargin)
    
    if nargin < 1
        help plotmri;
        return;
    end;
    if ~isstr(mri)
        
        g = finputcheck(varargin, { ...
                            'actgamma'  'real'      [0 Inf]    1;
                            'mrigamma'  'real'      [0 Inf]    1;
                            'actfactor' 'real'      [0 Inf]    0.4;
                            'transform' 'real'      []         eye(4);
                            'colchan'   'integer'   [0 Inf]    1;
                            'cmap'      'integer'   [0 Inf]    hot;
                            'plot'      'real'      []         [];
                            'scroll'    'string'    { 'on' 'off' } 'on' });
        
        if isstr(g), error(g); end;
        disp('Slice number in MNI coordinates (mm)');
        
        % find point with highest activation value
        % ----------------------------------------
        if isempty(g.plot)
            [tmp ind] = max(activations(:));        
            activations(ind) = 2;
            [i j] = find( activations == 2 );        
            activations(ind) = tmp;
            [sx sy sz] = size(activations);
            g.plot(1) = i;
            g.plot(2) = mod(j-1,sy)+1;
            g.plot(3) = (j-g.plot(2))/sy+1;
            plotinmricoord = g.plot;
            g.plot = g.transform*[g.plot+1 1]';
            g.plot = g.plot(1:3)';
            disp('Finding maximum of activity');
        end;
        
        % set the activation to true color
        % --------------------------------
        g.act    = activations;
        tmpact   = g.act.^g.actgamma;
        ncolors = size(g.cmap,1);
        tmpact    = round(tmpact/max(tmpact(:))*(ncolors-1))+1; % -> range: colors 1:ncolors
        newprob3d = zeros(size(tmpact,1), size(tmpact,2), size(tmpact,3), 3);
        
        for ix = 1:size(newprob3d,1)  % could somehow use matrix ops here???
            for iy = 1:size(newprob3d,2)
                for iz = 1:size(newprob3d,3)
                    newprob3d(ix, iy, iz, :) = g.cmap(tmpact(ix, iy, iz),:);
                end;
            end;
        end;
        
        % compute MRI after gamma factor etc...
        % -------------------------------------
        g.invtransf = pinv(g.transform);
        g.mri    = mri;
        try
            g.curmri = g.mri.^g.mrigamma;
        catch,
            g.curmri = g.mri.anatomy.^g.mrigamma;
        end;            
        
        % make true color if necessary
        % ----------------------------
        if ndims(g.curmri) == 3
            g.curmri(:,:,:,2) = g.curmri(:,:,:,1);
            g.curmri(:,:,:,3) = g.curmri(:,:,:,1);
        end;
        
        % normalize
        % ---------
        newprob3d  = newprob3d / max(newprob3d(:));
        g.curmri   = g.curmri  / max(g.curmri(:));        
        g.curmri   = 0.5*g.curmri + 0.5*newprob3d;
        g.curmri(plotinmricoord(1), plotinmricoord(2), plotinmricoord(3), :) = 0;
        
        % make scrolling buttons
        % ----------------------
        g.fid = figure( 'position', [60 705 1010 335]);
        set(g.fid, 'userdata', g);
        [sx sy sz tmp] = size(g.mri);
        h1 = axes('unit', 'pixel', 'position', [58  57 290 270], 'tag', 'view1');
        h2 = axes('unit', 'pixel', 'position', [358 57 290 270], 'tag', 'view2');
        h3 = axes('unit', 'pixel', 'position', [658 57 290 270], 'tag', 'view3');
        makescroll(h1, g.plot(3), 'view1', 2, g.scroll);
        makescroll(h2, g.plot(2), 'view2', 2, g.scroll);
        makescroll(h3, g.plot(1), 'view3', 2, g.scroll);
        set(gcf, 'color', 'k');
        
        plotmri('redraw', 'all');
        return;
    else
        g = get(gcf, 'userdata');
        options = { 'FaceColor','texturemap', 'EdgeColor','none', 'CDataMapping', ...
                    'direct','tag','img', 'facelighting', 'none' };
        coord(3) = str2num(get(findobj(g.fid, 'tag', 'textview1'), 'string'));
        coord(2) = str2num(get(findobj(g.fid, 'tag', 'textview2'), 'string'));
        coord(1) = str2num(get(findobj(g.fid, 'tag', 'textview3'), 'string'));
        slice = round(g.invtransf*[ coord 1]')';        
        
        % string command (for now just redraw)
        % ------------------------------------
        if strcmpi(activations, 'all') | strcmpi(activations, 'view1'), redraw(1, slice, g); end; 
        if strcmpi(activations, 'all') | strcmpi(activations, 'view2'), redraw(2, slice, g); end; 
        if strcmpi(activations, 'all') | strcmpi(activations, 'view3'), redraw(3, slice, g); end; 
    end;

% redraw different slices
% -----------------------
function redraw(viewnb, slicenb, g);
    
    options = { 'FaceColor','texturemap', 'EdgeColor','none', 'CDataMapping', ...
                'scaled','tag','img', 'facelighting', 'none' };
    [sx sy sz tmp] = size(g.curmri);
    maxs = max([ sx sy sz ]);
    switch viewnb
     case 1, img1 = squeeze(g.curmri(:,:,slicenb(3),:));
             axes(findobj(g.fid, 'tag', 'view1'));
             delete(findobj(g.fid, 'tag', 'surfview1'));
             surface([0 0; sx sx], [0 sy; 0 sy], [0 0; 0 0], img1, options{:}, 'tag', 'surfview1'); 
             xlim([0 maxs]); ylim([0 maxs]); zlim([0 maxs]);
             axis off; axis equal;
             view( [0 0 1] );
     case 2, img2 = squeeze(g.curmri(:,slicenb(2),:,:));
             axes(findobj(g.fid, 'tag', 'view2'));
             delete(findobj(g.fid, 'tag', 'surfview2'));
             surface([0 0; sx sx], [0 0; 0 0], [0 sz; 0 sz], img2, options{:}, 'tag', 'surfview2'); 
             xlim([0 maxs]); ylim([0 maxs]); zlim([0 maxs]);
             axis off; axis equal;
             view( [0 1 0] );
     case 3, img3 = squeeze(g.curmri(slicenb(1),:,:,:));
             axes(findobj(g.fid, 'tag', 'view3'));
             delete(findobj(g.fid, 'tag', 'surfview3'));
             surface([0 0; 0 0], [0 0; sy sy], [0 sz; 0 sz], img3, options{:}, 'tag', 'surfview3'); 
             xlim([0 maxs]); ylim([0 maxs]); zlim([0 maxs]);
             axis off; axis equal;
             view( [1 0 0] );
    end;            
    return;
    
    
% make GUI for each plot
% ----------------------
function makescroll(fid, curcoords, tag, coordinc, visible)
    
    set(fid, 'unit', 'normalized');
    pos = get(fid, 'position');
    set(fid, 'tag', tag);
    s = [pos(1) pos(2) 0 0]+[0.04 0 0 0];
    q = [pos(3) pos(4) pos(3) pos(4)];
    ht = 0.05;
    wd = 0.03;
    
    coordminus = [0.03 -0.05 wd   ht]+s;
    coordtext  = [0.08 -0.05 0.04 ht]+s;
    coordcom   = [0.11 -0.05 wd   ht]+s;
    coordplus  = [0.15 -0.05 wd   ht]+s;
    
    cb_minus = [ 'tmpobj = findobj(gcbf, ''tag'', ''text' tag ''');' ...
                 'tmpval = str2num(get(tmpobj,''string'')) - ' num2str(coordinc) ';' ...
                 'set(tmpobj, ''string'', num2str(tmpval));' ...
                 'clear tmpval, tmpobj;' ...
                 'plotmri(''update'', ''' tag ''');' ];
    cb_plus  = [ 'tmpobj = findobj(gcbf, ''tag'', ''text' tag ''');' ...
                 'tmpval = str2num(get(tmpobj,''string'')) + ' num2str(coordinc) ';' ...
                 'set(tmpobj, ''string'', num2str(tmpval));' ...
                 'clear tmpval, tmpobj;' ...
                 'plotmri(''update'', ''' tag ''');' ];
    
    % make buttons
    % ------------
    h(1) = uicontrol( 'unit', 'normalized', 'style', 'edit', 'tag', [ 'text' tag], 'position', ...
               coordtext, 'string', num2str(curcoords), 'callback', ['plotmri(''update'', ''' tag ''');' ]);
    h(2) = uicontrol( 'unit', 'normalized', 'style', 'text', 'position', ...
               coordcom, 'string', 'mm', 'backgroundcolor', 'k', 'foregroundcolor', 'w' );
    h(3) = uicontrol( 'unit', 'normalized', 'style', 'pushbutton', 'position', ...
               coordminus, 'string', '-', 'callback', cb_minus);
    h(4) = uicontrol( 'unit', 'normalized', 'style', 'pushbutton', 'position', ...
               coordplus, 'string', '+', 'callback', cb_plus);
    set(h, 'visible', visible);
