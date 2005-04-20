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
%    'actfactor' - [float] activity factor from 0 to 1. Default is 0.4.
%                 A factor too high might cause the voxel color to go over
%                 the color limits.
%    'plot'      - [X Y Z] plot one MRI slices at MNI coordinate [X Y Z].
%
% 
% Example: see loreta_importcomp()
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $

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
                            'plot'      'real'      []         [] });
        if isstr(g), error(g); end;
    
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
            g.plot = g.transform*[g.plot+1 1]';
            g.plot = g.plot(1:3)';
            disp('Finding maximum of activity');
        end;
        
        % compute MRI after gamma factor etc...
        % -------------------------------------
        g.invtransf = pinv(g.transform);
        g.mri    = mri;
        g.act    = activations;
        g.curmri = g.mri.^g.mrigamma;
        tmpact   = g.act.^g.actgamma;
        for ic = g.colchan
            g.curmri(:,:,:,ic) = g.curmri(:,:,:,ic) + tmpact*g.actfactor;
        end;
       
        % make scrolling buttons
        % ----------------------
        g.fid = figure( 'position', [60 705 1366 395]);
        set(g.fid, 'userdata', g);
        [sx sy sz tmp] = size(g.mri);
        subplot(1,3,1); set(gca, 'tag', 'view1'); makescroll(gca, g.plot(3), 'view1', 2);
        subplot(1,3,2); set(gca, 'tag', 'view2'); makescroll(gca, g.plot(2), 'view2', 2);
        subplot(1,3,3); set(gca, 'tag', 'view3'); makescroll(gca, g.plot(1), 'view3', 2);
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
                'direct','tag','img', 'facelighting', 'none' };
    [sx sy sz tmp] = size(g.curmri);
    switch viewnb
     case 1, img1 = squeeze(g.curmri(:,:,slicenb(3),:));
             axes(findobj(g.fid, 'tag', 'view1'));
             delete(findobj(g.fid, 'tag', 'surfview1'));
             surface([0 0; sx sx], [0 sy; 0 sy], [0 0; 0 0], img1, options{:}, 'tag', 'surfview1'); 
             axis off; axis equal;
     case 2, img2 = squeeze(g.curmri(:,slicenb(2),:,:));
             axes(findobj(g.fid, 'tag', 'view2'));
             delete(findobj(g.fid, 'tag', 'surfview2'));
             surface([0 0; sx sx], [0 sz; 0 sz], [0 0; 0 0], img2, options{:}, 'tag', 'surfview2'); 
             axis off; axis equal;
     case 3, img3 = squeeze(g.curmri(slicenb(1),:,:,:));
             axes(findobj(g.fid, 'tag', 'view3'));
             delete(findobj(g.fid, 'tag', 'surfview3'));
             surface([0 0; sy sy], [0 sz; 0 sz], [0 0; 0 0], img3, options{:}, 'tag', 'surfview3'); 
             axis off; axis equal;
    end;            
    return;
    
    
% make GUI for each plot
% ----------------------
function makescroll(fid, curcoords, tag, coordinc)
    
    pos = get(fid, 'position');
    set(fid, 'tag', tag);
    s = [pos(1) pos(2) 0 0];
    q = [pos(3) pos(4) pos(3) pos(4)];
    ht = 0.05;
    wd = 0.03;
    
    coordminus = [0.05 -0.05 wd   ht]+s;
    coordtext  = [0.10 -0.05 0.04 ht]+s;
    coordcom   = [0.13 -0.05 wd   ht]+s;
    coordplus  = [0.17 -0.05 wd   ht]+s;
    
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
    uicontrol( 'unit', 'normalized', 'style', 'edit', 'tag', [ 'text' tag], 'position', ...
               coordtext, 'string', num2str(curcoords), 'callback', ['plotmri(''update'', ''' tag ''');' ]);
    uicontrol( 'unit', 'normalized', 'style', 'text', 'position', ...
               coordcom, 'string', 'mm', 'backgroundcolor', 'k', 'foregroundcolor', 'w' );
    uicontrol( 'unit', 'normalized', 'style', 'pushbutton', 'position', ...
               coordminus, 'string', '-', 'callback', cb_minus);
    uicontrol( 'unit', 'normalized', 'style', 'pushbutton', 'position', ...
               coordplus, 'string', '+', 'callback', cb_plus);
    axis off;
