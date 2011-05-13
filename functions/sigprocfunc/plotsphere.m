% plotsphere() - This function is used to plot a sphere and
%                project them onto specific surfaces. This may
%                be used for plotting dipoles for instance.
%
% Usage:
%   >> handle = plotsphere(pos, rad, 'key', 'val')
%
% Inputs: 
%   pos       - [x y z] 3-D position of the sphere center
%   rad       - [real] sphere radius
%
% Optional inputs:
%   'nvert'   - number of vertices. Default is 15.
%   'color'   - sphere color. Default is red.
%   'proj'    - [x y z] project sphere to 3-D axes. Enter NaN for not
%               projecting. For instance [-40 NaN -80] will project
%               the sphere on the y/z plane axis at position x=-40 and on 
%               the x/y plane at position z=-80.
%   'projcol' - color of projected spheres
%   'colormap' - [real] sphere colormap { default: jet }
%
% Output:
%   handle    - sphere graphic handle(s). If projected sphere are ploted
%               handles of plotted projected spheres are also returned.
%
% Example:
% figure; plotsphere([3 2 2], 1, 'proj', [0 0 -1]);
% axis off; axis equal; lighting phong; camlight left; view([142 48])
%
% Authors: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 2004 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) Arnaud Delorme, 2004
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

function [handles] = plotsphere(pos, rad, varargin);
    
    if nargin < 2
        help plotsphere;
        return;
    end;
    
    g = finputcheck(varargin, { 'color'    { 'real','string' } []         [1 0 0];
                                'nvert'    'integer'           [2 Inf]    15;
                                'proj'     'real'              []         [];
                                'colormap' 'real'              []         jet(64);
                                'projcol'  { 'real','string' } []         [0 0 0] }, 'plotsphere');
    if isstr(g), error(g); end;
    
    % decode color if necessary
    % -------------------------
    if ~isstr(g.color) & length(g.color) == 1
        g.color = g.colormap(g.color,:);
    elseif isstr(g.color)
        g.color = strcol2real(g.color);
    end;
    if ~isstr(g.projcol) & length(g.projcol) == 1
        g.projcol = g.colormap(g.projcol,:);
    elseif isstr(g.projcol)
        g.projcol = strcol2real(g.projcol);
    end;
            
    % ploting sphere
    % ==============
    [xstmp ystmp zs] = sphere(g.nvert);
    l=sqrt(xstmp.*xstmp+ystmp.*ystmp+zs.*zs);
    normals = reshape([xstmp./l ystmp./l zs./l],[g.nvert+1 g.nvert+1 3]);
    xs = pos(1) + rad*ystmp;
    ys = pos(2) + rad*xstmp;
    zs = pos(3) + rad*zs;
    colorarray = repmat(reshape(g.color , 1,1,3), [size(zs,1) size(zs,2) 1]);
    hold on;
    handles = surf(xs, ys, zs, colorarray, 'tag', 'tmpmov', 'EdgeColor','none', 'VertexNormals', normals, ...
                   'backfacelighting', 'lit', 'facelighting', 'phong', 'facecolor', 'interp', 'ambientstrength', 0.3);
    %axis off; axis equal; lighting phong; camlight left; rotate3d

    % plot projections
    % ================
    if ~isempty(g.proj)
        colorarray  = repmat(reshape(g.projcol, 1,1,3), [size(zs,1) size(zs,2) 1]);
        if ~isnan(g.proj(1)), handles(end+1) = surf(g.proj(1)*ones(size(xs)), ys, zs, colorarray, ...
                                                    'edgecolor', 'none', 'facelighting', 'none'); end;
        if ~isnan(g.proj(2)), handles(end+1) = surf(xs, g.proj(2)*ones(size(ys)), zs, colorarray, ...
                                                    'edgecolor', 'none', 'facelighting', 'none'); end;
        if ~isnan(g.proj(3)), handles(end+1) = surf(xs, ys, g.proj(3)*ones(size(zs)), colorarray, ...
                                                    'edgecolor', 'none', 'facelighting', 'none'); end;
    end;
    
function color = strcol2real(color)
    switch color
     case 'r', color = [1 0 0];
     case 'g', color = [0 1 0];
     case 'b', color = [0 0 1];
     case 'c', color = [0 1 1];
     case 'm', color = [1 0 1];
     case 'y', color = [1 1 0];
     case 'k', color = [0 0 0];
     case 'w', color = [1 1 1];
     otherwise, error('Unknown color'); 
    end;
