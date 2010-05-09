% colorizemri() - Overlay activation on MNI normalized brain. 
%
% Usage:
%  >> [anat act mri] = colorizemri( mri, locations, act, 'key', val, ... );
%
% Inputs:
%      mri   - mri structure normalized to MNI brain. Supposed to contain 
%              a field anatomy and a field transform.\
%      locations - 3-D locations (one per row) of the center of activated 
%                  activated voxels in MNI coordinates
%      act       - activation value, one per voxel
%
% optional inputs:
%    'voxsize'  - [float] voxel size (in MNI voxels (mm)). Default is 1.
%    'smooth'   - [integer] optional 3-D smoothing (in MNI voxels (mm)).
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
% Outputs:
%      anat  - mri anatomical 4-D array.
%      act   - activation 4-D array.
%      mri   - mri structure with a new field "anatomycol" containing
%              true color voxel values.
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

function [ anat, actarray, mri ] = colorizemri( mri, locations, activations, varargin)
    
    if nargin < 3
        help colorizemri;
        return;
    end;
    g = finputcheck(varargin, { 'smooth'    'integer'   [0 Inf]    0;
                                'actgamma'  'real'      [0 Inf]    1;
                                'mrigamma'  'real'      [0 Inf]    1;
                                'actfactor' 'real'      [0 Inf]    0.4;
                                'voxsize'   'real'      [0 Inf]    1;
                                'colchan'   'integer'   [0 Inf]    1;
                                'plot'      'real'      []         [] });
    if isstr(g), error(g); end;
    
    % find location in the MRI itself
    % -------------------------------
    invtransf = pinv( mri.transform );
    mrilocsbeg = round(invtransf*[ locations-g.voxsize/2 ones(size(locations,1),1)]')';
    mrilocsend = round(invtransf*[ locations+g.voxsize/2 ones(size(locations,1),1)]')';
    tmplocsbeg = mrilocsbeg(:,1:3);
    tmplocsend = mrilocsend(:,1:3);
    mrilocsbeg = min(tmplocsbeg, tmplocsend);
    mrilocsend = max(tmplocsbeg, tmplocsend);
    
    % scale activations from 0 to 1
    % -----------------------------
    activations = activations/max(activations);
    
    % set activation at selected points
    % ---------------------------------
    actarray = zeros(size(mri.anatomy));
    for index = 1:size(mrilocsbeg,1)
        actarray(mrilocsbeg(index,1):mrilocsend(index,1), ...
                 mrilocsbeg(index,2):mrilocsend(index,2), ...
                 mrilocsbeg(index,3):mrilocsend(index,3)) = activations(index);
    end;
    
    % smoothing 3x3
    % -------------
    if g.smooth
        disp('Smoothing...');
        actarray = smooth3(actarray, 'box', g.smooth);
    end;
    
    % make MRI true color
    % -------------------
    mri.anatomycol = zeros([ size( mri.anatomy ) 3]);
    mri.anatomycol(:,:,:,1) = mri.anatomy;
    mri.anatomycol(:,:,:,2) = mri.anatomy;
    mri.anatomycol(:,:,:,3) = mri.anatomy;
    mri.anatomycol = mri.anatomycol/max(mri.anatomycol(:));
    mri.anatomycol = mri.anatomycol.^g.mrigamma;
    anat = mri.anatomycol;
    
    % overlay activation with MRI
    % ---------------------------
    actarray = actarray.^g.actgamma;
    for ic = g.colchan
        mri.anatomycol(:,:,:,ic) = mri.anatomycol(:,:,:,ic) + actarray*g.actfactor;
    end;
    
    % plot
    % ----
    if ~isempty(g.plot)
        figure( 'position', [60 705 1366 395]);
        options = { 'FaceColor','texturemap', 'EdgeColor','none', 'CDataMapping', ...
                    'direct','tag','img', 'facelighting', 'none' };
        slicecoord = round(invtransf*[ g.plot 1]')';
        img1 = squeeze(mri.anatomycol(:,:,slicecoord(3),:));
        img2 = squeeze(mri.anatomycol(:,slicecoord(2),:,:));
        img3 = squeeze(mri.anatomycol(slicecoord(1),:,:,:));
        [sx sy sz tmp] = size(mri.anatomycol);
        subplot(1,3,1); surface([0 0; sx sx], [0 sy; 0 sy], [0 0; 0 0], img1, options{:}); axis off; axis equal;
        subplot(1,3,2); surface([0 0; sx sx], [0 sz; 0 sz], [0 0; 0 0], img2, options{:}); axis off; axis equal;
        subplot(1,3,3); surface([0 0; sy sy], [0 sz; 0 sz], [0 0; 0 0], img3, options{:}); axis off; axis equal;
        set(gcf, 'color', 'k');
    end;
    return;
    
    
