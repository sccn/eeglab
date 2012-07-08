% loreta_importcomp() - Overlay activation on MNI normalized brain. 
%
% Usage:
%  >> mri = loreta_importcomp( mrifile, actfile, 'key', 'val', ... );
%
% Inputs:
%      mri   - mri structure normalized to MNI brain. Supposed to contain 
%              a field anatomy and a field transform.
%      act       - activation value, one per voxel
%
% optional inputs:
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
%    'filelocs' - [string] location file for LORETA voxels in MNI space.
%                 This is contained in the file 'LORETAtalairach.xyz' from
%                 the LORETA distribution. This is also distributed with 
%                 this plugin.
%
% Outputs:
%      mri   - mri structure with a new field "anatomycol" containing
%              true color voxel values.
% 
% Example: 
% [anat act mri] = loreta_importcomp( [ '/data/common/matlab/eeglab/plugins/dipfit2.0/' ...
%                             'standard_BEM/standard_mri.mat' ], 'LORETAloc.lor3');
% plotmri(anat, act, 'transform', mri.transform);
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

function [mri, act, mristruct] = loreta_importcomp(mrifile, fileact, varargin)
    
    if nargin < 2
        help loreta_importcomp
        return;
    end;
    g = finputcheck(varargin, { 'smooth'    'integer'   [0 Inf]    0;
                                'actgamma'  'real'      [0 Inf]    1;
                                'mrigamma'  'real'      [0 Inf]    1;
                                'filelocs'  'string'    []         'LORETAtalairach.xyz';
                                'actfactor' 'real'      [0 Inf]    0.4 });
    if isstr(g), error(g); end;
    
    % get locations
    % -------------
    fid = fopen(g.filelocs, 'r'); % should be at the same location as
                                             % this file
    if fid == -1, error([ 'Cannot find or open ''' g.filelocs ''' file' ]); end;
    nlines = fscanf(fid, '%d',1 ); % should be equal to 2394
    locations = fscanf(fid, '%d', [2394*3]);
    locations = reshape(locations, [3 2394])';
    fclose(fid);
    
    % get values
    % ----------
    %fid = fopen('LORETAloc.lor3', 'r')
    fid = fopen(fileact, 'r');
    if fid == -1, error([ 'Cannot find or open ''' fileact ''' file' ]); end;
    Xact = fread(fid, 2394, 'float');
    Yact = fread(fid, 2394, 'float');
    Zact = fread(fid, 2394, 'float');
    activations = sqrt(Xact.^2 + Yact.^2 + Zact.^2);
    fclose(fid);
    
    % read MRI if necessary
    % ---------------------
    if isstr(mrifile)
        try, 
            mrifile = load('-mat', mrifile);
            mrifile = mrifile.mri;
        catch,
            warning backtrace off;
            mrifile = read_fcdc_mri(mrifile);
            warning backtrace on;
        end;
    end;

    % colorize MRI
    % ------------
    [mri act mristruct] = colorizemri( mrifile, locations, activations, 'voxsize', 7, 'smooth', g.smooth, ...
                       'actfactor', g.actfactor, 'actgamma', g.actgamma, 'mrigamma', g.mrigamma);
    
    
