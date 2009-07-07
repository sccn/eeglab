function [vol, cfg] = prepare_singleshell(cfg, mri);

% PREPARE_SINGLESHELL creates a simple and fast method for the MEG forward
% calculation for one shell of arbitrary shape. This is based on a
% correction of the lead field for a spherical volume conductor by a
% superposition of basis functions, gradients of harmonic functions
% constructed from spherical harmonics.
%
% Use as
%   [vol] = prepare_singleshell(cfg)
%   [vol] = prepare_singleshell(cfg, seg)
%
% If you do not use a segmented MRI, the configuration should contain
%   cfg.headshape   = a filename containing headshape, a structure containing a
%                     single triangulated boundary, or a Nx3 matrix with surface
%                     points
%   cfg.spheremesh  = number, to retriangulate the mesh with a sphere (default = 3000)
%                     instead of specifying a number, you can specify 'same' to keep the
%                     vertices of the mesh identical to the original headshape points
%
% The following options are relevant if you use a segmented MRI
%   cfg.smooth      = 'no' or the FWHM of the gaussian kernel in voxels (default = 5)
%   cfg.mriunits    = 'mm' or 'cm' (default is 'mm')
%   cfg.sourceunits = 'mm' or 'cm' (default is 'cm')
%   cfg.threshold   = 0.5, relative to the maximum value in the segmentation
%
% This function implements
%   G. Nolte, "The magnetic lead field theorem in the quasi-static
%   approximation and its use for magnetoencephalography forward calculation
%   in realistic volume conductors", Phys Med Biol. 2003 Nov 21;48(22):3637-52.

% TODO the spheremesh option should be renamed consistently with other mesh generation cfgs
% TODO shape should contain pnt as subfield and not be equal to pnt (for consistency with other use of shape)

% Copyright (C) 2006-2007, Robert Oostenveld
%
% $Log: not supported by cvs2svn $
% Revision 1.21  2009/05/29 12:31:04  roboos
% only convert cfg.headshape from config to struct in case it is present
%
% Revision 1.20  2009/05/25 08:05:18  roboos
% ensure that cfg.headshape is a sturct and not a config object (in case tracking is on)
%
% Revision 1.19  2009/05/18 07:34:04  marvger
% fixed bugs introduced by the previous update
%
% Revision 1.18  2009/05/14 19:22:13  roboos
% consistent handling of cfg.headshape in code and documentation
%
% Revision 1.17  2009/04/01 12:29:26  roboos
% added checkconfig
%
% Revision 1.16  2008/09/22 20:17:43  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.15  2008/08/13 21:02:20  roboos
% use general read_headshape instead of specific subfunctions
%
% Revision 1.14  2008/07/31 16:10:33  roboos
% corrected documentation for default smooth=5
%
% Revision 1.13  2008/04/10 08:03:11  roboos
% renamed the fieldtrip/private/prepare_vol_sens function into prepare_headmodel
%
% Revision 1.12  2007/08/06 09:20:14  roboos
% added support for bti_hs
%
% Revision 1.11  2007/07/26 07:11:47  roboos
% remove double vertices when possible (twice)
% implemented cfg.sphremesh=same
% updated documentation
%
% Revision 1.10  2007/04/19 17:15:15  roboos
% retriangulate headshape to the desired number of vertices
%
% Revision 1.9  2006/08/16 10:52:30  marsie
% fixed use of spm_smooth()
%
% Revision 1.8  2006/08/01 10:31:03  marsie
% fixed bug in using spm_smooth
%
% Revision 1.7  2006/07/27 08:29:38  roboos
% use spm_smooth instead of spm_conv, updated documentation
%
% Revision 1.6  2006/06/08 07:50:53  roboos
% updated the conversion between the source and MRI units (support mm,cm,dm,m for both)
%
% Revision 1.5  2006/06/07 15:50:36  roboos
% changed checktoolbox into hastoolbox
%
% Revision 1.4  2006/04/18 19:04:35  roboos
% changed hard-coded 10000 vertex points for brain surface into cfg.spheremesh with default value of 4000
%
% Revision 1.3  2006/04/05 16:09:48  roboos
% forgot to add cfg to output in previous commit
%
% Revision 1.2  2006/04/05 15:07:37  roboos
% return the configuration as second argument, store the headshape points in the cfg
%
% Revision 1.1  2006/03/21 09:41:46  roboos
% new implementation, mainly copy and paste from prepare_localspheres
%

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

% set the defaults
if ~isfield(cfg, 'smooth');        cfg.smooth = 5;          end % in voxels
if ~isfield(cfg, 'mriunits');      cfg.mriunits = 'mm';     end
if ~isfield(cfg, 'sourceunits'),   cfg.sourceunits = 'cm';  end
if ~isfield(cfg, 'threshold'),     cfg.threshold = 0.5;     end % relative
if ~isfield(cfg, 'spheremesh'),    cfg.spheremesh = 4000;   end % approximate number of vertices in spere

if isfield(cfg, 'headshape') && isa(cfg.headshape, 'config')
  % convert the nested config-object back into a normal structure
  cfg.headshape = struct(cfg.headshape);
end

if nargin>1 && (~isfield(cfg,'headshape') || isempty(cfg.headshape))
  basedonmri       = 1;
  basedonheadshape = 0;
elseif nargin==1 && isfield(cfg,'headshape') && ~isempty(cfg.headshape)
  basedonmri       = 0;
  basedonheadshape = 1;
else
  error('inconsistent configuration, cfg.headshape should not be used in combination with an mri input')
end

if basedonmri
  % obtain the head shape from the segmented MRI
  seg = zeros(mri.dim);
  if isfield(mri, 'gray')
    fprintf('including gray matter in segmentation for brain compartment\n')
    seg = seg | (mri.gray>(cfg.threshold*max(mri.gray(:))));
  end
  if isfield(mri, 'white')
    fprintf('including white matter in segmentation for brain compartment\n')
    seg = seg | (mri.white>(cfg.threshold*max(mri.white(:))));
  end
  if isfield(mri, 'csf')
    fprintf('including CSF in segmentation for brain compartment\n')
    seg = seg | (mri.csf>(cfg.threshold*max(mri.csf(:))));
  end
  if ~strcmp(cfg.smooth, 'no'),
    % check whether the required SPM2 toolbox is available
    hastoolbox('spm2', 1);
    fprintf('smoothing the segmentation with a %d-pixel FWHM kernel\n',cfg.smooth);
    seg = double(seg);
    spm_smooth(seg, seg, cfg.smooth);
  end
  % threshold for the last time
  seg = (seg>(cfg.threshold*max(seg(:))));
  % determine the center of gravity of the segmented brain
  xgrid = 1:mri.dim(1);
  ygrid = 1:mri.dim(2);
  zgrid = 1:mri.dim(3);
  [X, Y, Z] = ndgrid(xgrid, ygrid, zgrid);
  ori(1) = mean(X(seg));
  ori(2) = mean(Y(seg));
  ori(3) = mean(Z(seg));
  pnt = triangulate_seg(seg, cfg.spheremesh, ori);
  pnt(:,4) = 1;
  pnt = (mri.transform * pnt')';
  % convert the MRI surface points into the same units as the source/gradiometer
  scale = 1;
  switch cfg.sourceunits
    case 'mm'
      scale = scale * 1000;
    case 'cm'
      scale = scale * 100;
    case 'dm'
      scale = scale * 10;
    case 'm'
      scale = scale * 1;
    otherwise
      error('unknown physical dimension in cfg.sourceunits');
  end
  switch cfg.mriunits
    case 'mm'
      scale = scale / 1000;
    case 'cm'
      scale = scale / 100;
    case 'dm'
      scale = scale / 10;
    case 'm'
      scale = scale / 1;
    otherwise
      error('unknown physical dimension in cfg.mriunits');
  end
  if scale~=1
    fprintf('converting MRI surface points from %s into %s\n', cfg.sourceunits, cfg.mriunits);
    shape = pnt(:,1:3) * scale;
  else
    shape = pnt(:,1:3);
  end

  % compute triangulation
  shape = unique(shape, 'rows');
  pnt  = shape;
  Npnt = size(pnt,1);
  avg  = mean(pnt, 1);
  pnt  = pnt - repmat(avg, Npnt, 1);  % shift center of points towards the origin
  dist = sqrt(sum(pnt.^2,2));
  pnt  = pnt ./ repmat(dist, 1, 3);   % normalize to a unit sphere
  tri  = convhulln(pnt);              % construct the triangulation using a convex hull
  pnt  = shape;  
  
  fprintf('placed %d points on the brain surface\n', length(shape));
elseif basedonheadshape
  % get the surface describing the head shape
  if isstruct(cfg.headshape) && isfield(cfg.headshape, 'pnt')
    % use the headshape surface specified in the configuration
    headshape = cfg.headshape;
  elseif isnumeric(cfg.headshape) && size(cfg.headshape,2)==3
    % use the headshape points specified in the configuration
    headshape.pnt = cfg.headshape;
  elseif ischar(cfg.headshape)
    % read the headshape from file
    headshape = read_headshape(cfg.headshape);
  else
    error('cfg.headshape is not specified correctly')
  end
  if ~isfield(headshape, 'tri')
    % generate a closed triangulation from the surface points
    headshape.pnt = unique(headshape.pnt, 'rows');
    headshape.tri = projecttri(headshape.pnt);
  end
  
  pnt  = headshape.pnt;
  tri  = headshape.tri;
  
end % basedonmri or basedonheadshape

if nargin<2
  % the triangulation is based on the shape,
  if isequal(cfg.spheremesh, 'same')
    % keep the same triangulation
    tri = projecttri(pnt); % the triangulation is not guaranteed closed, reconstruct it
  else
    [tri1, pnt1] = reducepatch(tri, pnt, 3*cfg.spheremesh);
    % remove double vertices
    pnt1 = unique(pnt1, 'rows');
    % reconstruct the triangulation
    tri1 = projecttri(pnt1);
    % replace the probably unevenly distributed triangulation with a regular one
    % and retriangulate it to the desired accuracy
    [pnt2, tri2] = msphere(cfg.spheremesh); % this is a regular triangulation
    [pnt, tri] = retriangulate(pnt1, tri1, pnt2, tri2, 2);
  end
end

% construct the geometry of the volume conductor model, containing a single boundary
% the initialization of the forward computation code is done later in prepare_headmodel
vol = [];
vol.bnd.pnt = pnt;
vol.bnd.tri = tri;
vol.type    = 'nolte';

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

