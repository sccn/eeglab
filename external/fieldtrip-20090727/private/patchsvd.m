function [grid] = patchsvd(cfg, grid);

% This function does something
%

% Copyright (c) 2006, Jan-Mathijs Schoffelen & Robert Oostenveld, F.C. Donders Centre
%
% $Log: not supported by cvs2svn $
% Revision 1.1  2006/09/19 11:01:59  jansch
% first implementation
%


% set the defaults
if ~isfield(cfg, 'patchsvd'),    cfg.patchsvd = 3;          end
if ~isfield(cfg, 'patchsvdnum'), cfg.patchsvdnum = 5;       end
if ~isfield(cfg, 'feedback'),    cfg.feedback = 'text';     end


Ndipoles = size(grid.pos,1);
Ninside  = length(grid.inside);
Nchans   = size(grid.leadfield{grid.inside(1)}, 1);
lfall    = cell(1,Ndipoles);

progress('init', cfg.feedback, 'computing patchsvd');
for dipindx=1:Ninside
  % renumber the loop-index variable to make it easier to print the progress bar
  i = grid.inside(dipindx);
  
  % compute the distance from this dipole to each other dipole
  dist = sqrt(sum((grid.pos-repmat(grid.pos(i,:), [Ndipoles 1])).^2, 2));

  % define the region of interest around this dipole
  sel  = find(dist<=cfg.patchsvd);
  sel  = intersect(sel, grid.inside);
  Nsel = length(sel);

  progress(dipindx/Ninside, 'computing patchsvd %d/%d, Nsel=%d\n', dipindx, Ninside, Nsel);

  % concatenate the leadfield of all dipoles that are inside the ROI into one matrix
  lfr = cell2mat(grid.leadfield(sel(:)'));
  % svd of leadfields of dipoles inside the ROI
  [U,S,V] = svd(lfr);
  
  lfall{i} = U(:,1:cfg.patchsvdnum);
end
progress('close');

% update the leadfields
grid.leadfield = lfall;

