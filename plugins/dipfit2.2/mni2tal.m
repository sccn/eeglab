function outpoints = mni2tal(inpoints)

% MNI2TAL - MNI to Talairach coordinates (best guess)
% 
% outpoints = mni2tal(inpoints)
% 
% inpoints - Nx3 or 3xN matrix of coordinates
%            (N being the number of points)
% 
% outpoints - the coordinate matrix with Talairach points
% 
% See also, TAL2MNI, MNI2TAL_MATRIX &
% http://www.mrc-cbu.cam.ac.uk/Imaging/mnispace.html
% 

% $Revision: 1.1 $ $Date: 2009-01-30 02:48:24 $

% Licence:  GNU GPL, no express or implied warranties
% Matthew Brett 10/8/99, matthew.brett@mrc-cbu.cam.ac.uk
% modified 02/2003, Darren.Weber_at_radiology.ucsf.edu
%                   - removed dependence on spm_matrix and
%                     abstracted the matrix in case it changes
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dimdim = find(size(inpoints) == 3);
if isempty(dimdim),
    error('input must be a Nx3 or 3xN matrix')
end
if dimdim == 2,
    inpoints = inpoints';
end

% Transformation matrices, different zooms above/below AC
M2T = mni2tal_matrix;

inpoints = [inpoints; ones(1, size(inpoints, 2))];

tmp = inpoints(3,:) < 0;  % 1 if below AC

inpoints(:,  tmp) = (M2T.rotn * M2T.downZ) * inpoints(:,  tmp);
inpoints(:, ~tmp) = (M2T.rotn * M2T.upZ  ) * inpoints(:, ~tmp);

outpoints = inpoints(1:3, :);
if dimdim == 2,
    outpoints = outpoints';
end

return
