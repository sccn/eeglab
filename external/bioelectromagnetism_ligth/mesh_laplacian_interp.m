function [int, keepindex, repindex] = mesh_laplacian_interp(lap, index);

% mesh_laplacian_interp - Computes the zero Laplacian interpolation matrix
% 
% Useage: [int, keepindex, repindex] = mesh_laplacian_interp(lap, index)
% 
% This function calculates an interpolation matrix that provides
% the coefficients for the calculation of potential values at all
% unknown vertices of a mesh, given known potential values at
% a subset of the mesh vertices (at 'index').  The interpolation
% solution is constrained by a minimal norm of the Laplacian
% of the mesh.  See the reference below for details.
% 
% 'lap' is the laplacian matrix for the full mesh (see mesh_laplacian)
% 'int' is the matrix which interpolates from the points in 'index'
% to the full mesh.  'index' is a row vector of indices into a 
% subset of the vertices used to calculate 'lap'.  This subset 
% is where the electric potential is known and usually corresponds 
% to the given electrode vertices, eg:
% 
% index = dsearchn(scalpvert,elecvert)';
% 
% If 'index' contains repeated indices, this is a serious problem
% but this function will try to continue with only the unique indices. 
% The 'keepindex' array can be used to select these.
% The 'repindex' array is the repeated indices.
% 
% Interpolations can be done using matrix 'int', eg:
% 
% [int, keepindex, repindex] = mesh_laplacian_interp(lap,index);
% if isempty(repindex),
%   Vint = int * Vknown;
% else
%   Vint = int * Vknown(keepindex);
% end
% 
% This attempts to implement interpolation method B (p. 336) of 
% Oostendorp T, Oosterom A & Huiskamp G (1989),
% Interpolation on a triangulated 3D surface.
% Journal of Computational Physics, 80: 331-343.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:57 $

% Licence:  GNU GPL, no implied or express warranties
% History:  (c) 04/2002 Robert Oostenveld
%           - agreed to release 'lapint' under GNU GPL
%           04/2002, Darren.Weber_at_radiology.ucsf.edu
%           - introduced check for index replications and
%             adjusted calculations/output accordingly
%           - converted lap to sparse matrix and solution
%             of interpolation matrix with \ operator
%           - accepts sparse lap input and returns sparse int
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(lap,1)~=size(lap,2), error('MESH_LAPLACIAN_INTERP: lap matrix is not square'); end

if issparse(lap), lap = full(lap); end

tic

% -- Check for any duplicate values in the known index

if size(index,1)~=1, index = index'; end

% Remove any replicate indices from 'index'
[KnownIndex, i, i] = union(index,index);
if ~isequal(length(KnownIndex),length(index)),
    fprintf('\a\n\nMESH_LAPLACIAN_INTERP: Warning:\n\nTrimming duplicate values from index, use keepindex!\n\n');
end
keepindex = sort(i);
repindex = setdiff(1:length(index),sort(i));

KnownIndex = index(keepindex); % unsort KnownIndex
clear index

k = length(KnownIndex);
n = length(lap);

fprintf('MESH_LAPLACIAN_INTERP: Calc Interpolation matrix for %5d to %5d vertices...',k,n);

% find 'unknown' indices of lap matrix
UnknownIndex = setdiff(1:n, KnownIndex);

% reshuffle rows & columns of lap matrix
lapi = [KnownIndex, UnknownIndex];
lap = lap(lapi, :); % rows
lap = lap(:, lapi); % columns

% Segregate known/unknown portions of lap
k = length(KnownIndex);
n = length(lap);

L11 = lap(1:k    ,1:k    );
L12 = lap(1:k    ,(k+1):n);
L21 = lap((k+1):n,1:k    );
L22 = lap((k+1):n,(k+1):n);

clear lap lapi; % tidy up some memory

% Convert to sparse for quicker computation
A = sparse([L12; L22]);
B = sparse([L11; L21]);

clear L11 L12 L21 L22;  % tidy up some memory

%int = -pinv(A) * B;  % cannot use pinv with sparse matrix
int = -A \ B;

% Convert result back to full matrix
int = full(int);

% append the interpolating piece to the identity matrix 
% these take care of the known potentials
int = [eye(k); int];

% reshuffle the columns of the interpolating matrix
[tmp, order] = sort(KnownIndex);
int = int(:,order);

% reshuffle the rows of the interpolating matrix
[tmp, order] = sort([KnownIndex, UnknownIndex]);
int = int(order, :);

int = sparse(int);

t = toc;
fprintf('done (%6.2f sec).\n',t);

return
