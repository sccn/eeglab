function FV = mesh_smooth(FV,origin,attract),

% mesh_smooth - smooth vertices of a mesh
%
% FV = mesh_smooth(FV,origin,attract)
% 
% FV is a struct with fields:
%
% FV.vertices   - Nx3 matrix of Cartesian vertex coordindates (X,Y,Z)
% FV.faces      - Mx3 matrix of triangulation of FV.vertices
%
% origin        - 1x3 row vector, usually (0,0,0)
%
% attract       - how much to align the mesh vertices, percent 0:1;
%                 currently not implemented, effectively 1
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:58 $

% Licence:  GNU GPL, no implied or express warranties
% History:  10/2002, Darren.Weber_at_radiology.ucsf.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~isfield(FV,'edge'),
    FV.edge = mesh_edges(FV);
end

Nvert = size(FV.vertices,1);

xo = origin(1); yo = origin(2); zo = origin(3);

fprintf('...mesh smoothing...'); tic;

% Check every vertex
for v = 1:Nvert,
    
    x = FV.vertices(v,1);
    y = FV.vertices(v,2);
    z = FV.vertices(v,3);
    
    % Find direction cosines for line from centre to vertex
    d = sqrt( (x-xo)^2 + (y-yo)^2 + (z-zo)^2 );
    l = (x-xo)/d; % cos alpha
    m = (y-yo)/d; % cos beta
    n = (z-zo)/d; % cos gamma
    
    % Calc distance of neighbour vertices
    vi = find(FV.edge(v,:));  % the indices of the neighbours
    % remove duplicate vertices
    %[vert, i, j] = unique(FV.vertices(vi,:),'rows');
    X = FV.vertices(vi,1);
    Y = FV.vertices(vi,2);
    Z = FV.vertices(vi,3);
    D = sqrt( (X-xo).^2 + (Y-yo).^2 + (Z-zo).^2 );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Maybe check for outliers in D at this point
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Dmean = mean(D,1);
    
    % now modify d above using Dmean of neighbour vertices
    if     d > Dmean,  d = d - (d - Dmean);
    elseif d < Dmean,  d = d + (Dmean - d);
    end
    
    % locate vertex at the new distance
    x = (l * d) + xo;
    y = (m * d) + yo;
    z = (n * d) + zo;
    FV.vertices(v,:) = [ x y z ];
    
end

t = toc; fprintf('done (%5.2f sec).\n',t);

return
