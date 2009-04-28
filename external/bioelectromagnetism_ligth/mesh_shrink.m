function [FV] = mesh_shrink(FV,origin,dist),

% mesh_shrink - implode vertices of mesh by specific distance
%
% FV = mesh_shrink(FV,origin,dist)
%
% FV is a struct with fields:
%
% FV.vertices   - Nx3 matrix of Cartesian vertex coordindates (X,Y,Z)
% FV.faces      - Mx3 matrix of triangulation of FV.vertices
%
% origin        - 1x3 row vector, usually (0,0,0)
%
% dist          - how far to implode the mesh toward the origin;
%                 this distance is relative to current distance from
%                 the origin, not the total distance from the origin.
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:57 $

% Licence:  GNU GPL, no implied or express warranties
% History:  10/2002, Darren.Weber_at_radiology.ucsf.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    xo = origin(1); yo = origin(2); zo = origin(3);
    
    Nvert = size(FV.vertices,1);
    
    fprintf('...mesh implosion...'); tic;
    
    for v = 1:Nvert,
        
        x = FV.vertices(v,1);
        y = FV.vertices(v,2);
        z = FV.vertices(v,3);
        
        % Find direction cosines for line from centre to vertex
        d = sqrt( (x-xo)^2 + (y-yo)^2 + (z-zo)^2 );
        
        l = (x-xo)/d; % cos alpha
        m = (y-yo)/d; % cos beta
        n = (z-zo)/d; % cos gamma
        
        % now decrease d by dist
        d = d - dist;
        
        % locate vertex at this new distance
        x = (l * d) + xo;
        y = (m * d) + yo;
        z = (n * d) + zo;
        
        FV.vertices(v,:) = [ x y z ];
    end
    
    t = toc; fprintf('...done (%5.2f sec)\n',t);
    
return
    
