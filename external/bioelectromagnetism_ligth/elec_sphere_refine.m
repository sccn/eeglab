function [ p ] = elec_sphere_refine(p)

% elec_sphere_refine - creates spherical vertices between input vertices
%
% [ p ] = elec_sphere_refine( p )
%
% p is the eeg_toolbox struct (see eeg_toolbox_defaults)
% 
% This function works with the p.elec.data struct,
% especially the spherical projections created by 
% elec_sphere_project.
% 
% The spherical electrode positions are triangulated with
% convhulln.  Then, for each face, 4 new vertices are created 
% at the triangle edge midpoints and the triangle centroid.  Each
% face is divided into 6 faces.
% 
% The vertices created in the plane of each face are then projected
% to the spherical surface.
%
% See also mesh_refine, mesh_laplacian, mesh_laplacian_interp
% 


% This can be done until some minimal distance (D) of the mean 
% distance between vertices of all triangles is achieved.  If
% no D argument is given, the function refines the mesh once.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:55 $

% Licence:  GNU GPL, no implied or express warranties
% History:  05/2002, Darren.Weber_at_radiology.ucsf.edu, created
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Check the input struct
if ~exist('p','var'),
   [p] = elec_open;
end


tic;
fprintf('ELEC_SPHERE_REFINE: processing faces...');

% NOTE
% The centroid is located one third of the way from each vertex to 
% the midpoint of the opposite side. Each median divides the triangle 
% into two equal areas; all the medians together divide it into six 
% equal parts, and the lines from the median point to the vertices 
% divide the whole into three equivalent triangles.

V = [p.elec.data.Xsp, p.elec.data.Ysp, p.elec.data.Zsp];
F = convhulln(V);

Zmin = min(V(:,3));

% Initialise a new vertices and faces matrix
Nvert = size(V,1);
Nface = size(F,1);
V2 = zeros(Nface*4,3);
F2 = zeros(Nface*6,3);

for r=1:size(F,1),
    
    % Get the triangle vertices
    A = V(F(r,1),:);
    B = V(F(r,2),:);
    C = V(F(r,3),:);
    
    % Now find the midpoint between all vertices
    ABmid = (A + B) ./ 2;
    BCmid = (B + C) ./ 2;
    CAmid = (C + A) ./ 2;
    
    % Now find the centroid length of the medians
    Ac = BCmid + ( (A - BCmid) ./3 );
    %Bc = CAmid + ( (B - CAmid) ./3 );  % Bc = Ac
    %Cc = ABmid + ( (C - ABmid) ./3 );  % Cc = Ac
    
    % Store the midpoints and the centroid vertices
    NABmid = Nvert + (r*4-3);
    V2(r*4-3,:) = ABmid;
    NBCmid = Nvert + (r*4-2);
    V2(r*4-2,:) = BCmid;
    NCAmid = Nvert + (r*4-1);
    V2(r*4-1,:) = CAmid;
    Nc     = Nvert + (r*4-0);
    V2(r*4-0,:) = Ac;
    
    % Create new faces for centroid
    F2(r*6-5,:) = [F(r,1), NABmid, Nc     ];
    F2(r*6-4,:) = [F(r,1), Nc    , NCAmid ];
    F2(r*6-3,:) = [F(r,2), Nc    , NABmid ];
    F2(r*6-2,:) = [F(r,2), NBCmid, Nc     ];
    F2(r*6-1,:) = [F(r,3), NCAmid, Nc     ];
    F2(r*6-0,:) = [F(r,3), Nc    , NBCmid ];
    
    %figure; patch('vertices',[A;B;C],'faces',[1 2 3],'facecolor',[.7 .7 .7]); hold on;
    %plot3(A(1),A(2),A(3),'ro');
    %plot3(BCmid(1),BCmid(2),BCmid(3),'ro');
    %plot3(Ac(1),Ac(2),Ac(3),'bo')
    %plot3(Bc(1),Bc(2),Bc(3),'bd')
    %plot3(Cc(1),Cc(2),Cc(3),'b.')
    %if isequal(r,2), return; end
    
end

% Add the new vertices to the old ones
V = [V;V2];

% All of the above is plane geometry,
% So now project all points to the sphere
V = proj_sph(V,p.elec.data.Rsp(1),p.elec.data.centroid);

% Reset all vertices below Zmin to Zmin
% This is a crude modification - needs work
Vlow = find(V(:,3) < Zmin);
V(Vlow,3) = Zmin;

% Store the results in new elements of p.elec.data
p.elec.data.Vsp = V;
p.elec.data.Fsp = F2;


t=toc;
fprintf('...done (%5.2f sec)\n',t);

return


% Find the length of each median
%A2BClen = sqrt ( sum( (A - BCmid).^2, 2 ) );
%B2CAlen = sqrt ( sum( (B - CAmid).^2, 2 ) );
%C2ABlen = sqrt ( sum( (C - ABmid).^2, 2 ) );


function V = proj_sph(v,r,c)
    
    % Find the projection point of X,Y,Z to the fitted sphere radius r
    % Cartesian inputs:
    % v is the vertex matrix
    % r is the sphere radius
    % c is the sphere centroid
    
    X = v(:,1);
    Y = v(:,2);
    Z = v(:,3);
    
    xo = c(1);
    yo = c(2);
    zo = c(3);
        
    % Convert Cartesian X,Y,Z to spherical (radians)
    theta = atan2( (Y-yo), (X-xo) );
    phi = atan2( sqrt( (X-xo).^2 + (Y-yo).^2 ), (Z-zo) );
    % do not recalc: r = sqrt( (X-xo).^2 + (Y-yo).^2 + (Z-zo).^2);
	
	%   Recalculate X,Y,Z for constant r, given theta & phi.
	R = ones(size(phi)) * r;    
	x = R .* sin(phi) .* cos(theta);
	y = R .* sin(phi) .* sin(theta);
	z = R .* cos(phi);

    V = [x y z];
    
return
