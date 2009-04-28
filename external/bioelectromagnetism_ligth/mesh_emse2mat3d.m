function [faces,vertices] = mesh_emse2mat3d(vertex,patch)

% mesh_emse2mat3d - Convert emse vertex/patches to matlab vertices/faces
%
% Useage: [faces,vertices] = mesh_emse2mat3d(vertex,patch)
%
%           vertex & patch are generated from emse files with
%           mesh_emse2matlab (see this file for more help).
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:57 $

% Licence:  GNU GPL, no implied or express warranties
% History:  12/98 Abbas Kouzani
%           09/01 Darren.Weber_at_radiology.ucsf.edu
%                 - converted to function, rather than script
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%S=sprintf('load %s',file);eval(S);

vertices = zeros(length(vertex),3);
for i=1:length(vertex),
    vertices(i,1) = vertex(i).location1;
    vertices(i,2) = vertex(i).location2;
    vertices(i,3) = vertex(i).location3;
end

faces = zeros(length(patch),3);
for i=1:length(patch),

    v = patch(i).vertex1;   flag = 0;
    for j=1:length(vertex),
        if isequal(v, vertex(j).address)
            faces(i,1) = j; flag = 1; break;
        end
    end
    if (flag==0) fprintf('...error in patch(%d).vertex1!\n',i); break;end
    
    v = patch(i).vertex2;   flag = 0;
    for j=1:length(vertex),
        if isequal(v, vertex(j).address)
            faces(i,2) = j; flag = 1; break;
        end
    end
    if (flag==0) fprintf('...error in patch(%d).vertex2!\n',i); break;end
    
    v = patch(i).vertex3;   flag = 0;
    for j=1:length(vertex),
        if isequal(v, vertex(j).address)
            faces(i,3) = j; flag = 1; break;
        end
    end
    if (flag==0) fprintf('...error in patch(%d).vertex3!\n',i); break;end
end

%disp('...Saving');
%GET1=input(' Enter output filename: ','s');
%S=sprintf('save %s faces vertices',GET1);
%eval(S);
