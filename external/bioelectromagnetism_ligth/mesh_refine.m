function [ FV ] = mesh_refine(FV,Nface)

% mesh_refine - creates smaller triangles from a triangle mesh
%
% [ FV ] = mesh_refine( FV, Nface )
%
% FV.vertices   - vertex matrix (Nx3)
% FV.faces      - face matrix (Mx3), indices into vertex matrix rows
%
% Nface         - subdivide faces into 4 or 6 faces,
%                 the default 4 provides an even subdivision
% 
% This function calls mesh_refine_tri4 or mesh_refine_tri6.  See
% these for more details.
% 


% This can be done until some minimal distance (D) of the mean 
% distance between vertices of all triangles is achieved.  If
% no D argument is given, the function refines the mesh once.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:57 $

% Licence:  GNU GPL, no implied or express warranties
% History:  08/2002, Darren.Weber_at_radiology.ucsf.edu, created
%                    adapted this function as a wrapper to
%                    mesh_refine_tri4 & mesh_refine_tri6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~exist('FV','var'),
    error('MESH_REFINE: NO input FV struct');
elseif isempty(FV),
    error('MESH_REFINE: NO input FV struct');
end

if ~exist('Nface','var'),
    Nface = 4;
elseif isempty(Nface),
    Nface = 4;
end


switch Nface,
    
case 4,
    FV = mesh_refine_tri4(FV);
case 6,
    FV = mesh_refine_tri6(FV);
otherwise
    FV = mesh_refine_tri4(FV);
end


return
 
