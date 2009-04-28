function [vertices,faces] = mesh_freesurferTRI2matlab(file)

% mesh_freesurferTRI2matlab - Read FreeSurfer tesselation (.tri)
% 
% USEAGE: [vertices,faces] = mesh_freesurferTRI2matlab(file)
% 
% This function will load an ascii file that contains a one
% line specification of the number of vertices followed 
% by rows of vertex points.  It then reads a one line
% specification of the number of faces followed by rows
% of face indices into the vertex rows.  Each vertex row 
% contains a vertex index number and 3 x,y,z coordinates. 
% Each face row contains a face index and three vertex 
% indices.  Vertices in the .tri file are indexed from one
% and those returned are indexed from one.
% 
% See also the mesh_freesurfer2matlab function to load
% the tesselations that are created by the mris_convert 
% function of freesurfer, which have a different text format 
% from those of the BEM .tri files.
% 
% The freesurfer tesselations may contain too many faces
% for efficient computations.  If so, try 'reducepatch'.
%
% The returned matrices can be input to the patch command, like so:
%
%    Hpatch = patch('Vertices',vertices,'Faces',faces,...
%                   'EdgeColor',[.8 .8 .8],'FaceColor',[0.9 0.9 0.9]);
%
% This will plot the mesh as a patch object.  See the patch command
% and matlab help for more information on coloring this object.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:57 $

% Licence:  GNU GPL, no implied or express warranties
% History:  03/02 Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(file,'r');

if isequal(fid,-1),
    S=sprintf('Could not open file: "%s"',file);
    error(S);
else
    
    fprintf('...Reading FreeSurfer Tesselation (.tri)\n');
    
    tic;
    
    % Check for comment on first line of file
    frewind(fid); temp = fscanf(fid,'%s',1); frewind(fid);
    if findstr(temp,'#'), temp = fgetl(fid); end
    
    % Read vertices
    Nvertices = fscanf(fid,'%d',1);
    fprintf('...Reading %d Vertices\n',Nvertices);
    vertices = fscanf(fid,'%f',[4,Nvertices]);
    % remove first row (index) and translate
    vertices = vertices(2:4,:)';
    
    % Read faces
    Nfaces    = fscanf(fid,'%d',1);
    fprintf('...Reading %d Faces\n',Nfaces);
    faces = fscanf(fid,'%d',[4,Nfaces]);
    % remove first row (index) & translate
    faces = faces(2:4,:)';
    
    fclose(fid);
    
    t = toc;
    fprintf('...done (%6.2f sec).\n',t);
    
end

return
