function mesh_write_fs_surf(p)

% mesh_write_fs_surf - Save mesh to FreeSurfer (.asc) file
%
% USEAGE: mesh_write_fs_surf(p)
% 
% Write a binary surface file for each mesh in p.mesh.data.  If any 
% cells of p.mesh.data.meshtype are empty, they will be skipped.
% 
% All output is handled by fs_write_surf.  This is just a wrapper
% to ouput all p.mesh.data, converting from meters into mm.
% 
% See also fs_write_surf, fs_read_surf, and
% the FreeSurfer website at http://surfer.nmr.mgh.harvard.edu/.
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:58 $

% Licence:  GNU GPL, no implied or express warranties
% History:  10/2002, Darren.Weber_at_radiology.ucsf.edu
%           01/2003, Darren.Weber_at_radiology.ucsf.edu
%                    now using fs_write_surf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nMESH_WRITE_FS_SURF...\n');

if ~exist('p','var'),
    error('...no input p struct.\n');
elseif isempty(p),
    error('...input p struct is empty.\n');
elseif isempty(p.mesh.data),
    error('...input p struct has no mesh data.\n');
end

[path,name,ext] = fileparts(strcat(p.mesh.path,filesep,p.mesh.file));
file = fullfile(path,[name ext]);

fprintf('...writing FreeSurfer binary surfaces to:\n\t%s\n',path);

tic;

Meshes = p.mesh.data.meshtype;

[path,name,ext] = fileparts(file);

for i=1:size(Meshes,2),
    
    if Meshes{i},
        
        % Convert vertices from meters to mm
        fprintf('...converting vertices from meters to mm.\n');
        vertices = p.mesh.data.vertices{i} .* 1000;
        faces    = p.mesh.data.faces{i};
        
        outputName = strcat(name,'.',Meshes{i});
        outputFile = fullfile(path,outputName);
        
        fprintf('...writing tesselation: %s\n',outputName);
        fs_write_surf(outputFile,vertices,faces);
        
    end
end

t=toc; fprintf('...done (%5.2f sec).\n\n',t);

return
