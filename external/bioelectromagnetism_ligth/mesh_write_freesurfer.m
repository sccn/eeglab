function mesh_write_freesurfer(p)

% mesh_write_freesurfer - Save mesh to FreeSurfer (.asc) file
%
% USEAGE: mesh_write_freesurfer(p)
%
% Write a binary surface file for each mesh in p.mesh.data.  If any 
% cells of p.mesh.data.meshtype are empty, they will be skipped.
% 
% This function will also output an *.w overlay file from the
% vertex scalar data contained in p.mesh.data.Cdata - if it 
% contains a timeseries, the output values will be taken from 
% the column specified in p.mesh.samplePoint.
% 
% All output is handled by fs_write_surf and fs_write_wfile.
% 
% See the FreeSurfer website at http://surfer.nmr.mgh.harvard.edu/
% for more information.
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:58 $

% Licence:  GNU GPL, no implied or express warranties
% History:  10/2002, Darren.Weber_at_radiology.ucsf.edu
%           01/2003, Darren.Weber_at_radiology.ucsf.edu
%                    now using fs_write_surf and fs_write_wfile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nMESH_WRITE_FREESURFER...\n');

if ~exist('p','var'),
    error('...no input p struct.\n');
elseif isempty(p),
    error('...input p struct is empty.\n');
elseif isempty(p.mesh.data),
    error('...input p struct has no mesh data.\n');
end

[path,name,ext] = fileparts(strcat(p.mesh.path,filesep,p.mesh.file));
file = fullfile(path,[name ext]);

fprintf('...writing FreeSurfer/ASCII meshes to:\n\t%s\n',path);

tic;

Meshes = p.mesh.data.meshtype;

for i=1:size(Meshes,2),
    
    if Meshes{i},
        
        vertices = p.mesh.data.vertices{i};
        faces    = p.mesh.data.faces{i};
        
        [path,name,ext] = fileparts(file);
        name = strcat(name,'.',Meshes{i});
        file = fullfile(path,[name ext]);
        
        % Convert vertices from meters to mm
        vertices = vertices .* 1000;
        fs_write_surf(file,vertices,faces);
        
        % check for mesh Cdata
        if size(p.mesh.data.Cdata,2) >= i,
            if p.mesh.data.Cdata{i},
                if size(p.mesh.data.Cdata{i},2) > 1,
                    % Obtain the Cdata at the selected time point
                    w = p.mesh.data.Cdata{i}(:,p.mesh.samplePoint);
                else
                    w = p.mesh.data.Cdata{i};
                end
            end
        end
        
        % Output associated overlay file (*.w)
        file = fullfile(path,[name '.w']);
        fs_write_wfile(file,w);
        
        
        % ascii output format replaced by binary option 01/2003
        % write_freesurfer(file,vertices,faces,Meshes{i});
        
    end
end

t=toc; fprintf('...done (%5.2f sec).\n\n',t);

return


% following asc output is obsolete, given fs_write_surf, 01/2003


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_freesurfer(file,vertex,face,meshtype)

    [path,name,ext] = fileparts(file);
    ext = [ext,'.asc'];
    
    name = strcat(name,'.',meshtype);
    file = fullfile(path,[name ext]);
    
    fid = fopen(file,'w','ieee-le');
    
    if(fid == -1),
        msg = sprintf('...could not open file:\n\t%s\n',file);
        error(msg);
    else
        
        fprintf('...writing tesselation: %s\n',[name ext]);
        
        %-------------------------------------------------
        % Output header
        
        Nvertex = size(vertex,1);
        Nfaces  = size(face,1);
        
        fprintf(fid,'%d %d\n',Nvertex,Nfaces);
        
        %-------------------------------------------------
        % Output vertices
        
        % Convert vertices from meters to mm
        vertex(:,1) = vertex(:,1) .* 1000;
        vertex(:,2) = vertex(:,2) .* 1000;
        vertex(:,3) = vertex(:,3) .* 1000;
        
        if size(vertex,2) < 4,
            vertex(:,4) = zeros(size(vertex,1),1);
        end
        
        vertex(:,4) = zeros(size(vertex,1),1);
        
        % Write vertex matrix
        for v = 1:Nvertex,
            fprintf(fid,'%f %f %f %g\n',vertex(v,1),vertex(v,2),vertex(v,3),vertex(v,4));
            %fprintf(fid,'%12.6f %12.6f %12.6f %+20.8g\n',vertex(v,1),vertex(v,2),vertex(v,3),vertex(v,4));
            %fprintf(fid,'%12.6f %12.6f %12.6f %12.6f\n',vertex(v,1),vertex(v,2),vertex(v,3),vertex(v,4));
        end
        
        %-------------------------------------------------
        % Output faces
        
        % Check for last column of face
        if size(face,2) < 4,
            face(:,4) = ones(size(face,1),1); % -1 below
        end
        % matlab vertex indices start at one, so
        % subtract 1 because FreeSurfer vertices start at zero
        face = face - 1;
        % Write face matrix
        for t = 1:Nfaces,
            fprintf(fid,'%d %d %d %d\n',face(t,1),face(t,2),face(t,3),face(t,4));
        end
        
        
        fclose(fid);
        
	end
    
return
