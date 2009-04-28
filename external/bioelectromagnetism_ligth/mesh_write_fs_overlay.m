function mesh_write_fs_overlay(p)

% mesh_write_fs_overlay - Save mesh CData to FreeSurfer overlay file
%
% USEAGE: mesh_write_fs_overlay(p)
%
% Write a binary overlay file (*.w) for each mesh in p.mesh.data.  
% If any cells of p.mesh.data.meshtype are empty, they are skipped.
% 
% This function outputs an *.w overlay file from the
% vertex scalar data contained in p.mesh.data.Cdata - if it 
% contains a timeseries, the output values are taken from 
% the column specified in p.mesh.samplePoint.
% 
% All output is handled by fs_write_wfile.
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

fprintf('\nMESH_WRITE_FS_OVERLAY...\n');

if ~exist('p','var'),
    error('...no input p struct.\n');
elseif isempty(p),
    error('...input p struct is empty.\n');
elseif isempty(p.mesh.data),
    error('...input p struct has no mesh data.\n');
end

[path,name,ext] = fileparts(strcat(p.mesh.path,filesep,p.mesh.file));
file = fullfile(path,[name ext]);

fprintf('...writing FreeSurfer overlays to:\n\t%s\n',path);

tic;

Meshes = p.mesh.data.meshtype;

[path,name,ext] = fileparts(file);

for i=1:size(Meshes,2),
    
    if Meshes{i},
        
        % check for mesh Cdata
        if size(p.mesh.data.Cdata,2) >= i,
            if ~isempty(p.mesh.data.Cdata{i}),
                if size(p.mesh.data.Cdata{i},2) > 1,
                    % Obtain the Cdata at the selected time point
                    w = p.mesh.data.Cdata{i}(:,p.mesh.samplePoint);
                else
                    w = p.mesh.data.Cdata{i};
                end
                
                outputExt = [ext, '.w'];
                outputName = strcat(name,'.',Meshes{i});
                outputFile = fullfile(path,[outputName outputExt]);
                
                fprintf('...writing overlay: %s\n',[outputName outputExt]);
                fs_write_wfile(outputFile,w);
            else
                msg = sprintf('\np.mesh.data.Cdata{%d} is empty\n',i);
                warning(msg);
            end
        else
            msg = sprintf('\nno p.mesh.data.Cdata{%d}\n',i);
            warning(msg);
        end
    end
end

t=toc; fprintf('...done (%5.2f sec).\n\n',t);

return
