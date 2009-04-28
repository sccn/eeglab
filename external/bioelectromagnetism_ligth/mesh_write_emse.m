function mesh_write_emse(p)

% mesh_write_emse - Save mesh to EMSE (.wfr) file
% 
% USEAGE: mesh_write_emse(p)
% 
% Write a .wfr file, in minor revision 3 format (ascii), 
% for each mesh in p.mesh.data.  If any cells of 
% p.mesh.data.meshtype are empty, these cells will 
% be skipped.
% 
% See the EMSE website at http://www.sourcesignal.com
% for more information on file formats.
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:58 $

% Licence:  GNU GPL, no implied or express warranties
% History:  06/2002 Darren.Weber_at_radiology.ucsf.edu
%                 - created function from mesh_emse2matlab
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nMESH_WRITE_EMSE...\n');

if ~exist('p','var'),
  error('...no input p struct.\n');
elseif isempty(p),
  error('...input p struct is empty.\n');
elseif isempty(p.mesh.data),
  error('...input p struct has no mesh data.\n');
end

[path,name,ext] = fileparts(strcat(p.mesh.path,filesep,p.mesh.file));
file = fullfile(path,[name ext]);

fprintf('...writing EMSE meshes to:\n\t%s\n',fullfile(path,name));

tic;

Meshes = p.mesh.data.meshtype;

for i=1:size(Meshes,2),
  
  if Meshes{i},
    write_emse(file,p.mesh.data.vertices{i},p.mesh.data.faces{i},Meshes{i});
  end
  
end

t=toc; fprintf('...done (%5.2f sec).\n\n',t);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_emse(file,vertex,face,meshtype)

[path,name,ext] = fileparts(file);
ext = '.wfr';

name = strcat(name,'_',meshtype);

file = fullfile(path,[name ext]);


fid = fopen(file,'w','ieee-le');

if(fid == -1),
  fprintf('...could not open file: %s',file);
  return;
else
  
  fprintf('...writing tesselation: %s\n',[name ext]);
  
  % Write prolog
  fprintf(fid,'3\t4000\n');
  fprintf(fid,'3\n');
  
  % Write mesh type
  type = lower(meshtype);
  switch type,
    case 'unknown',     fprintf(fid,  '0\n');
    case 'scalp',       fprintf(fid, '40\n');
    case 'outer skull', fprintf(fid, '80\n');
    case 'outer_skull', fprintf(fid, '80\n');
    case 'inner skull', fprintf(fid,'100\n');
    case 'inner_skull', fprintf(fid,'100\n');
    case 'cortex',      fprintf(fid,'200\n');
    case 'pial',        fprintf(fid,'200\n'); % cortex variant
    case 'white',       fprintf(fid,'200\n'); % cortex variant
    case 'smoothwm',    fprintf(fid,'200\n'); % cortex variant
    otherwise,          fprintf(fid,  '0\n');
      fprintf('\n...WARNING, unknown meshtype!\n\n');
  end
  
  
  % EMSE Voxel Coordinates
  % Voxel coordinates measure location in terms of the voxels inherent in 
  % the given volumetric set. The origin is the bottom (inferior) axial 
  % slice, the posterior row and in the rightmost column. This coordinate 
  % system is right-handed (although, internally, the origin is in the 
  % anterior row, and thus is left-handed; this representation is not 
  % available to the user). The order of the displayed coordinates is 
  % (slice#, row#, column#).
  %
  % EMSE MRI Coordinates
  % MRI coordinates share the same origin as internal voxel coordinates, 
  % but differ from the latter in two ways: first, the coordinates 
  % are measured in millimeters, not voxels. Secondly, the origin is that 
  % of the internal representation; that is, the inferior slice, anterior 
  % row and rightmost column. As mentioned above, this internal representation 
  % is left-handed. To correct for this, the row axis is numbered in the 
  % opposite direction, making the displayed coordinate system right-handed. 
  % The order of the displayed coordinates is (x, y, z).
  
  % Given a point P(x,y,z) in head frame (the activation point on the 
  % cortical mesh) and you want to find the corresponding voxel in the 
  % vmi file.  Symbolically you have P(head) and you want to find P(voxel).
  % 
  % 1.  The registration file contains the matrix HeadToImage,
  %     so P(MRI-mm) = HeadToImage*P(head), where P(MRI-mm) is the 
  %     point in MRI coordinates.
  % 2.  From the voxel size, you can find P(MRI-voxel), which 
  %     is the MRI coordinates expressed in voxels
  % 3.  Use the offset between the MRI coordinate frame and 
  %     the Image coordinate frame to find P(voxel).
  %
  %Demetrios Voreades, Ph.D.
  %Applications Engineer, Source Signal Imaging
  %
  
  
  % Rotate -90 degrees around Z, given that emse coordinates
  % have +X through Nasion and +Y through left ear.
  fprintf('...rotating coordinate axes so +X anterior, +Y left, +Z superior\n');
  vertex = rz(vertex,-90,'degrees');
  
  % Write vertex data
  for v = 1:size(vertex,1),
    fprintf(fid,'v\t%12.8f\t%12.8f\t%12.8f\n',vertex(v,1),vertex(v,2),vertex(v,3));
  end
  
  % matlab vertex indices start at one,
  % not zero, so we subtract one from matlab values
  face = face - 1;
  for t = 1:size(face,1),
    fprintf(fid,'t\t%d\t%d\t%d\t\n',face(t,1),face(t,2),face(t,3));
  end
  
  
  fclose(fid);
  
end

return
