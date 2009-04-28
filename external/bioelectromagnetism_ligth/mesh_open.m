function [p] = mesh_open(p)

% mesh_open - calls functions to read a triangulation file
% 
% Usage: [p] = mesh_open(p)
% 
% p is a parameter structure (see eeg_toolbox_defaults for
% more details). In this function, it should contain at least
% the following string fields:
%
%       p.mesh.path - the directory location of the file to load
%       p.mesh.file - the name of the file to load
%       p.mesh.type - the file format:
%
%      'emse'
%      'brainstorm'
%      'ascii' or 'freesurfer_ascii' for *.asc/.tri files
%      'freesurfer_surf' for freesurfer binary surface
%      'freesurfer_curv' for freesurfer binary curvature
%      'freesurfer_overlay' for freesurfer binary overlay (.w)  
% 
% The last two file types are associated with the current cortex,
% if it contains the same number of vertices (eg, load 'freesurfer_surf'
% first, then 'freesurfer_curv' or 'freesurfer_overlay').
% 
% The file formats supported here are described in more detail 
% at their respective websites:
% 
%       FreeSurfer: http://surfer.nmr.mgh.harvard.edu
%       EMSE:       http://www.sourcesignal.com
%       BrainStorm: http://neuroimage.usc.edu/brainstorm
% 
% The return structure creates or updates p.mesh.data, which 
% contains cell arrays:
%
%       p.mesh.data.meshtype     type of surface, strings
%       p.mesh.data.vertices     Mx3 (x,y,z) vertices
%       p.mesh.data.faces        Mx3 vertex indices
%       p.mesh.data.Cdata        scalar overlay, M vert x N overlays
%
% For example, having loaded a scalp and an inner skull mesh, 
% p.mesh.data.meshtype could be:
% 
%       p.mesh.data.meshtype{1} = 'scalp'
%       p.mesh.data.meshtype{2} = 'inner skull'
%
% To plot the data returned, see mesh_plot or try
% 
% Hpatch = patch('Vertices',p.mesh.data.vertices{1}',...
%                'Faces',p.mesh.data.faces{1},...
%                'Property','PropertyValue',...);
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:57 $

% Licence:  GNU GPL, no express or implied warranties
% History:  02/2002 Darren.Weber_at_radiology.ucsf.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nMESH_OPEN...\n'); tic;

if ~exist('p','var'),
 [p] = eeg_toolbox_defaults;
  fprintf('...creating default p structure.\n');
end

[path,name,ext] = fileparts(strcat(p.mesh.path,filesep,p.mesh.file));
file = fullfile(path,[name ext]);

if exist(file) ~= 2,
  msg = sprintf('\n\n...file does not exist:\n\t%s\n',file);
  error(msg);
end

type = lower(p.mesh.type);

switch type,
  
  case 'emse',
    
    fprintf('...loading EMSE mesh from:\n\t%s\n',file);
    
    % Get EMSE .wfr data (in meters)
    options = {'vertex','face'};
    [vertices,faces,edges,meshtype] = mesh_emse2matlab(file,options);
    
    vertex_matrix = [vertices.x; vertices.y; vertices.z]';
    face_matrix = [faces.vertex1;faces.vertex2;faces.vertex3]';
    
    % Rotate 90 degrees around Z for EMSE data
    vertex_matrix = Rz(vertex_matrix,90,'degrees');
    
    % Is this a new or replacement mesh?
    [p.mesh.current,meshExists] = mesh_check(p,meshtype);
    
    p.mesh.data.meshtype{p.mesh.current}    = meshtype;
    p.mesh.data.vertices{p.mesh.current}    = vertex_matrix;
    p.mesh.data.faces{p.mesh.current}       = face_matrix;
    p.mesh.data.lapmat{p.mesh.current}      = [];
    p.mesh.data.lapint{p.mesh.current}      = [];
    p.mesh.data.timeseries{p.mesh.current}  = [];
    p.mesh.data.Cdata{p.mesh.current}       = [];
    
  case 'brainstorm',
    
    fprintf('...loading BrainStorm data from:\n\t%s\n',file);
    load(file);
    
    p.mesh.data = [];
    
    for i=1:size(Comment,2),
      if isempty(Comment{i}), continue; end
      fprintf('...converting tesselation: %s\n',Comment{i});
      p.mesh.data.meshtype{i}   = Comment{i};
      p.mesh.data.vertices{i}   = Vertices{i}';  % transpose Vertices
      p.mesh.data.faces{i}      = Faces{i};
      p.mesh.data.lapmat{i}     = [];
      p.mesh.data.lapint{i}     = [];
      p.mesh.data.timeseries{i} = [];
      p.mesh.data.Cdata{i}      = [];
      p.mesh.current = i;
    end
    
  case {'ascii','freesurfer_ascii'},
    
    fprintf('...loading ASCII or FreeSurfer data from:\n\t%s\n',file);
    
    % Get Freesurfer data
    if findstr(file,'.tri'),
      [vertex_matrix,face_matrix] = freesurfer_read_tri(file);
    else
      [vertex_matrix,face_matrix] = freesurfer_read_ascii(file);
    end
    
    fprintf('...converting surface coordinates (mm to meters)\n');
    vertex_matrix = vertex_matrix ./ 1000;
    
    meshtype = freesurfer_meshtype(file); % see function below
    
    [p.mesh.current,meshExists] = mesh_check(p,meshtype);
    
    p.mesh.data.meshtype{p.mesh.current}   = meshtype;
    p.mesh.data.vertices{p.mesh.current}   = vertex_matrix;
    p.mesh.data.faces{p.mesh.current}      = face_matrix;
    p.mesh.data.lapmat{p.mesh.current}     = [];
    p.mesh.data.lapint{p.mesh.current}     = [];
    p.mesh.data.timeseries{p.mesh.current} = [];
    p.mesh.data.Cdata{p.mesh.current}      = [];
    
  case 'freesurfer_surf',
    
    fprintf('...loading FreeSurfer binary surface from:\n\t%s\n',file);
    
    meshtype = freesurfer_meshtype(file); % see function below
    
    [p.mesh.current,meshExists] = mesh_check(p,meshtype);
    
    [vertex_matrix, face_matrix] = freesurfer_read_surf(file);
    
    fprintf('...converting surface coordinates (mm to meters)\n');
    vertex_matrix = vertex_matrix ./ 1000;
    
    p.mesh.data.meshtype{p.mesh.current}   = meshtype;
    p.mesh.data.vertices{p.mesh.current}   = vertex_matrix;
    p.mesh.data.faces{p.mesh.current}      = face_matrix;
    p.mesh.data.lapmat{p.mesh.current}     = [];
    p.mesh.data.lapint{p.mesh.current}     = [];
    p.mesh.data.timeseries{p.mesh.current} = [];
    p.mesh.data.Cdata{p.mesh.current}      = [];
    
  case 'freesurfer_curv',
    
    fprintf('...loading FreeSurfer binary curvature from:\n\t%s\n',file);
    
    [curv, Nfaces] = freesurfer_read_curv(file);
    
    allocated = 0;
    for meshN = 1:length(p.mesh.data.meshtype),
      
      type = p.mesh.data.meshtype{meshN};
      
      if size(p.mesh.data.vertices{meshN},1) == size(curv,1),
        fprintf('...allocating overlay to Cdata for ''%s'' surface.\n',type);
        p.mesh.data.Cdata{meshN} = curv;
        p.mesh.current = meshN;
        allocated = 1;
      end
      
    end
    if allocated < 1,
      fprintf('...failed to allocate curvature to any surface, incompatible Nfaces!\n');
    end
    
  case 'freesurfer_overlay',
    
    fprintf('...loading FreeSurfer binary overlay from:\n\t%s\n',file);
    
    [w,vert] = freesurfer_read_wfile(file);
    
    Nvert = length(vert); clear vert;
    
    allocated = 0;
    for meshN = 1:length(p.mesh.data.meshtype),
      
      type = p.mesh.data.meshtype{meshN};
      
      if size(p.mesh.data.vertices{meshN},1) == Nvert,
        fprintf('...allocating overlay to Cdata for ''%s'' surface.\n',type);
        p.mesh.data.Cdata{meshN} = w;
        p.mesh.current = meshN;
        allocated = 1;
      end
      
    end
    if allocated < 1,
      fprintf('...failed to allocate overlay to any surface, incompatible vertices!\n');
    end
    
  otherwise,
    fprintf('...mesh format: %s\n', p.mesh.type);
    fprintf('...sorry, cannot load this data format at present.\n');
    return;
end

t=toc; fprintf('...done (%5.2f sec).\n\n',t);

return
