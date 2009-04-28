function [p] = mesh_write(p)

% mesh_write - a switch yard for several mesh writing functions
% 
% Usage: [p] = mesh_write(p)
% 
% p is a parameter structure (see eeg_toolbox_defaults for
% more details). In this function, it should contain mesh data
% in p.mesh.data and the following string fields:
% 
%       p.mesh.path - the directory location of the file to load
%       p.mesh.file - the name of the file to load
%       p.mesh.type - the file formats (case insensitive) are:
% 
%                    'emse'
%                    'brainstorm'
%                    'ascii' or 'freesurfer' for *.asc file
%                    'fs_surf' for freesurfer binary surface
%                    'fs_curv' for freesurfer binary curvature
%                    'fs_overlay' for freesurfer binary overlay (.w)
% 
% If you have a file format you would like supported, please email
% the author with an example and description of the format.  The 
% file formats supported here are described in more detail at 
% their respective websites:
% 
%       EMSE:       http://www.sourcesignal.com
%       BrainStorm: http://neuroimage.usc.edu/brainstorm/
%       FreeSurfer: http://surfer.nmr.mgh.harvard.edu/
% 
% See also, mesh_open, mesh_write_emse, 
%           mesh_write_brainstorm, mesh_write_freesurfer
%           mesh_write_fs_surf, mesh_write_fs_curv, 
%           mesh_write_fs_overlay
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:58 $

% Licence:  GNU GPL, no express or implied warranties
% History:  02/2002 Darren.Weber_at_radiology.ucsf.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~exist('p','var'),
    msg = sprintf('\nMESH_WRITE...no input p struct.\n');
    error(msg);
elseif isempty(p),
    msg = sprintf('\nMESH_WRITE...input p struct is empty.\n');
    error(msg);
elseif isempty(p.mesh.data),
    msg = sprintf('\nMESH_WRITE...input p struct has no mesh data.\n');
    error(msg);
end

type = lower(p.mesh.type);

switch type,
    
case 'emse',
    
    mesh_write_emse(p);
    
case 'brainstorm',
    
    mesh_write_brainstorm(p);
    
case {'ascii','freesurfer'}, % for *.asc file
    
    mesh_write_freesurfer(p);
    
case 'fs_surf', % for freesurfer binary surface
    
    mesh_write_fs_surf(p);
    
case 'fs_curv', % for freesurfer binary curvature
    
    mesh_write_fs_curv(p);
    
case 'fs_overlay', % for freesurfer binary overlay (.w)
    
    mesh_write_fs_overlay(p);
    
otherwise,
    fprintf('...mesh format: %s\n', p.mesh.type);
    fprintf('...sorry, cannot write this data format at present.\n');
    return;
end

return
