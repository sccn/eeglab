function mesh_write_brainstorm(p)

% mesh_write_brainstorm - Save eeg toolbox meshes to brainstorm format
%
% USEAGE: mesh_write_brainstorm(p)
% 
% This routine saves each mesh in p.mesh.data into the correct 
% structure format required in brainstorm.
% 
% The workspace variables created are saved, in the brainstorm 
% file format, into p.mesh.file with a '.mat' extension.
% 
% The variables from this file can be loaded using the matlab 'load' 
% command. To plot the default cortex, skull & scalp BrainStorm data:
%
% Hpatch1 = patch('Vertices',Vertices{1}','Faces',Faces{1},...
%                 'EdgeColor',[.6 .6 .6],'FaceColor',[0.9 0.9 0.9]);
% Hpatch2 = patch('Vertices',Vertices{2}','Faces',Faces{2},...
%                 'EdgeColor',[.6 .6 .6],'FaceColor',[0.9 0.9 0.9]);
% Hpatch3 = patch('Vertices',Vertices{3}','Faces',Faces{3},...
%                 'EdgeColor',[.6 .6 .6],'FaceColor',[0.9 0.9 0.9]);
%
% See the BrainStorm website at http://neuroimage.usc.edu/brainstorm/
% for more information about the BrainStorm toolbox and the format
% and content of the subjecttess variables.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:58 $

% Licence:  GNU GPL, no implied or express warranties
% History:  04/2002, Darren.Weber_at_radiology.ucsf.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nMESH_WRITE_BRAINSTORM...\n');

if ~exist('p','var'),
    error('...no input p struct.\n');
else
    if isempty(p.mesh.data),
        error('...input p.mesh.data is empty.\n');
    end
end

[path,name,ext] = fileparts(strcat(p.mesh.path,filesep,p.mesh.file));
ext = '.mat';
brainstormfile = fullfile(path,[name ext]);

if ~exist('brainstormfile','var'),
    error('...no input brainstormfile.\n');
end

tic;

Comment  = p.mesh.data.meshtype;
Faces    = cell(size(Comment));
Vertices = cell(size(Comment));

for i=1:size(Comment,2),
        fprintf('...converting tesselation: %s\n',Comment{i});
        Vertices{i} = p.mesh.data.vertices{i}';  % transpose vertices
        Faces{i}    = p.mesh.data.faces{i};
end

fprintf('...saving BrainStorm data to:\n\t%s\n',brainstormfile);
save(brainstormfile, 'Comment', 'Faces', 'Vertices');

t=toc; fprintf('...done (%5.2f sec).\n\n',t);

return
