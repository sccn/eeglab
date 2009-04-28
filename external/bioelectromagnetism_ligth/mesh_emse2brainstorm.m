function [Comment,Faces,Vertices] = mesh_emse2brainstorm(file,meshes)

% mesh_emse2brainstorm - Convert EMSE meshes (.wfr) to brainstorm format
%
% USEAGE: [Comment,Faces,Vertices] = mesh_emse2brainstorm(fileprefix,meshes)
%
% This routine calls the mesh_emse_wfr2matlab routine for
% each layer of the subject tesselation required in brainstorm
% (ie, cortex, inner skull, scalp).  It returns the variables
% required by BrainStorm for source modelling.  The EMSE coordinate
% framework is different from that of Neuroscan and Brainstorm, so
% all vertices and faces are rearranged to suit.
%
% fileprefix should be a string containing the path and subject
% code for the emse meshes to load and convert.  For example:
% 'c:\data\c01_'.  If the cell array 'meshes' is given, it should 
% contain the names of the meshes to load.  Only those meshes will 
% be loaded. For example: meshes = {'cortex','innerskull','scalp'} 
% is the default.  This routine then loads the following EMSE meshes:
%
% c:\data\c01_cortex.wfr
% c:\data\c01_innerskull.wfr
% c:\data\c01_scalp.wfr
%
% The workspace variables created are saved, in the brainstorm 
% file format, into the file called 'c:\data\c01_subjecttess.mat'.
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
% and content of the subjecttess variables.  Also, see the EMSE/MRVU
% website at http://www.sourcesignal.com.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:57 $

% Licence:  GNU GPL, no implied or express warranties
% History:  01/02 Darren.Weber_at_radiology.ucsf.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check for existence of saved file
brainstormfile = strcat(file,'subjecttess');

fid = fopen(strcat(brainstormfile,'.mat'),'r');

if (fid ~= -1),
    fclose(fid);
    fprintf('\nLoading BrainStorm data from %s',brainstormfile);
    fprintf('\nTo reload data from EMSE, first delete/rename BrainStorm data.\n\n');
    load(brainstormfile);
	return;
else
	
	if ~exist('meshes','var'),
        meshes = {'cortex','innerskull','scalp'};
	end
	
	Comment = meshes;
	Faces = cell(size(meshes));
	Vertices = cell(size(meshes));
	
	% Load EMSE data
	options = {'vertex','patch'}; % only load vertex & patch data
	for i = 1:max(size(meshes)),
        fprintf('\nLoading EMSE mesh: %s\n',meshes{i});
        inputfile = strcat(file,meshes{i},'.wfr');
        [vertices,patches] = mesh_emse2matlab(inputfile,options);
        % EMSE vertices are rearranged here, swapping x & y to
        % match the coordinate frame of the NeuroScan 3Dspace electrodes.
        Vertices{i} = [ vertices.y; vertices.x; vertices.z ]; % 3 x n
        % The faces remain the same order as they refer only to the
        % vertices as a set of 3D coordinates.
        Faces{i} = [ patches.vertex1; patches.vertex2; patches.vertex3 ]'; % m x 3
	end
	
	fprintf('\nSaving BrainStorm data to %s',brainstormfile);
	save(brainstormfile, 'Comment', 'Faces', 'Vertices');
end
