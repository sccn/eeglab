function [vertex,face,edge,mesh] = mesh_emse2matlab(file,options)

% mesh_emse2matlab - Convert EMSE mesh (.wfr) to matlab format
%
% USEAGE: [vertices,faces,edges,meshtype] = mesh_emse2matlab(file,[options])
%
% All values returned are in meters.  With the returned structures, 
% create a vertex & face matrix:
% 
%   vertex_matrix = [vertices.x; vertices.y; vertices.z]';
%   face_matrix = [faces.vertex1;faces.vertex2;faces.vertex3]';
%
% These can be input to the patch command:
%
%    Hpatch = patch('Vertices',vertex_matrix,'Faces',face_matrix,...
%                   'EdgeColor',[.6 .6 .6],'FaceColor',[0.9 0.9 0.9]);
%
% See the patch and light commands to colour this object.
% 
% 'options' ... a cell array of strings.  By default it contains
% options = {'vertex','face','edge'}.  By default, this routine
% reads all available data from the emse file.  If 'options' is given
% and it doesn't contain one of these strings, that data will not be
% read or returned.
%
% meshtype is: 'unknown','scalp','outer skull','inner skull', or 'cortex'
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:57 $

% Licence:  GNU GPL, no implied or express warranties
% History:  10/98 Abbas Kouzani (egazk@flinders.edu.au)
%           09/01 Darren.Weber_at_radiology.ucsf.edu
%                 - created function, rather than script
%                 - added functionality to handle different
%                   minor revisions.
%           03/02 Darren.Weber_at_radiology.ucsf.edu
%                 - optimised fscanf processes for matlab, the
%                   function now operates in less than 1/2 time,
%                   especially for revision 3 data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('options','var'),
    options = {'vertex','face','edge'};
end

[path,name,ext] = fileparts(file);
file = fullfile(path,[name ext]);

[fid,msg] = fopen(file,'r');
if ~isempty(msg), error(msg); end

tic;

% Read prolog
version   =fscanf(fid,'%f',1);
file_type =fscanf(fid,'%f',1);
minor_rev =fscanf(fid,'%f',1);

fprintf('...WFR Version = %d, Minor_Revision = %d, File-Type = %d\n',...
    version,minor_rev,file_type);

if ~(file_type==8000 | file_type==4000)
    S=sprintf('Could not convert WFR file type: %d',file_type);
    error(S);
end

% Read header (format depends on minor revision)
if(minor_rev==3)
    type      =fscanf(fid,'%f',1);
else
	radius    =fscanf(fid,'%f',1);
	vertex_no =fscanf(fid,'%d',1);
	face_no   =fscanf(fid,'%d',1);
	edge_no   =fscanf(fid,'%d',1);
    
    fprintf('...Radius = %f meters\n',radius);
    %...Vertices = %d\n...Patches = %d\n...Edges = %d\n',...
    %    vertex_no,face_no,edge_no);

    if(minor_rev==1), type = 0;
    else,             type = fscanf(fid,'%f',1);
    end
end

if(minor_rev==1)
    mesh = 'unknown';
else
	switch type
        case 0,             mesh = 'unknown';
        case { 64,  40},    mesh = 'scalp';
        case {128,  80},    mesh = 'outer skull';
        case {256, 100},    mesh = 'inner skull';
        case {512, 200},    mesh = 'cortex';
        otherwise,          mesh = 'unknown';
	end
end

fprintf('...Mesh type: %s\n',mesh);

% Read data (format depends on minor revision)
vertex = struct([]);
face = struct([]);
edge = struct([]);

if(minor_rev==3)
    
    % Read the whole file
    fprintf('...Reading Minor Revision 3 Data File');
    Tmp = fscanf(fid,'%s%f%f%f',[4,inf]);
    fprintf('...done\n');
    
    % Vertices
    if(max(strcmp('vertex',options)) > 0),
        fprintf('...Creating Vertices Struct');
        % first get the numeric code for 'v'
        vcode = Tmp(1,1);
        % now find all columns of Tmp with vcode in 1st row
        vindex = find(Tmp(1,:) == vcode);
        
        tmp = Tmp(2:4,vindex);
        tmp = num2cell(tmp);
        vertex = struct(...
            'index', [1:length(vindex)],...
            'x',     tmp( 1,:),...
            'y',     tmp( 2,:),...
            'z',     tmp( 3,:));
        clear tmp;
        fprintf('...done\n');
    else
        fprintf('...Skipping Vertices Struct\n');
    end
    
    % Faces
    if(max(strcmp('face',options)) > 0),
        fprintf('...Creating Faces Struct');
        % first get the numeric code for 't'
        tcode = Tmp(1,length(vindex)+1);
        % now find all columns of Tmp with tcode in 1st row
        tindex = find(Tmp(1,:) == tcode);
        
        % matlab vertex indices start at one,
        % not zero, so we add one to these emse values
        tmp = Tmp(2:4,tindex) + 1;
        
        tmp = num2cell(tmp);
        face = struct(...
            'index',   [1:length(tindex)],...
            'vertex1', tmp( 1,:),...
            'vertex2', tmp( 2,:),...
            'vertex3', tmp( 3,:));
        clear tmp;
        fprintf('...done\n');
    else
        fprintf('...Skipping Faces Struct\n');
    end
    
    % Edges
    fprintf('...No edges for minor revision 3\n');
    clear Tmp;
    
elseif(minor_rev~=4)
    % minor revision 1 & 2 format
    disp('...Reading Minor Revision 1 or 2 Data');
    if(max(strcmp('vertex',options)) > 0),
        fprintf('...Reading %d Vertices',vertex_no);
        
        tmp = zeros(13,vertex_no);
        tmp = fscanf(fid,'%d%x%d%d%g%g%g%g%g%g%g%g%g',[13,vertex_no]);
        
        tmp(1,:) = tmp(1,:) + 1;
        tmp      = num2cell(tmp);
        vertex   = struct(...
            'index',        tmp( 1,:),...
            'address',      tmp( 2,:),...
            'channel_index',tmp( 3,:),...
            'x',            tmp( 5,:),...
            'y',            tmp( 6,:),...
            'z',            tmp( 7,:),...
            'xnormal',      tmp( 9,:),...
            'ynormal',      tmp(10,:),...
            'znormal',      tmp(11,:),...
            'potential',    tmp(12,:),...
            'curvature',    tmp(13,:));
        clear tmp;
        fprintf('...done\n');
    else
        fprintf('...Skipping %d Vertices\n',vertex_no);
        tmp = zeros(13,vertex_no);
        tmp = fscanf(fid,'%d%x%d%d%g%g%g%g%g%g%g%g%g',[13,vertex_no]);
        clear tmp;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(max(strcmp('face',options)) > 0),
        fprintf('...Reading %d Faces',face_no);
        
        tmp = zeros(18,face_no);
        tmp = fscanf(fid,'%d%x%g%g%g%g%g%g%g%g%g%g%x%x%x%x%x%x',[18,face_no]);
        tmp(1,:) = tmp(1,:) + 1;
        tmp      = num2cell(tmp);
        face     = struct(...
            'index',        tmp( 1,:),...
            'address',      tmp( 2,:),...
            'solid_angle',  tmp( 3,:),...
            'magnitude',    tmp( 4,:),...
            'potential',    tmp( 5,:),...
            'area',         tmp( 6,:),...
            'center_x',     tmp( 7,:),...
            'center_y',     tmp( 8,:),...
            'center_z',     tmp( 9,:),...
            'normal_x',     tmp(10,:),...
            'normal_y',     tmp(11,:),...
            'normal_z',     tmp(12,:),...
            'vertex1',      tmp(13,:),...
            'vertex2',      tmp(14,:),...
            'vertex3',      tmp(15,:),...
            'edge1',        tmp(16,:),...
            'edge2',        tmp(17,:),...
            'edge3',        tmp(18,:));
        clear tmp;
        fprintf('...done\n');
        
        % In minor rev4, the face vertex and edges
        % refer to the vertex.address field, so this
        % is corrected here.  Not sure how to avoid 'for'
        fprintf('...Converting Face vertices from address to index (this takes a while)');
        for i=1:face_no,
            face(i).vertex1 = find([vertex.address] == face(i).vertex1);
            face(i).vertex2 = find([vertex.address] == face(i).vertex2);
            face(i).vertex3 = find([vertex.address] == face(i).vertex3);
        end
        fprintf('...done\n');
        if(max(strcmp('edge',options)) > 0),
            fprintf('...Converting Face edges from address to index (this takes a while)');
            for i=1:face_no,
                face(i).edge1 = find([vertex.address] == face(i).edge1);
                face(i).edge2 = find([vertex.address] == face(i).edge2);
                face(i).edge3 = find([vertex.address] == face(i).edge3);
            end
            fprintf('...done\n');
        end
    else
        disp('...Skipping Faces');
        tmp = zeros(18,face_no);
        tmp = fscanf(fid,'%d%x%g%g%g%g%g%g%g%g%g%g%x%x%x%x%x%x',[18,face_no]);
        clear tmp;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(max(strcmp('edge',options)) > 0),
        fprintf('...Reading %d Edges',edge_no);
        
        tmp = zeros(4,edge_no);
        tmp = fscanf(fid,'%d%x%x%x',[4,edge_no]);
        
        tmp(1,:) = tmp(1,:) + 1;
        tmp      = num2cell(tmp);
        edge     = struct(...
            'index',        tmp(1,:),...
            'address',      tmp(2,:),...
            'vertex1',      tmp(3,:),...
            'vertex2',      tmp(4,:));
        clear tmp;
        fprintf('...done\n');
        fprintf('...Converting Edge vertices from address to index (this takes a while)');
        for i=1:edge_no,
            edge(i).vertex1 = find([vertex.address] == edge(i).vertex1);
            edge(i).vertex2 = find([vertex.address] == edge(i).vertex2);
        end
        fprintf('...done\n');
    else
        disp('...Skipping Edges');
        %for i=1:edge_no,
        %    tmp=fscanf(fid,'%*d',1);
        %    tmp=fscanf(fid,'%*x',3);
        %end
    end
else
    % minor revision 4 format
    disp('...Reading Minor Revision 4 Data');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(max(strcmp('vertex',options)) > 0),
        fprintf('...Reading %d Vertices',vertex_no);
        
        index = meshgrid(1:1:vertex_no,1);
        index = num2cell(index);
        
        tmp = zeros(11,vertex_no);
        tmp = fscanf(fid,'%d%d%g%g%g%d%g%g%g%g%g',[11,vertex_no]);
        
        tmp      = num2cell(tmp);
        vertex   = struct(...
            'index',        index,...
            'channel_index',tmp( 1,:),...
            'x',            tmp( 3,:),...
            'y',            tmp( 4,:),...
            'z',            tmp( 5,:),...
            'xnormal',      tmp( 7,:),...
            'ynormal',      tmp( 8,:),...
            'znormal',      tmp( 9,:),...
            'potential',    tmp(10,:),...
            'curvature',    tmp(11,:));
        clear tmp index;
        fprintf('...done\n');
    else
        disp('...Skipping Vertices');
        tmp = zeros(11,vertex_no);
        tmp = fscanf(fid,'%d%d%g%g%g%d%g%g%g%g%g',[11,vertex_no]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(max(strcmp('face',options)) > 0),
        fprintf('...Reading %d Faces',face_no);
        
        index = meshgrid(1:1:face_no,1);
        index = num2cell(index);
        
        tmp = zeros(16,face_no);
        tmp = fscanf(fid,'%g%g%g%g%g%g%g%g%g%g%d%d%d%d%d%d',[16,face_no]);
        
        tmp(11:16,:) = tmp(11:16,:) + 1;
        
        tmp      = num2cell(tmp);
        face     = struct(...
            'index',        index,...
            'solid_angle',  tmp( 1,:),...
            'magnitude',    tmp( 2,:),...
            'potential',    tmp( 3,:),...
            'area',         tmp( 4,:),...
            'center_x',     tmp( 5,:),...
            'center_y',     tmp( 6,:),...
            'center_z',     tmp( 7,:),...
            'normal_x',     tmp( 8,:),...
            'normal_y',     tmp( 9,:),...
            'normal_z',     tmp(10,:),...
            'vertex1',      tmp(11,:),...
            'vertex2',      tmp(12,:),...
            'vertex3',      tmp(13,:),...
            'edge1',        tmp(14,:),...
            'edge2',        tmp(15,:),...
            'edge3',        tmp(16,:));
        clear tmp index;
        fprintf('...done\n');
    else
        disp('...Skipping Faces');
        tmp = zeros(16,face_no);
        tmp = fscanf(fid,'%g%g%g%g%g%g%g%g%g%g%d%d%d%d%d%d',[16,face_no]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(max(strcmp('edge',options)) > 0),
        fprintf('...Reading %d Edges',edge_no);
        
        index = meshgrid(1:1:edge_no,1);
        address = num2cell(index * 0);
        index = num2cell(index);
        
        tmp = zeros(2,edge_no);
        tmp = fscanf(fid,'%d%d',[2,edge_no]);
        tmp = tmp + 1;
        tmp = num2cell(tmp);
        
        edge = struct(...
            'index',   index,...
            'address', address,...
            'vertex1', tmp(1,:),...
            'vertex2', tmp(2,:));
        clear tmp index address;
        fprintf('...done\n');
    else
        fprintf('...Skipping Edges\n');
        %tmp=fscanf(fid,'%*d',2*edge_no);
    end
end

fclose(fid);
t=toc;
fprintf('...done (%5.2f sec)\n',t);

return
