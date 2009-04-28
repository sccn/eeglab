function [FV] = mesh_bet2matlab(fileprefix)

% mesh_bet2matlab - Read FSL BET tesselation (.coo/.dat)
% 
% USEAGE: [FV] = mesh_bet2matlab(fileprefix)
% 
% This function will load ascii files created by the FSL
% BET function, see http://www.fmrib.ox.ac.uk/fsl
% 
% When used with the -x option, the FSL BET function can 
% output the tesselation of the brain surface (ie, only 
% a fairly smooth, non-detailed, surface), eg:
% 
% bet T1 T1_brain -x
% 
% The brain tesselation is output into two files, eg:
% 
%   T1_brain.coo
%   T1_brain.dat
% 
% The T1_brain.coo file contains the co-ordinates of the
% vertices, and the .dat file contains the details of 
% the links between vertices (each line corresponds to 
% each vertex in the .coo file, and the 1's correspond 
% to which other vertices they are connected to).
% 
% The returned FV struct contains the 'vertices' and
% 'faces' matrices, which can be input to the patch 
% command, eg:
% 
% Hpatch = patch('Vertices',FV.vertices,'Faces',FV.faces,...
%                'EdgeColor',[.8 .8 .8],...
%                'FaceColor',[0.9 0.9 0.9]);
% 
% This will plot the mesh as a patch object.  See the patch 
% command and matlab help for more information on coloring 
% this object.
% 
% See also: mesh_freesurfer2matlab
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:57 $

% Licence:  GNU GPL, no implied or express warranties
% History:  03/02 Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('fileprefix','var'),
    error('No input fileprefix');
end

if findstr('.coo',fileprefix),
    fprintf('MESH_BET2MATLAB: Removing .coo extension from ''%s''\n',fileprefix);
    fileprefix = strrep(fileprefix,'.coo','');
end
if findstr('.dat',fileprefix),
    fprintf('MESH_BET2MATLAB: Removing .dat extension from ''%s''\n',fileprefix);
    fileprefix = strrep(fileprefix,'.dat','');
end

FV.vertices = readCOO(fileprefix);

Nvertices   = size(FV.vertices,1);

edgematrix  = readDAT(fileprefix,Nvertices);

FV.faces    = edge2face(FV.vertices,edgematrix);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function faces = edge2face(vertices,edgematrix)
    
    Nvertices = size(vertices,1);
    
    plot3(vertices(:,1),vertices(:,2),vertices(:,3),'b.'), hold on
    
    for i = 1:Nvertices,
        
        % Get current vertex XYZ coordinates
        v  = vertices(i,:);
        
        plot3(v(:,1),v(:,2),v(:,3),'go');
        
        % Get neighbouring vertex XYZ coordinates
        vn = vertices(find(edgematrix(i,:)),:);
        
        plot3(vn(:,1),vn(:,2),vn(:,3),'ro');
        
        
        faces = tri([v;vn]);
        
        
        
        
        return
        
        
        
        
        
    end
    
    faces = [];
    
return


function faces = tri(V)
    
    % Must identify angle of rotation from
    % central vertex, V(1,:), to all other vertices
    % in a counterclockwise direction
    
    % First translate all vertices to an origin
    % given by the central vertex
    Vo = V - repmat(V(1,:),size(V,1),1);
    
    % Find direction cosines for line from centre to vertex
    d = zeros(size(Vo,1),1);
    l = zeros(size(Vo,1),1);
    m = zeros(size(Vo,1),1);
    n = zeros(size(Vo,1),1);
    
    for i = 1:length(Vo),
        
        x = Vo(i,1);  y = Vo(i,2);  z = Vo(i,3);
        
        d(i) = sqrt( (x)^2 + (y)^2 + (z)^2 );
        
        if d(i) > 0,
            l(i) = x/d(i); % cos alpha
            m(i) = y/d(i); % cos beta
            n(i) = z/d(i); % cos gamma
        end
    end
    
    
    faces = [];
    
    
    %hold off
    %plot3(Vo(:,1),Vo(:,2),Vo(:,3),'b.'), view(2), hold on
    
    
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vertices = readCOO(fileprefix)
	
	[path,name,ext] = fileparts(strcat(fileprefix,'.coo'));
	file = fullfile(path,[name ext]);
	
	fid = fopen(file,'r');
	
	if isequal(fid,-1),
        S=sprintf('Could not open file: "%s"',file);
        error(S);
	else
        fprintf('...Reading BET vertices...');
        tic;
        
        % Read vertices
        vertices = fscanf(fid,'%s%f%f%f',[4,inf]);
        fclose(fid);
        
        % remove last row (all zeros) and translate
        vertices = vertices(2:4,:)';
        
        t = toc;
        fprintf('...done (%6.2f sec).\n',t);
	end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edgematrix = readDAT(fileprefix,Nvertices)
	
	[path,name,ext] = fileparts(strcat(fileprefix,'.dat'));
	file = fullfile(path,[name ext]);
	
	fid = fopen(file,'r');
	
	if isequal(fid,-1),
        S=sprintf('Could not open file: "%s"',file);
        error(S);
	else
        fprintf('...Reading BET edge matrix...');
        tic;
        
        % Read faces
        edgematrix = fscanf(fid,'%1d',[Nvertices,Nvertices]);
        fclose(fid);
        
        t = toc;
        fprintf('...done (%6.2f sec).\n',t);
	end
return
