function [p] = mesh_fit_elec(p)

% mesh_fit_elec - find mesh vertices nearest to electrodes
% 
% [p] = mesh_fit_elec(p)
% 
% This function is in development.  A transformation matrix,
% with rotations and translations, is in development to
% facilitate coregistration of elec vertices to the scalp 
% vertices.  This needs to be refined by a minimisation 
% function for the difference between the electrode 
% locations and the nearest scalp vertices. The scalp 
% mesh should be sufficiently refined to provide unique 
% vertices for each electrode vertex.
% 
% At present, it works with the BrainStorm meshes and the
% elec_124_cart.txt file in the eeg_example_data folder.
% The BrainStorm meshes contain a 'scalp' mesh with vertices
% in meters and the electrode file contains position 
% vectors in centimeters.  The electrodes are converted
% to meters by ELEC_LOAD.  In the example data, they were 
% previously coregistered and fitted to the scalp mesh.
% 
% The return vertices of the scalp that are nearest to
% the electrodes are returned in the following:
% 
% p.mesh.data.meshtype{p.mesh.current}
% p.mesh.data.vertices{p.mesh.current}
% p.mesh.data.faces{p.mesh.current}
% 
% where p.mesh.current is also returned.
% 
% See also: fiducial_coregister
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:57 $

% Licence:  GNU GPL, no express or implied warranties
% History:  04/2002, Darren.Weber_at_radiology.ucsf.edu
%           09/2002, Darren.Weber_at_radiology.ucsf.edu
%                    - bug fixing location of scalp mesh vertices
%                    that correspond to electrode vertices.
%                    - Still working with example dataset, not
%                    a genuine coregistration algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% debugging on/off
global DB;
DB = 0;
test = 0;  % coregistration testing



fprintf('MESH_FIT_ELEC: In development, not yet a true coregistration\n');

p.mesh.current = mesh_check(p,'scalp');

if isempty(p.mesh.current),
    error('MESH_FIT_ELEC: No ''scalp'' mesh in p struct.\n');
end

% Extract scalp triangulation
FV.vertices = p.mesh.data.vertices{p.mesh.current};
FV.faces    = p.mesh.data.faces{p.mesh.current};

if isempty(p.elec.data),
    error('MESH_FIT_ELEC: ''p.elec.data'' is empty.\n');
else
    x = p.elec.data.x;
    y = p.elec.data.y;
    z = p.elec.data.z;
end

fprintf('MESH_FIT_ELEC: ''Coregistering'' Electrodes...');
tic

if test,
	% Electrode fiducial points
	En = p.elec.data.nasion;
	Er = p.elec.data.rpa;
	El = p.elec.data.lpa;
	reg.fiducials.head = [ En; Er; El ];
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% NEED AN INTERACTIVE PROCESS TO SELECT THESE
	% Mesh fiducial points
	if isempty(p.mri.fiducials),
        if ~isempty(p.mri.data),
            fprintf('\n\nMESH_FIT_ELEC: Select Fiducial Points in MRI...\n\n');
            avw_view(p.mri.data);
            waitfor(gcf);
            p.mri.fiducials = evalin('base',p.mri.fiducials);
		else
            fprintf('\n\nMESH_FIT_ELEC: Select Fiducial Points in MRI...\n\n');
           [p] = mri_open(p);
            avw_view(p.mri.data);
            waitfor(gcf);
            p.mri.fiducials = evalin('base',p.mri.fiducials);
		end
	end
	reg.fiducials.mri = p.mri.fiducials;
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	T = fiducial_coregister(reg.fiducials.head,reg.fiducials.mri);
    
    E = [p.elec.data.x p.elec.data.y p.elec.data.z];
    E2MRI = [E,ones(size(E,1),1)] * T;
end


% perform least squares fit of elec to scalp???
%options = optimset('fminsearch');
%[sumd,fval,exitflag,output] = fminsearch('mesh_fit_elec_optim',...
%                              0, options, vertices, [x,y,z]);
%fprintf('\n%s%f\n', 'Iterations = ', output.iterations);
%fprintf('\n%s%d\n', 'Exit = ', exitflag);




% Try to adjust the height of the electrodes to
% those of the scalp mesh - really need a full
% translation/rotation/scaling least squares solution
elecmax = max(z);
vertmax = max(FV.vertices(:,3));
z = z + (vertmax - elecmax);


if DB,
    figure; plot3(x,y,z,'ro'); hold on;
    patch('vertices',FV.vertices,'faces',FV.faces,...
          'FaceColor',[0.0 0.0 0.6],'Edgecolor',[.6 .6 .6]);
    set(gca,'Projection','perspective'); 
    set(gca,'DataAspectRatio',[1 1 1]); 
    axis off tight vis3d
end


% Find nearest scalp vertices to electrode positions.
% In the process, reorder the scalp mesh vertices 
% and correct the scalp mesh faces so that all scalp
% vertices that lie near the electrodes
% are numbered 1-N, in the same order as the
% electrodes.  This is done here so that the
% mesh_laplacian_interp function will work
% correctly.

% This mess might be better handled by adding
% electrode vertices into the scalp mesh, rather
% than finding the scalp mesh vertices nearest
% to those of the electrodes.  However, this would
% require retesselation of the scalp mesh, which
% is a tricky process

[FV,k] = mesh_nearest_reorder(FV,[x,y,z]);


% Assign the adjusted triangulation back
% into the original structure
p.mesh.data.vertices{p.mesh.current} = FV.vertices;
p.mesh.data.faces{p.mesh.current}    = FV.faces;


p.mesh.data.elecindex{p.mesh.current} = k;

x = FV.vertices(k,1);
y = FV.vertices(k,2);
z = FV.vertices(k,3);

% A quick triangulation of the scalp vertices
% that lie nearest to the electrode vertices
Efaces = convhulln(FV.vertices(k,:));

% Now store the vertices
[p.mesh.current,meshExists] = mesh_check(p,'elec');

p.mesh.data.meshtype{p.mesh.current} = 'elec';
p.mesh.data.vertices{p.mesh.current} = [x,y,z];
p.mesh.data.faces{p.mesh.current} = Efaces;


t = toc; fprintf('done (%6.2f sec).\n',t);

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FV,k] = mesh_nearest_reorder(FV,elecverts)
    
    % Find the indices of the nearest vertices in
    % FV to those in vertices
    k = dsearchn(FV.vertices,elecverts);
    
    if size(k,1) ~= size(elecverts,1),
        msg = sprintf('Sorry, duplicate vertices found in electrode/mesh coregistration');
        error(msg);
    end
    
    
    % Debug visualization
    global DB;
    if DB,
        figure; hold on;
        patch('vertices',FV.vertices,'faces',FV.faces,...
              'FaceColor',[0.8 0.8 0.8],'Edgecolor',[.6 .6 .6]);
        set(gca,'Projection','perspective');
        set(gca,'DataAspectRatio',[1 1 1]);
        axis off tight vis3d
    end
    
    
    % OK, now reorder the triangulation in FV
    % so that all the mesh vertices are ordered
    % 1:N
    
    Nelec = size(elecverts,1);
    
    for Eindex = 1:Nelec,
        
        % It's important in this process that the
        % mesh indices are found over again, in case
        % they are modified for electrode Eindex-X
        elecpos = elecverts(Eindex,:);
        Kindex = dsearchn(FV.vertices,elecpos);
        
        if Eindex == Kindex,
            continue;  % great, do nothing
        end
        
        % Get vertex coordinates at Eindex
        Evertex = FV.vertices(Eindex,:);
        % Get vertex coordinates at Kindex
        Kvertex = FV.vertices(Kindex,:);
        
        if DB,
            fprintf('Fitting Electrode: %3d\n',Eindex);
            H1 = plot3(Evertex(1),Evertex(2),Evertex(3),'ro');
            H2 = plot3(elecpos(1),elecpos(2),elecpos(3),'rd');
            H3 = plot3(Kvertex(1),Kvertex(2),Kvertex(3),'bo');
            pause
            delete(H1)
            delete(H2)
            delete(H3)
        end
        
        % Now swap them
        FV.vertices(Kindex,:) = Evertex;
        FV.vertices(Eindex,:) = Kvertex;
        
        % Swapping the vertices requires
        % rearranging the faces
        
        % find all faces that contain Vindex
        EfaceIndices = find(FV.faces == Eindex);
        % find all faces that contain Kindex
        KfaceIndices = find(FV.faces == Kindex);
        
        % Where faces refer to Kvertex, make them
        % refer to Vvertex and vice versa
        FV.faces(KfaceIndices) = Eindex;
        FV.faces(EfaceIndices) = Kindex;
        
        % Debug testing....
        if DB,
            % Get vertex coordinates at Vindex
            Evertex = FV.vertices(Eindex,:);
            % Get vertex coordinates at Kindex
            Kvertex = FV.vertices(Kindex,:);
            H1 = plot3(Evertex(1),Evertex(2),Evertex(3),'ro');
            H2 = plot3(elecpos(1),elecpos(2),elecpos(3),'rd');
            H3 = plot3(Kvertex(1),Kvertex(2),Kvertex(3),'bo');
            pause
            
            delete(H1)
            delete(H2)
            delete(H3)
        end
        
    end
    
    % Now the nearest vertices in FV to vertices
    % should be numbered 1:Nelec
    k = dsearchn(FV.vertices,elecverts);
    
    if ~isequal(k,[1:Nelec]'),
        msg = sprintf('MESH_FIT_ELEC: Mesh vertices not in order of electrodes!');
        warning(msg);
    end
    
return
