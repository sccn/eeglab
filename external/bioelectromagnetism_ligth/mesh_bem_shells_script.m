
% mesh_bem_shells_script - create bem shells for N subjects
% 
% Assumes preprocessing of MRI with 
% FreeSurfer                http://surfer.nmr.mgh.harvard.edu/
% FSL tools (BET & FAST)    http://fmrib.ox.ac.uk/

% 
% Licence:  GNU GPL, no implied or express warranties
% History:  08/2002, Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



S = 'c01';

% These options are control flow points below
meshplot    = 0;
elecplot    = 0;
getfid      = 0;
coregister  = 0;

% Create and initial correction of BEM
% meshes from MRI and FreeSurfer data
create      = 0;
% Refinement of BEM correction process?
correct     = 0;

% default output is brainstorm format,
% but these can be output also:
% output BEM in EMSE format?
emse        = 0;
% output BEM in FreeSurfer format?
freesurfer  = 1;



% Fiducial points in RAS volume, obtained using avw_view (in meters)
% 
%                  Nasion                 Right                Left
%
mriFID.sub{ 1} = 'c01';
mriFID.xyz{ 1} = [  0.005  0.097 -0.022;  0.071  0.023 -0.052; -0.072  0.016 -0.050 ];
mriFID.sub{ 2} = 'c02';
mriFID.xyz{ 2} = [  0.007  0.088 -0.007;  0.081  0.015 -0.046; -0.070  0.012 -0.051 ];
mriFID.sub{ 3} = 'c03';
mriFID.xyz{ 3} = [ -0.002  0.098 -0.008;  0.055  0.040 -0.038; -0.095  0.030 -0.046 ];
mriFID.sub{ 4} = 'c04';
mriFID.xyz{ 4} = [  0.000  0.092  0.000;  0.070  0.024 -0.036; -0.076  0.031 -0.036 ];
mriFID.sub{ 5} = 'c05';
mriFID.xyz{ 5} = [  0.010  0.094 -0.010;  0.075  0.024 -0.048; -0.072  0.024 -0.052 ];
mriFID.sub{ 6} = 'c06';
mriFID.xyz{ 6} = [ -0.002  0.088  0.000;  0.077  0.008 -0.048; -0.075 -0.008 -0.048 ];
mriFID.sub{ 7} = 'c07';
mriFID.xyz{ 7} = [ -0.001  0.097 -0.018;  0.078  0.008 -0.055; -0.076  0.008 -0.055 ];
mriFID.sub{ 8} = 'c08';
mriFID.xyz{ 8} = [  0.000  0.096 -0.019;  0.087  0.005 -0.039; -0.065  0.004 -0.038 ];
mriFID.sub{ 9} = 'c09';
mriFID.xyz{ 9} = [  0.010  0.086 -0.020;  0.100 -0.020 -0.080; -0.060  0.020 -0.058 ];
mriFID.sub{10} = 'c10';
mriFID.xyz{10} = [  0.000  0.091 -0.018;  0.076  0.011 -0.052; -0.068  0.006 -0.052 ];

mriFID.sub{11} = 'p02';
mriFID.xyz{11} = [  0.002  0.091  0.016;  0.080  0.025 -0.036; -0.080  0.015 -0.043 ];
mriFID.sub{12} = 'p04';
mriFID.xyz{12} = [ -0.001  0.086 -0.009;  0.086  0.002 -0.062; -0.070  0.003 -0.065 ];
mriFID.sub{13} = 'p05';
mriFID.xyz{13} = [ -0.002  0.094 -0.001;  0.068  0.025 -0.040; -0.071  0.005 -0.042 ];
mriFID.sub{14} = 'p06';
mriFID.xyz{14} = [ -0.002  0.084  0.000;  0.082  0.013 -0.050; -0.066  0.013 -0.052 ];
mriFID.sub{15} = 'p07';
mriFID.xyz{15} = [  0.001  0.092  0.015;  0.080  0.003 -0.033; -0.070  0.004 -0.032 ];
mriFID.sub{16} = 'p08';
mriFID.xyz{16} = [ -0.003  0.095 -0.002;  0.070  0.018 -0.040; -0.074  0.022 -0.035 ];
mriFID.sub{17} = 'p09';
mriFID.xyz{17} = [ -0.002  0.100  0.004;  0.100  0.002 -0.028; -0.050  0.000 -0.036 ];





if meshplot & ~correct,
    
    patch('vertices',p.mesh.data.vertices{4},'faces',p.mesh.data.faces{4},...
          'FaceColor',[1 0 0],'Edgecolor','none','FaceAlpha',.3);
    
    daspect([1 1 1]); axis tight; hold on
    
    patch('vertices',p.mesh.data.vertices{3},'faces',p.mesh.data.faces{3},...
          'FaceColor',[0 1 0],'Edgecolor','none','FaceAlpha',.5);
    
    patch('vertices',p.mesh.data.vertices{2},'faces',p.mesh.data.faces{2},...
          'FaceColor',[0 0 1],'Edgecolor','none','FaceAlpha',1.0);
    
    material dull
    camlight headlight
    
    mouse_rotate
    
    return
    
    % Check that vertices for scalp, oskull, iskull are coincident
    v = 1;
    
    vert = p.mesh.data.vertices{4};
    x = vert(v,1);
    y = vert(v,2);
    z = vert(v,3);
    plot3(x,y,z,'ro')
    
    vert = p.mesh.data.vertices{3};
    x = vert(v,1);
    y = vert(v,2);
    z = vert(v,3);
    plot3(x,y,z,'go')
    
    vert = p.mesh.data.vertices{2};
    x = vert(v,1);
    y = vert(v,2);
    z = vert(v,3);
    plot3(x,y,z,'bo')
    
    clear vert v x y z
    
    return
end


if elecplot,
    
    
    patch('vertices',p.mesh.data.vertices{4},'faces',p.mesh.data.faces{4},...
          'FaceColor',[1 0 0],'Edgecolor','none','FaceAlpha',.6);
    
    lighting phong
    material dull
    camlight headlight
    
    hold on
    
    x = p.elec.data.x;
    y = p.elec.data.y;
    z = p.elec.data.z;
    plot3(x,y,z,'bo')
    
%     x = E2MRI(:,1);
%     y = E2MRI(:,2);
%     z = E2MRI(:,3);
%     plot3(x,y,z,'bo')
    
    daspect([1 1 1]); axis tight; hold on
    
    mouse_rotate
    
    return
    
end


% 'c01', 'c02', 'c03', 'c04', 'c05', 'c06', 'c07', 'c08', 'c09', 'c10', 'p02', 'p04', 'p05', 'p06', 'p07', 'p08', 'p09'

for sub = { S },
    
    % Image files to process
    
    % These images were created with the freesurfer mri_convert command, eg:
    % mri_convert -oid 1 0 0 -ojd 0 1 0 -okd 0 0 1 orig mri\analyze\c01_orig_axial_ras.img
    % These Analyze files were then processed with FSL tools to find the skull.
    IMG.path   = sprintf('\\\\POTZII\\data\\freesurfer\\subjects\\ptsdpet-%s\\mri\\analyze\\',char(sub));
    IMG.scalp  = sprintf('%s_orig_axial_ras',char(sub));
    IMG.oskull = sprintf('%s_orig_axial_ras_skull',char(sub));
    IMG.iskull = sprintf('%s_orig_axial_ras_bet',char(sub));
    IMG.intensity.scalp  = 100;
    IMG.intensity.oskull = 0.05;
    IMG.intensity.iskull = 150;
    IMG.tolerance.scalp  = 20;
    IMG.tolerance.oskull = 0.01;
    IMG.tolerance.iskull = 80;
    
    % Freesurfer surface to process
    % This freesurfer surface contains the whole brain surface.
    FS.path   = sprintf('\\\\POTZII\\data\\freesurfer\\subjects\\ptsdpet-%s\\surf\\',char(sub));
    FS.file   = 'rh.pial.asc';
    
    % Create the meshes and correct them (use pial cortex for this)
    if create,
       [p] = mesh_bem_shells(IMG,FS);
        p.mesh.path = FS.path;
        p.mesh.file = 'BS_subjecttess.mat';
        mesh_write(p);
       [p] = mesh_bem_correct(p);
        p.mesh.file = 'BS_corrected_subjecttess.mat';
        mesh_write(p);
        
        clear p
    end
    
    
    
    % Open already processed surfaces
   [p] = mesh_open; % init p struct
%     p.mesh.path = FS.path;
%     p.mesh.file = 'BS_subjecttess.mat';
%     p.mesh.file = 'BS_corrected_subjecttess.mat';
    p.mesh.path = sprintf('\\\\POTZII\\data\\data_source\\%s\\meshes\\',char(sub));
    p.mesh.file = sprintf('%s_subjecttess.mat',char(sub));
   [p] = mesh_open(p);
    
    
% 	% LOAD SMOOTHWM AS CORTEX
% 	FS.file   = 'rh.smoothwm.asc';
% 	cd(FS.path)
% 	FS.surf = mesh_freesurfer2matlab(FS.file);
% 	Nsurf = 1;
% 	if isfield(FS.surf,'vertices'),
%         p.mesh.data.meshtype{Nsurf} = 'cortex';
%         p.mesh.data.vertices{Nsurf} = FS.surf.vertices;
%         p.mesh.data.faces{Nsurf}    = FS.surf.faces;
%         %patch('vertices',FS.surf.vertices,'faces',FS.surf.faces,'FaceColor',[.8 .0 .0],'Edgecolor','none');
% 	end
% 	
% 	% Convert from mm to meters
% 	for m = 1:size(p.mesh.data.meshtype,2),
%         p.mesh.data.vertices{m} = p.mesh.data.vertices{m} ./ 1000;
% 	end
% 	p.mesh.path = sprintf('\\\\POTZII\\data\\data_source\\%s\\meshes\\',char(sub));
% 	p.mesh.file = sprintf('%s_subjecttess.mat',char(sub));
% 	mesh_write(p);
    
    
    
    % Use manual revision to refine the dist struct and recorrect BEM
    if correct,
        dist.four2three = 5;    % scalp to outer skull
        dist.four2two   = 4;    % scalp to inner skull
        dist.three2two  = 3;    % outer skull to inner skull
        dist.two2one    = 1;    % inner skull to cortex
       [p] = mesh_bem_correct(p,[],dist,[],0);
    end
    
    
    
    
    % COREGISTER ELECTRODES TO BEM
    if coregister,
        
        % Use avw_view to update p.mri.fiducials in the base workspace
        % Then copy those values to the mriFID struct above and
        % adjust as necessary
        if getfid,
            cd(IMG.path)
            avw = avw_img_read(IMG.scalp);
            avw_view(avw);
            return
        end
        
        % Load the electrode data
        p.elec.path = sprintf('\\\\POTZII\\data\\data_source\\%s\\meshes\\',char(sub));
        p.elec.file = sprintf('%s_124fit.txt',char(sub));
        p.elec.plot = 0;
       [p] = elec_open(p);
        
        % Get surface fiducials from mriFID struct above
        Nfid = strmatch(sub,mriFID.sub);
        
        p.mri.fiducials = mriFID.xyz{Nfid}; % 3x3, nasion, right, left fiducial points in rows
        % Get electrode fiducials
        Efid = [p.elec.data.nasion; p.elec.data.rpa; p.elec.data.lpa];
        % Calculate coregistration transform
        T = fiducial_coregister(Efid,p.mri.fiducials);
        
        [path,file,ext] = fileparts([p.elec.path p.elec.file]);
        file = fullfile(path,[file,'_reg.mat']);
        save(file,'T')
        
        % Transform all the electrode coordinates to MRI space
        E = [p.elec.data.x p.elec.data.y p.elec.data.z];
        E2MRI = [E,ones(size(E,1),1)] * T;
        
        p.elec.data.x = E2MRI(:,1);
        p.elec.data.y = E2MRI(:,2);
        p.elec.data.z = E2MRI(:,3);
        
        REF2MRI = [p.elec.data.ref 1] * T;
        p.elec.data.ref    = REF2MRI(:,1:3);
        
        ORIGIN2MRI = [p.elec.data.centroid 1] * T;
        p.elec.data.centroid = ORIGIN2MRI(:,1:3);
        
        Efid2MRI = [ Efid [ 1 1 1]' ] * T;
        p.elec.data.nasion = Efid2MRI(1,1:3);
        p.elec.data.rpa    = Efid2MRI(2,1:3);
        p.elec.data.lpa    = Efid2MRI(3,1:3);
        
        
        clear E REF2MRI ORIGIN2MRI Efid2MRI
        
        %elec_write_3dspace(p);
        %elec_write_brainstorm(p);
        
        p.elec.file = sprintf('%s_124fit_new.elp',char(sub));
        elec_write_emse(p);
        
    end
    
    % Output BEM surfaces as EMSE format
    if emse,
        p.mesh.file = sprintf('%s.wfr',char(sub));
        p.mesh.type = 'emse';
        mesh_write(p);
    end
    
    % Output BEM surfaces as freesurfer format
    if freesurfer,
        p.mesh.path = FS.path;
        p.mesh.file = 'rh.asc';
        p.mesh.type = 'freesurfer';
        p.mesh.data.meshtype{1} = ''; % do not write cortex
        mesh_write(p);
    end
    
end

return
