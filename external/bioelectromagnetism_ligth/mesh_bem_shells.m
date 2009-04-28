function[p] = mesh_bem_shells(IMG,FS)

% mesh_bem_shells - Returns BEM scalp, skull and cortex meshes
% 
%[p] = mesh_bem_shells(IMG,FS)
% 
% MRI volumes to process, eg:
% IMG.path   = '\data\freesurfer\subjects\c01\mri\analyze\';
% IMG.scalp  = 'c01_orig_axial_ras';
% IMG.oskull = 'c01_orig_axial_ras_bet_skull';
% IMG.iskull = 'c01_orig_axial_ras_bet';
% IMG.intensity.scalp = 
% 
% The *_bet_skull volume should be converted to a binary mask
% and a mask of the *_bet volume should be subtracted from it.
% In addition, consider dilating the *_bet_skull volume and
% multiplying it by a FAST segmentation of the skull/csf (seg0)
% to remove extraneous data.
% 
% Freesurfer surface to process, eg:
% FS.path   = '\data\freesurfer\subjects\c01\surf\';
% FS.file   = 'rh.white.asc';
% 
% Except for the scalp, this function assumes extensive 
% MRI preprocessing with:
% 
% FreeSurfer         http://surfer.nmr.mgh.harvard.edu/
% FSL tools (BET)    http://fmrib.ox.ac.uk/
% 
% The return meshes are in p.mesh.data
% 
% If you only want a scalp surface from an MRI volume,
% you can copy the scalp section of the code from this 
% function into your own script.
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:56 $

% Licence:  GNU GPL, no implied or express warranties
% History:  10/2002, Darren.Weber_at_radiology.ucsf.edu, created
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = eeg_toolbox_defaults;
p.mesh.path = IMG.path;
p.mesh.file = 'BrainStorm_subjecttess.mat';

fprintf('\n\n**************************************\nProcessing %s\n\n',IMG.path);


% LOAD FREESURFER CORTEX
cd(FS.path)
FS.surf = mesh_freesurfer2matlab(FS.file);
Nsurf = 1;
if isfield(FS.surf,'vertices'),
    p.mesh.data.meshtype{Nsurf} = 'cortex';
    p.mesh.data.vertices{Nsurf} = FS.surf.vertices;
    p.mesh.data.faces{Nsurf}    = FS.surf.faces;
    %patch('vertices',FS.surf.vertices,'faces',FS.surf.faces,'FaceColor',[.8 .0 .0],'Edgecolor','none');
end

% TESSELATE SKULL/SCALP SURFACES

cd(IMG.path)

% --- Tesselate Inner skull
fprintf('\n\n**************************************\nProcessing "Inner Skull"\n\n');
Nsurf = 2;
intensity = IMG.intensity.iskull;
tolerance = IMG.tolerance.iskull;
avw = avw_img_read(IMG.iskull);
ISkull = mesh_shrinkwrap(avw.img,[],0,0,intensity,tolerance,50,0.5,0.4);
if isfield(ISkull,'vertices'),
    p.mesh.data.meshtype{Nsurf} = 'inner_skull';
    p.mesh.data.vertices{Nsurf} = ISkull.vertices;
    p.mesh.data.faces{Nsurf}    = ISkull.faces;
    %patch('vertices',ISkull.vertices,'faces',ISkull.faces,'FaceColor',[.6 .0 .0],'Edgecolor',[.8 .8 .8],'FaceAlpha',.2); axis tight; hold on
end

% --- Tesselate Outer skull
fprintf('\n\n**************************************\nProcessing "Outer Skull"\n\n');
Nsurf = 3;
intensity = IMG.intensity.oskull;
tolerance = IMG.tolerance.oskull;
avw = avw_img_read(IMG.oskull);
OSkull = mesh_shrinkwrap(avw.img,[],1,0,intensity,tolerance,10,0.005,0.1);
if isfield(OSkull,'vertices'),
    p.mesh.data.meshtype{Nsurf} = 'outer_skull';
    p.mesh.data.vertices{Nsurf} = OSkull.vertices;
    p.mesh.data.faces{Nsurf}    = OSkull.faces;
    %patch('vertices',OSkull.vertices,'faces',OSkull.faces,'FaceColor',[.6 .0 .0],'Edgecolor',[.8 .8 .8],'FaceAlpha',.2); axis tight; hold on
end

% --- Tesselate Scalp ---
fprintf('\n\n**************************************\nProcessing Scalp\n\n');
Nsurf = 4;
intensity = IMG.intensity.scalp;
tolerance = IMG.tolerance.scalp;
avw = avw_img_read(IMG.scalp);
scalp = mesh_shrinkwrap(avw.img,[],0,0,intensity,tolerance,50,0.5,0.4);
if isfield(scalp,'vertices'),
    p.mesh.data.meshtype{Nsurf} = 'scalp';
    p.mesh.data.vertices{Nsurf} = scalp.vertices;
    p.mesh.data.faces{Nsurf}    = scalp.faces;
    %patch('vertices',scalp.vertices,'faces',scalp.faces,'FaceColor',[.6 .0 .0],'Edgecolor',[.8 .8 .8],'FaceAlpha',.2); axis tight; hold on
end


return
