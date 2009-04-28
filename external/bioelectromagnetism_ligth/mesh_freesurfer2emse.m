function [FV] = mesh_freesurfer2emse(sfile,rfile)

% mesh_freesurfer2emse - Convert FreeSurfer tesselation to EMSE
%
% USEAGE: [FV] = mesh_freesurfer2emse(sfile,rfile)
%
% This function will call several functions to (a) read a
% freesurfer surface (sfile), (b) read an EMSE registration
% file (rfile) and (c) output the surface as an EMSE wireframe
% file (*.wfr) in the registration space.
% 
% The freesurfer surface is read with mesh_freesurfer2matlab.
% The EMSE registration file is read with emse_read_reg.
% The EMSE wireframe is writen with mesh_writ
% Websites:
% 
% EMSE:       http://www.sourcesignal.com/
% Freesurfer: http://surfer.nmr.mgh.harvard.edu/
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:57 $

% Licence:  GNU GPL, no implied or express warranties
% History:  08/2002, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



msg = sprintf('MESH_FREESURFER2EMSE: In development');
%error(msg);



datapath = '\\potzii\data\freesurfer\subjects\ptsdpet-c01keep\';
cd(datapath)

sfile = 'surf\rh.pial.asc';
rfile = 'mri\analyze\c01_orig_axial_ras.reg';
mfile = 'mri\analyze\c01_orig_axial_ras.img';

% First read the freesurfer surface
[FV] = mesh_freesurfer2matlab(sfile);

% Convert from mm to meters
%FV.vertices = FV.vertices ./ 1000;

avw = avw_img_read(mfile);



keyboard




% Now read the EMSE registration file
%REG = emse_read_reg(rfile);
%scale = ones(size(FV.vertices,1),1);
%vert = [FV.vertices, scale];
%vert = vert * REG.mri2elec;
%vert = Rz(vert(:,1:3),180,'degrees');
%FV.vertices = vert(:,1:3);

%EMSE.file = '\\POTZII\data\freesurfer\subjects\ptsdpet-c01keep\mri\analyze\c01_orig_axial_ras_scalp.wfr';
%[vertices,faces,edges,meshtype] = mesh_emse2matlab(EMSE.file,{'vertex','face'});
%EMSE.vertices = [vertices.x; vertices.y; vertices.z]';
%EMSE.faces = [faces.vertex1;faces.vertex2;faces.vertex3]';
% Rotate 90 degrees around Z for EMSE data
%EMSE.vertices = Rz(EMSE.vertices,90,'degrees');

%clear vertices faces edges meshtype;

% Plot EMSE scalp surface
%figure; hold on
%patch('vertices',EMSE.vertices,'faces',EMSE.faces,'FaceColor',[.6 .0 .0],'Edgecolor','none','FaceAlpha',.4);

% Plot freesurfer cortex and transformed vertices
%patch('vertices',FV.vertices,'faces',FV.faces,'FaceColor',[.0 .0 .6],'Edgecolor','none');
%patch('vertices',vert,'faces',FV.faces,'FaceColor',[.6 .0 .0],'Edgecolor','none');
%set(gca,'Projection','perspective');
%set(gca,'DataAspectRatio',[1 1 1]);
%axis off tight vis3d
%light
%mouse_rotate



%p.mesh.path = '\\potzii\data\freesurfer\subjects\ptsdpet-c01keep\mri\analyze\';
%p.mesh.file = 'c01_orig_axial_ras.wfr';
%p.mesh.data.meshtype{1} = 'pial';
%p.mesh.data.vertices{1} = FV.vertices;
%p.mesh.data.faces{1}    = FV.faces;
%mesh_write_emse(p);

return




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
