
% mesh_3shell_script
% 
% Assumes preprocessing of MRI with 
% FreeSurfer                http://surfer.nmr.mgh.harvard.edu/
% FSL tools (BET & FAST)    http://fmrib.ox.ac.uk/
% 
% An example FSL script can be found at the end of this script
%

% Licence:  GNU GPL, no implied or express warranties
% History:  08/2002, Darren.Weber_at_radiology.ucsf.edu
%           08/2003, Darren.Weber_at_radiology.ucsf.edu
%                    added more control variables and the FSL script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% the plot option is exclusive of the following commands, 
% it should be done after initial attempt to get scalp etc.
plot = 0; 

% the initial results may not be very pretty, these tools
% are not ideal!  It usually takes some time to tweak the
% shrinkwrap results and then finally use bem_shells_correct to 
% get reasonable BEM surfaces.  All development work needs to
% be done on modifying the shrinkwrap method (it's not great!).
doScalp  = 0;
doIskull = 0;
doOskull = 0;
doCortex = 0;

% only use this after initial processing above, use the plot
% option to check the initial results.  This option will use
% bem_shells_correct, which can do a neat job of correcting
% the surfaces for intersections.
doCorrection = 0;

doOutputEMSE = 0;
doOutputBrainStorm = 0;

%IMG.path   = 'D:\freesurfer\subjects\ucsf_morgan\mri\analyze';
IMG.path     = '/data/ucsf/mri/ucsf_morgan/mri/analyze';
IMG.meshpath = '/data/ucsf/mri/ucsf_morgan/surf';

IMG.scalp  = 'morgan_orig_axial_las_smooth';
IMG.oskull = 'morgan_orig_axial_las_skull';
IMG.iskull = 'morgan_orig_axial_las_brain';

IMG.intensity.scalp  =  60;
IMG.tolerance.scalp  =  20;
IMG.intensity.oskull =  15;
IMG.tolerance.oskull =  10;
IMG.intensity.iskull =  80;
IMG.tolerance.iskull =  40;

if ~exist('p','var'),
   [p] = eeg_toolbox_defaults;
end

if plot,
    
    % plot scalp surface (index 4), in red
    if ~isempty(p.mesh.data.vertices{4}),
        patch('vertices',p.mesh.data.vertices{4},'faces',p.mesh.data.faces{4},...
            'FaceColor',[.6 0 0],'Edgecolor','none','FaceAlpha',.2);
    end
    % plot outer skull surface (index 3), in green
    if ~isempty(p.mesh.data.vertices{3}),
        patch('vertices',p.mesh.data.vertices{3},'faces',p.mesh.data.faces{3},...
            'FaceColor',[.0 .6 .0],'Edgecolor',[.8 .8 .8],'FaceAlpha',.4);
    end
    % plot inner skull surface (index 2), in blue
    if ~isempty(p.mesh.data.vertices{2}),
        patch('vertices',p.mesh.data.vertices{2},'faces',p.mesh.data.faces{2},...
            'FaceColor',[.0 .0 .6],'Edgecolor',[.8 .8 .8],'FaceAlpha',.6);
    end
    % plot cortex surface (index 1), in black
    if ~isempty(p.mesh.data.vertices{1}),
        patch('vertices',p.mesh.data.vertices{1},'faces',p.mesh.data.faces{1},...
            'FaceColor',[.0 .0 .0],'Edgecolor',[.8 .8 .8],'FaceAlpha',1);
    end
    
    camlight; mouse_rotate;
    axis tight; hold on
    
% 	FV = ISkull;
% 	patch('vertices',FV.vertices,'faces',FV.faces,'FaceColor',[.6 .0 .0],'Edgecolor','none','FaceAlpha',.2); axis tight; hold on
% 	x = FV.vertices(1,1);
% 	y = FV.vertices(1,2);
% 	z = FV.vertices(1,3);
% 	plot3(x,y,z,'ro')
%   
    return
end



% TESSELATE SKULL/SCALP SURFACES


% Image files to process

% These images were created with the freesurfer mri_convert command, eg:
% mri_convert -oid -1 0 0 -ojd 0 1 0 -okd 0 0 1 orig mri\analyze\morgan_orig_axial_las.img
% These Analyze files were then processed with the FSL tools to find the
% skull.

cd(IMG.path)

if doScalp,
    
    % --- Tesselate Scalp ---
    fprintf('\n\n**************************************\nProcessing Scalp\n\n');
    Nsurf = 4;
    intensity = IMG.intensity.scalp;
    tolerance = IMG.tolerance.scalp;
    avw = avw_img_read(IMG.scalp);
% [FV, Edges] = avw_shrinkwrap(avw,FV,smooth,vthresh,interpVal,...
%               fitval,fittol,fititer,fitchange,fitvattr)
    scalp = avw_shrinkwrap(avw,[],0,0,[],intensity,tolerance,50,0.5,0.4);
    if isfield(scalp,'vertices'),
        p.mesh.data.meshtype{Nsurf} = 'scalp';
        p.mesh.data.vertices{Nsurf} = scalp.vertices ./ 1000; % convert to meters
        p.mesh.data.faces{Nsurf}    = scalp.faces;
    end
end


if doIskull,
    
    % --- Tesselate Inner skull
    fprintf('\n\n**************************************\nProcessing "Inner Skull"\n\n');
    Nsurf = 2;
    intensity = IMG.intensity.iskull;
    tolerance = IMG.tolerance.iskull;
    avw = avw_img_read(IMG.iskull);
% [FV, Edges] = avw_shrinkwrap(avw,FV,smooth,vthresh,interpVal,...
%               fitval,fittol,fititer,fitchange,fitvattr)
    ISkull = avw_shrinkwrap(avw,[],0,0,[],intensity,tolerance,50,0.5,0.4);
    if isfield(ISkull,'vertices'),
        p.mesh.data.meshtype{Nsurf} = 'inner_skull';
        p.mesh.data.vertices{Nsurf} = ISkull.vertices ./ 1000; % convert to meters
        p.mesh.data.faces{Nsurf}    = ISkull.faces;
    end
end


if doOskull,
    
    % --- Tesselate Outer skull
    fprintf('\n\n**************************************\nProcessing "Outer Skull"\n\n');
    Nsurf = 3;
    intensity = IMG.intensity.oskull;
    tolerance = IMG.tolerance.oskull;
    avw = avw_img_read(IMG.oskull);
% [FV, Edges] = avw_shrinkwrap(avw,FV,smooth,vthresh,interpVal,...
%               fitval,fittol,fititer,fitchange,fitvattr)
    OSkull = avw_shrinkwrap(avw,[],0,0,0,intensity,tolerance,50,0.05,0.1);
    if isfield(OSkull,'vertices'),
        p.mesh.data.meshtype{Nsurf} = 'outer_skull';
        p.mesh.data.vertices{Nsurf} = OSkull.vertices ./ 1000; % convert to meters
        p.mesh.data.faces{Nsurf}    = OSkull.faces;
        %patch('vertices',OSkull.vertices,'faces',OSkull.faces,...
        %      'FaceColor',[.6 .0 .0],'Edgecolor',[.8 .8 .8],...
        %      'FaceAlpha',.2); axis tight; hold on
    end
end

if doCortex,
    
    % LOAD FREESURFER CORTEX
    % Freesurfer surface to process
    % This freesurfer surface contains the whole brain surface, after
    % altering the freesurfer pons/corpus callosum cutting positions.
    p.mesh.path   = IMG.meshpath;
    p.mesh.file   = 'rh.pial';
    p.mesh.type   = 'fs_surf';
    
   [p] = mesh_open(p);
    
end

if doCorrection,
    
    dist.four2three = 6;    % scalp to outer skull
    dist.four2two   = 8;    % scalp to inner skull
    dist.three2two  = 5;    % outer skull to inner skull
    dist.two2one    = 2;    % inner skull to cortex
    
    cortex = 1; % process cortex (1 = yes, 0 = no)
   [p] = mesh_bem_correct(p, [], dist, 1, cortex);
    
end


if doOutputEMSE,
    mesh_write_emse(p);
end
if doOutputBrainStorm,
    mesh_write_brainstorm(p);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example FSL script
% Copy this text into a text file and chmod a+x that file.
% Remove the % matlab comments and everything before the next line,
% including these instructions.
% #! /bin/bash
% 
% 
% imagepath="."
% imagetype="_orig_axial_las"
% 
% sub="subject";
% img=`echo $sub$imagetype`;
% 
% printf "\nprocessing $img\n";
% 
% printf "\n5mm Guassian spatial smoothing: ${img}_smooth.*\n";
% /usr/local/fsl/bin/ip ${img} ${img}_smooth 0 -s 5
% 
% printf "\nUsing BET for brain extraction: ${img}_brain.*\n";
% bet ${img} ${img}_brain -s
% 
% # Binarise the bet skull
% #avwmaths ${img}_brain_skull -bin BETskull
% # Dilate bet skull and add a constant
% #avwmaths_32R BETskull -dil -add 10000 BETskull
% 
% printf "\nSkull extraction and 2mm Guassian spatial smoothing: ${img}_skull.*\n";
% /usr/local/fsl/bin/ip ${img}_brain_skull ${img}_skull 0 -s 2
% 
% printf "\nSpatial normalization\n";
% flirt -in ${img}_brain.hdr -ref /usr/local/fsl/etc/standard/avg152T1_brain.hdr -out ${img}_brain_reg.hdr -omat %{img}_brain_reg.mat -bins 256 -cost corratio -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp trilinear
% convert_xfm -matonly -omat ${img}_brain_reginv.mat -inverse ${img}_brain_reg.mat
% 
% printf "\nBrain segmentation\n";
% fast -t 1 -c 3 -op -e -ov -od $img ${img}_brain
% 
% exit
