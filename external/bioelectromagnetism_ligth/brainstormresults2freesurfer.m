% brainstormresults2freesurfer - script to extract brainstorm inverse to freesurfer

% Loop over all controls and patients

% Search ImageGridTime for selected time point
% Extract Cdata from ImageGridAmp at selected time (column)

%----------------------------------------------
% Define component timing

sa=[100 % first positive peak P100 (OT)
    150 % positive & negative peak P150 (SF) / N150 (OT)
    180 % negative N180 @ PT,AT and SP,SC
    250 % positive P250 (SP)
    400 ]; % P400 (SC,SP)

wm=[ 80
    90 
    150
    260
    300
    400
    550 ];

ea=[ 80 % N80  @ SF        , P80  @ OT/PT
    145   % N150 @ OT/PT     , P150 @ SPF (SC)
    240   %                  , P240 @ SPF/IPF & OT/SP
    300   % N300 @ L_PT/OT/IP,
    350   % Series: P350 @ Frontal (IPF), P500 @ SP
    450   % Series: P350 @ Frontal (IPF), P500 @ SP
    550 ];% Series: P350 @ Frontal (IPF), P500 @ SP; Also N550 @ left frontal

%----------------------------------------------
% Define groups
groups = {'c','p'};

%----------------------------------------------
% Define subjects
subs = [1:10];

cond = {'ea','sa','wm'};

for g=1:length(groups),
    for s=1:length(subs),
        
        sub = sprintf('%s%02d',groups{g},subs(s));
        
        %----------------------------------------------------
        % Load cortical tesselation from subjecttess.mat
        
        subtessfile = ['E:\matlab\brainstorm_v1\subjects\',sub,filesep,sub,'_subjecttess.mat'];
        
        if exist(subtessfile) ~= 2,
            % skip this group/subject, as no data available
            fprintf('\nskipping subject: %s\n\n',sub);
            continue;
        else
            p.mesh.path = ['E:\matlab\brainstorm_v1\subjects\',sub,filesep];
            p.mesh.file = [sub,'_subjecttess.mat'];
            p.mesh.type = 'BrainStorm';
           [p] = mesh_open(p);
            
            p.mesh.current = mesh_check(p,'cortex');
            
            % effectively remove all extraneous surfaces for mesh_write_freesurfer below
            indices = strcmp(p.mesh.data.meshtype,'cortex');
            exclude = find(indices == 0);
            p.mesh.data.meshtype(exclude) = [];
            
        end
        
        
        %----------------------------------------------------
        % define cortical source activity matrix file names
        
        %c01_oac_volts_data_results_image_svd12.mat
        %c01_ouc_volts_data_results_image_svd12.mat
        oucfile = ['E:\matlab\brainstorm_v1\studies\',sub,filesep,sub,'_ouc_volts_data_results_image_svd12.mat'];
        oacfile = ['E:\matlab\brainstorm_v1\studies\',sub,filesep,sub,'_oac_volts_data_results_image_svd12.mat'];
        oatfile = ['E:\matlab\brainstorm_v1\studies\',sub,filesep,sub,'_oat_volts_data_results_image_svd12.mat'];
        tacfile = ['E:\matlab\brainstorm_v1\studies\',sub,filesep,sub,'_tac_volts_data_results_image_svd12.mat'];
        
        if exist(oucfile) ~= 2, error(['no ',oucfile]); end;
        if exist(oacfile) ~= 2, error(['no ',oacfile]); end;
        if exist(oatfile) ~= 2, error(['no ',oatfile]); end;
        if exist(tacfile) ~= 2, error(['no ',tacfile]); end;
        
        %----------------------------------------------------
        % Save Cdata to freesurfer curvature files
        
        % D:\freesurfer\subjects\ptsdpet-c01\surf
        p.mesh.path = sprintf('D:\\freesurfer\\subjects\\ptsdpet-%s\\surf\\',sub);
        
        p.mesh.type = 'FS_curv'; % output freesurfer binary curvature files
        
        % cortical values are scaled into nano A.m (dipole moment)
        scale = 10^9;
        
        
        for c = 1:length(cond),
            
            
            switch cond{c},
                
            case 'sa',
                
                %---------------------------------------------------------------------
                % SA condition
                
                vert = size(p.mesh.data.vertices{p.mesh.current},1);
                times = length(sa);
                
                Dif.oac = zeros(vert,times);
                Dif.ouc = Dif.oac;
                
                load(oacfile,'ImageGridAmp','ImageGridTime'); % Time is in sec here
                for t=1:length(sa),
                    
                    % find column of ImageGridAmp corresponding to wm(t)
                    col = find(ImageGridTime == sa(t)/1000);
                    if isempty(col), error('Cannot find time column in ImageGridTime\n'); end
                    
                    p.mesh.data.Cdata{p.mesh.current} = ImageGridAmp(:,col) .* scale;
                    Dif.oac(:,t) = p.mesh.data.Cdata{p.mesh.current};
                    
                    p.mesh.file = sprintf('rh.sa_%04d_oac_float',sa(t));
                    mesh_write(p);
                    
                end
                
                clear ImageGridAmp ImageGridTime; pack;
                
                load(oucfile,'ImageGridAmp','ImageGridTime'); % Time is in sec here
                for t=1:length(sa),
                    
                    % find column of ImageGridAmp corresponding to sa(t)
                    col = find(ImageGridTime == sa(t)/1000);
                    if isempty(col), error('Cannot find time column in ImageGridTime\n'); end
                    
                    p.mesh.data.Cdata{p.mesh.current} = ImageGridAmp(:,col) .* scale;
                    Dif.ouc(:,t) = p.mesh.data.Cdata{p.mesh.current};
                    
                    p.mesh.file = sprintf('rh.sa_%04d_ouc_float',sa(t));
                    mesh_write(p);
                end
                
                clear ImageGridAmp ImageGridTime; pack;
                
                for t=1:length(sa),
                    
                    p.mesh.data.Cdata{p.mesh.current} = Dif.oac(:,t) - Dif.ouc(:,t);
                    p.mesh.file = sprintf('rh.sa_%04d_dif_float',sa(t));
                    mesh_write(p);
                    
                end
                
                clear Dif; pack;
                
                
                
            case 'wm',
                
                %---------------------------------------------------------------------
                % WM condition
                
                vert = size(p.mesh.data.vertices{p.mesh.current},1);
                times = length(wm);
                
                Dif.tac = zeros(vert,times);
                Dif.oac = Dif.tac;
                
                load(tacfile,'ImageGridAmp','ImageGridTime'); % Time is in sec here
                for t=1:length(wm),
                    
                    % find column of ImageGridAmp corresponding to wm(t)
                    col = find(ImageGridTime == wm(t)/1000);
                    if isempty(col), error('Cannot find time column in ImageGridTime\n'); end
                    
                    p.mesh.data.Cdata{p.mesh.current} = ImageGridAmp(:,col) .* scale;
                    Dif.tac(:,t) = p.mesh.data.Cdata{p.mesh.current};
                    
                    p.mesh.file = sprintf('rh.wm_%04d_tac_float',wm(t));
                    mesh_write(p);
                    
                end
                
                clear ImageGridAmp ImageGridTime; pack;
                
                load(oacfile,'ImageGridAmp','ImageGridTime'); % Time is in sec here
                for t=1:length(wm),
                    
                    % find column of ImageGridAmp corresponding to sa(t)
                    col = find(ImageGridTime == wm(t)/1000);
                    if isempty(col), error('Cannot find time column in ImageGridTime\n'); end
                    
                    p.mesh.data.Cdata{p.mesh.current} = ImageGridAmp(:,col) .* scale;
                    Dif.oac(:,t) = p.mesh.data.Cdata{p.mesh.current};
                    
                    p.mesh.file = sprintf('rh.wm_%04d_oac_float',wm(t));
                    mesh_write(p);
                end
                
                clear ImageGridAmp ImageGridTime; pack;
                
                for t=1:length(wm),
                    
                    p.mesh.data.Cdata{p.mesh.current} = Dif.tac(:,t) - Dif.oac(:,t);
                    p.mesh.file = sprintf('rh.wm_%04d_dif_float',wm(t));
                    mesh_write(p);
                    
                end
                
                clear Dif; pack;
                
                
                
            case 'ea',
                
                %---------------------------------------------------------------------
                % EA condition
                
                vert = size(p.mesh.data.vertices{p.mesh.current},1);
                times = length(ea);
                
                Dif.oat = zeros(vert,times);
                Dif.oac = Dif.oat;
                
                load(oatfile,'ImageGridAmp','ImageGridTime'); % Time is in sec here
                for t=1:length(ea),
                    
                    % find column of ImageGridAmp corresponding to sa(t)
                    col = find(ImageGridTime == ea(t)/1000);
                    if isempty(col), error('Cannot find time column in ImageGridTime\n'); end
                    
                    p.mesh.data.Cdata{p.mesh.current} = ImageGridAmp(:,col) .* scale;
                    Dif.oat(:,t) = p.mesh.data.Cdata{p.mesh.current};
                    
                    p.mesh.file = sprintf('rh.ea_%04d_oat_float',ea(t));
                    mesh_write(p);
                    
                end
                
                clear ImageGridAmp ImageGridTime; pack;
                
                load(oacfile,'ImageGridAmp','ImageGridTime'); % Time is in sec here
                for t=1:length(ea),
                    
                    % find column of ImageGridAmp corresponding to sa(t)
                    col = find(ImageGridTime == ea(t)/1000);
                    if isempty(col), error('Cannot find time column in ImageGridTime\n'); end
                    
                    p.mesh.data.Cdata{p.mesh.current} = ImageGridAmp(:,col) .* scale;
                    Dif.oac(:,t) = p.mesh.data.Cdata{p.mesh.current};
                    
                    p.mesh.file = sprintf('rh.ea_%04d_oac_float',ea(t));
                    mesh_write(p);
                end
                
                clear ImageGridAmp ImageGridTime; pack;
                
                for t=1:length(ea),
                    
                    p.mesh.data.Cdata{p.mesh.current} = Dif.oat(:,t) - Dif.oac(:,t);
                    p.mesh.file = sprintf('rh.ea_%04d_dif_float',ea(t));
                    mesh_write(p);
                    
                end
                
                clear Dif; pack;
            end
        end
    end
end
