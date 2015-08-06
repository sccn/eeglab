classdef pop_limoresults < handle
    
    %class properties - access is private so nothing else can access these
    properties (Access = private)
        gui_h;
        study;
        limofiles_path;
        limofiles_filename;
        regnames;
        datorica_indx;
        reg_indx;
    end
    
    methods
        % =================================================================
        function obj = pop_limoresults(STUDY,analysis)
            
            if ~ismember(analysis,{'dat','ica'}), return; end
            
            % Check if STUDY.design.limo
            %-------------------------------------------------------------- 
            datorica_list = {'dat','ica'};
            obj.datorica_indx = find(strcmp(analysis,datorica_list));
            
            obj.gui_h    = guibuilder(STUDY,analysis);
            obj.gui_h    = guihandles(obj.gui_h.fig);
            obj.study    = STUDY;
             
            level_tmp     = get(obj.gui_h.popupmenu_level,'Value');
            m2plot_tmp    = get(obj.gui_h.popupmenu_measure2plot,'Value');
            var2plot_indx = get(obj.gui_h.popupmenu_modelvar2plot,'Value');
            [var2plot_list,filespath,limoindx] = getmeasures2plot(obj.study,level_tmp-1,m2plot_tmp,obj.datorica_indx);
           
            % Detecting Measure and updating limofiles_path and
            %--------------------------------------------------------------
            obj.limofiles_path     = filespath;
            obj.limofiles_filename = var2plot_list{var2plot_indx};
            
            % Updating Electrode list
            %--------------------------------------------------------------
            try
                load(fullfile(filespath,'LIMO.mat'),'LIMO');
                electoplot_list = ['All Channels';{LIMO.data.chanlocs.labels}'];
                electoplot_indx = 1;
                set(obj.gui_h.listbox_elect2plot,'String',electoplot_list);
                set(obj.gui_h.listbox_elect2plot,'Value',electoplot_indx);
            catch
            end;
            
            % Creating Name of Regressors
            %--------------------------------------------------------------
            try
                % --
                % Getting the Names of variables
                tmpregnames  = {};
                vartype_tmp  = {};
                reg_indx     = {};
                for i = 1: length(STUDY.design(STUDY.currentdesign).variable)
                    if strcmp(STUDY.design(STUDY.currentdesign).variable(i).vartype,'categorical')
                        if isempty(tmpregnames)
                            tmpregnames                          = STUDY.design(STUDY.currentdesign).variable(i).value;
                            [vartype_tmp{1:length(tmpregnames)}] = deal(STUDY.design(STUDY.currentdesign).variable(i).vartype);
                        else
                            for j = 1:length(STUDY.design(STUDY.currentdesign).variable(i).value)
                                tmpregnames{end+1} = STUDY.design(STUDY.currentdesign).variable(i).value{j};
                                vartype_tmp{end+1} = STUDY.design(STUDY.currentdesign).variable(i).vartype;
                            end
                        end
                    else
                        if isempty(tmpregnames)
                            tmpregnames{1} = STUDY.design(STUDY.currentdesign).variable(i).label;
                            vartype_tmp{1} = STUDY.design(STUDY.currentdesign).variable(i).vartype;

                        else
                            tmpregnames{end+1} = STUDY.design(STUDY.currentdesign).variable(i).label;
                            vartype_tmp{end+1} = STUDY.design(STUDY.currentdesign).variable(i).vartype;
                        end
                    end
                end
                
                % Set regressors order as LIMO asume it (first all categorical, then all continuos)
                %----------------------------------------------------------------------------------
                cont_indx   = find(strcmp('continuous',vartype_tmp));
                cat_indx    = find(strcmp('categorical',vartype_tmp));
                tmpregnames = tmpregnames([cat_indx cont_indx]);
                
                % Adding Baseline
                %----------------------------------------------------------
                tmpregnames{end+1} = 'Baseline';
                % --
                limo = load(fullfile(obj.limofiles_path,'LIMO.mat' ));
                Nreg = size(limo.LIMO.design.X,2);
                
                 % Rgressors names in case of split continous ones
                 %---------------------------------------------------------
                if ~isempty(cont_indx) && (Nreg ~= length(tmpregnames))
                    display('Split regressor detected. Giving generic names to continuos regressors ...' );
                    if ~isempty(cat_indx), inintcont = length(cat_indx)+1; else inintcont = 1 ; end;
                    
                    for i = inintcont:(Nreg-1)
                        tmpregnames{i} = ['Cont_' num2str(i - (inintcont - 1))];
                        display(['Naming Regressor: ' tmpregnames{i}]);
                    end
                end
                
                for i = 1:length(tmpregnames)
                    tmpregnames{i} = ['Reg_' tmpregnames{i}];
                    reg_indx{i}    = ['[' num2str(i) ']'];
                end
 
                % Adding Combinations of cat regressors
                %----------------------------------------------------------
                if ~isempty(cat_indx) && size(cat_indx,2) > 1
                    %for i = 2:size(cat_indx,2)
                    cat_comb = nchoosek(1:length(cat_indx),2);
                    for j =1:size(cat_comb,1)
                        strtmp = tmpregnames(cat_comb(j,:));
                        tmpregnames{end+1} = [strtmp{1} sprintf(' + %s', strtmp{2:end})];
                        clear strtmp;
                        reg_indx{end+1}    = ['[' num2str(cat_comb(j,:) ) ']'];
                    end
                    %end
                end
                
                % Updating fields
                %----------------------------------------------------------
                obj.regnames = tmpregnames;
                obj.reg_indx = reg_indx;
                
            catch
                % Case where model is not computed
                %----------------------------------------------------------
                display(['File LIMO.mat not founded in :'  fullfile(obj.limofiles_path,'LIMO.mat')]);
                set(obj.gui_h.popupmenu_modelvar2plot,'Enable','off');
                set(obj.gui_h.popupmenu_plottype,'Enable','off');
                set(obj.gui_h.pushbutton_plot,'Enable','off');
                eeglab_warning('Make sure to compute the model for this measure');
            end
           
            % CALLBACKS
            %--------------------------------------------------------------
            %set the callback functions for button_plot
            set(obj.gui_h.pushbutton_plot, 'Callback', @obj.callback_plot);
            
            %set the callback functions for popupmenu_level
            set(obj.gui_h.popupmenu_level,'Callback', @obj.callback_popupmenu_level);
            
            %set the callback functions for popupmenu_measure2plot
            set(obj.gui_h.popupmenu_measure2plot,'Callback', @obj.callback_popupmenu_level);
            
            %set the callback functions for popupmenu_modelvar2plot
            set(obj.gui_h.popupmenu_modelvar2plot,'Callback', @obj.callback_popupmenu_modelvar2plot);
            
            %set the callback functions for popupmenu_plottype
            set(obj.gui_h.popupmenu_plottype,'Callback', @obj.callback_popupmenu_plottype);
            
            %set the callback functions for checkbox_stats
            set(obj.gui_h.checkbox_stats,'Callback', @obj.callback_checkbox_stats);
            
            %set the callback functions for checkbox_stats
            set(obj.gui_h.listbox_elect2plot,'Callback', @obj.callback_listbox_elect2plot); 
            
        end
        % =================================================================
        function obj = callback_popupmenu_level(obj,~,~)
            % Getting value (subject from level)
            val_level = get( obj.gui_h.popupmenu_level,'Value');
            if val_level ~= 1
                % Getting value from popupmenu_measure2plot
                % =========================================================
                val_mplot    = get(obj.gui_h.popupmenu_measure2plot,'Value');
                
                % Updating popupmenu_modelvar2plot
                % =========================================================
                string_vplot  = get(obj.gui_h.popupmenu_modelvar2plot,'String');
                val_vplot     = get(obj.gui_h.popupmenu_modelvar2plot,'Value');
                val_typeplot  = get(obj.gui_h.popupmenu_plottype,'Value'); 
                string_elec   = get(obj.gui_h.listbox_elect2plot,'String'); 
                val_elec      = get(obj.gui_h.listbox_elect2plot,'Value');
                
                [var2plot_list,filespath] = getmeasures2plot(obj.study,val_level-1,val_mplot,obj.datorica_indx);
                
                % Updating Electrode list
                % =========================================================
                try
                    load(fullfile(filespath,'LIMO.mat'),'LIMO');
                    electoplot_list = ['All Channels';{LIMO.data.chanlocs.labels}'];
                    set(obj.gui_h.listbox_elect2plot,'String',electoplot_list);
                catch
                end;
                
                if ~strcmp('No Variables Computed',var2plot_list{1})
                    % Getting index for var2plot_list and setting electrode
                    % index
                    if val_typeplot == 1 || val_typeplot == 2
                        % Updating electrode list val
                        set(obj.gui_h.listbox_elect2plot,'Value',1);
                        if ~ismember(string_vplot(val_vplot),var2plot_list)
                            % Getting index of  'Condition_effect_1.mat' for default value
                            % 1 st pass
                            var2plot_indx = find(strcmp(var2plot_list,'Condition_effect_1.mat'), 1);
                            % 2nd pass
                            if isempty(var2plot_indx)
                                var2plot_indx = find(strcmp(var2plot_list,'Covariate_effect_1.mat'),1);
                            end
                            % Setting to 1 all fails
                            if isempty(var2plot_indx)
                                var2plot_indx = 1;
                            end
                            set(obj.gui_h.popupmenu_modelvar2plot,'String',var2plot_list);
                            set(obj.gui_h.popupmenu_modelvar2plot,'Value',var2plot_indx);
                        else
                        end
                    elseif val_typeplot == 3
                        var2plot_list = string_vplot;
                        
                        % Updating electrode list val
                        newval_elec = find(strcmp(string_elec(val_elec),electoplot_list));
                        if isempty(newval_elec)
                            newval_elec = 2;
                        end
                        set(obj.gui_h.listbox_elect2plot,'Value',newval_elec);
                    end
                    set(obj.gui_h.popupmenu_modelvar2plot,'Enable','on');
                    set(obj.gui_h.popupmenu_plottype,'Enable','on');
                    set(obj.gui_h.pushbutton_plot,'Enable','on');
                else
                    set(obj.gui_h.popupmenu_modelvar2plot,'Enable','off');
                    set(obj.gui_h.popupmenu_plottype,'Enable','off');
                    set(obj.gui_h.pushbutton_plot,'Enable','off');
                    eeglab_warning('Make sure to compute the model for this measure');
                end

                % Updating limofiles_path and limofiles_filename
                obj.limofiles_path     = filespath;
                if strcmp(get(obj.gui_h.popupmenu_modelvar2plot,'Enable'),'off')
                    obj.limofiles_filename = [];
                else
                    obj.limofiles_filename = var2plot_list{get(obj.gui_h.popupmenu_modelvar2plot,'Value')};
                end
                
            else
                % Just to show individual results
                eeglab_warning('Invalid selection for individual results. Selecting 1st subject instead');
                set( obj.gui_h.popupmenu_level,'Value',2);
                obj = callback_popupmenu_level(obj);
            end
        end
        % =================================================================
        function obj = callback_popupmenu_modelvar2plot(obj,~,~)
            plottype     = get(obj.gui_h.popupmenu_plottype,'Value');
            stringtmp    = get(obj.gui_h.popupmenu_modelvar2plot,'String');
            valtmp       = get(obj.gui_h.popupmenu_modelvar2plot,'Value');
            
            if valtmp == length(stringtmp) && (plottype ~= 3)
                limo_contrast_manager;
                % Get new generated file and add it to the list as current
                currentfiles = dir(obj.limofiles_path);
                currentfiles = {currentfiles.name}';
                if strcmp(currentfiles{1},'.')
                    currentfiles(1:2) = [];
                end
                stringtmp(end) = [];
                newfile = setdiff(currentfiles,stringtmp);
                 if ~isempty(newfile)
                     stringtmp(end+1) = newfile;
                     stringtmp(end+1) = {'Add New Var'};
                 end    
                 set(obj.gui_h.popupmenu_modelvar2plot,'String',stringtmp);
                 set(obj.gui_h.popupmenu_modelvar2plot,'Value',length(stringtmp)-1); 
                 obj.limofiles_filename = stringtmp{length(stringtmp)-1};
            else
                obj.limofiles_filename = stringtmp{valtmp};
            end
        end
        % =================================================================
        function obj = callback_popupmenu_plottype(obj,~,~)
            val_ptype = get( obj.gui_h.popupmenu_plottype,'Value');
            
            if val_ptype == 1 || val_ptype == 2
                obj = callback_popupmenu_level(obj);
                
                % Updating listbox_elect2plot
                set(obj.gui_h.listbox_elect2plot,'Value',1);
                set(obj.gui_h.listbox_elect2plot,'Enable','off');
            elseif val_ptype == 3
                set(obj.gui_h.popupmenu_modelvar2plot,'String',obj.regnames);
                set(obj.gui_h.popupmenu_modelvar2plot,'Value',1);
                set(obj.gui_h.listbox_elect2plot,'Enable','on');
                set(obj.gui_h.listbox_elect2plot,'Value',2);
            end
            
        end
        % =================================================================
        function obj = callback_checkbox_stats(obj,~,~)
            
            if get(obj.gui_h.checkbox_stats,'Value');
                set(obj.gui_h.popupmenu_mcc,'Enable','on');
                set(obj.gui_h.popupmenu_compmethod,'Enable','on');
            else
                set(obj.gui_h.popupmenu_mcc,'Enable','off');
                set(obj.gui_h.popupmenu_mcc,'Value',1);
                set(obj.gui_h.popupmenu_compmethod,'Enable','off');
                set(obj.gui_h.popupmenu_compmethod,'Value',1);
            end
        end
        % =================================================================
        function obj = callback_listbox_elect2plot(obj,~,~)

            valtmp_plottype = get( obj.gui_h.popupmenu_plottype,'Value');
            valtmp_chan     = get( obj.gui_h.listbox_elect2plot,'Value');            
            
            if valtmp_plottype == 3 && valtmp_chan == 1
                eeglab_warning('Invalid selection for requested plot. Selecting 1st electrode instead ');
                set(obj.gui_h.listbox_elect2plot,'Value',2);                
            end
        end
        % =================================================================
        function obj = callback_plot(obj,~, ~)
            
            % --- EEGLAB ADDED ---
            % DEFS
            checkbox_stat_status = get(obj.gui_h.checkbox_stats,'Value');
            if checkbox_stat_status
                handles.p = str2double(get(obj.gui_h.edit_pval,'string'));
                
                indxtmp     = get(obj.gui_h.popupmenu_mcc,'Value');
                %valtmp      = get(obj.gui_h.popupmenu_mcc,'string');
                handles.MCC = indxtmp;
                
                handles.tfce = 0;
                handles.dir = pwd;
                
                indxtmp           = get(obj.gui_h.popupmenu_compmethod,'Value');
                %valtmp            = get(obj.gui_h.popupmenu_mcc,'string');
                if indxtmp == 1
                    handles.bootstrap = 0;
                else
                    handles.bootstrap = 1;
                end
            else
                handles.p = 0.05;
                handles.MCC = 1;
                handles.dir = pwd;
                handles.bootstrap = 0;
                handles.tfce = 0;
            end
            
            % [FileName,PathName,FilterIndex]=uigetfile('*.mat','Select Univariate Results to display'); % eeglab commented out
            
            PathName = obj.limofiles_path;
            FileName = obj.limofiles_filename;
            
            [tmp,ext] = fileparts(FileName); % NOTE: GET EXTENSIOn!!!!!!!!
            FilterIndex = 1; 

            % ---------------------
            
            % Getting the plot type
            tmpval   = get(obj.gui_h.popupmenu_plottype,'Value');
            
            switch tmpval
                % --------------------------------------------------------
                %             IMAGE ALL (COMBINED PLOT)
                % --------------------------------------------------------
                case 1
                                        
                    % [FileName,PathName,FilterIndex]=uigetfile('*.mat','Select Univariate Results to display');
                    
                    if FilterIndex == 1
                        % cd(PathName); % eeglab commented out
                        handles.LIMO = load(fullfile(PathName,'LIMO.mat'));% eeglab mod
                        
                        % check if bootstrap or tfce should be computed
                        % ---------------------------------------------
                        % 1st level
                        if handles.LIMO.LIMO.Level == 1;
                            if handles.bootstrap == 1 && ~exist(sprintf('H0%sH0_%s', filesep, FileName), 'file') ...
                                    && strncmp(FileName,'con',3) == 0 && strncmp(FileName,'ess',3) ==0
                                if strcmp(questdlg('Level 1: compute all bootstraps?','bootstrap turned on','Yes','No','No'),'Yes');
                                    LIMO = handles.LIMO.LIMO;
                                    LIMO.design.bootstrap = 1;
                                    if handles.tfce == 1
                                        LIMO.design.tfce = 1;
                                    end
                                    save LIMO LIMO
                                    if strcmp(handles.LIMO.LIMO.Analysis,'Time-Frequency')
                                        limo_eeg_tf(4);
                                    else
                                        limo_eeg(4);
                                    end
                                end
                            end
                            
                            if handles.tfce == 1 && ~exist(sprintf('TFCE%stfce_%s', filesep, FileName), 'file') ...
                                    && exist(sprintf('H0%sH0_%s', filesep, FileName), 'file') && strncmp(FileName,'con',3) == 0 ...
                                    && strncmp(FileName,'ess',3) ==0
                                if strcmp(questdlg('Level 1: compute all tfce?','tfce turned on','Yes','No','No'),'Yes');
                                    LIMO = handles.LIMO.LIMO;
                                    LIMO.design.tfce = 1;
                                    save LIMO LIMO
                                    if strcmp(handles.LIMO.LIMO.Analysis,'Time-Frequency')
                                        limo_eeg_tf(4);
                                    else
                                        limo_eeg(4);
                                    end
                                end
                            end
                        end
                        
                        % contrasts stuff
                        if handles.bootstrap == 1 && ~exist(sprintf('H0%sH0_%s', filesep, FileName), 'file')
                            if strncmp(FileName,'con',3)
                                load Yr; cd H0; H0_Betas.mat;
                                result = limo_contrast(Yr, H0_Betas, handles.LIMO.LIMO, 0,3); clear Yr H0_Betas
                            elseif strncmp(FileName,'ess',3)
                                load Yr; cd H0; H0_Betas.mat;
                                result = limo_contrast(Yr, H0_Betas, handles.LIMO.LIMO, 1,3); clear Yr H0_Betas
                            end
                        end
                        
                        if handles.tfce == 1 && ~exist(sprintf('TFCE%stfce_%s', filesep, FileName), 'file') ...
                                && exist(sprintf('H0%sH0_%s', filesep, FileName), 'file')
                            if strncmp(FileName,'con',3)
                                load(FileName);
                                if strcmp(handles.LIMO.LIMO.Analysis,'Time-Frequency'); x = 3;
                                else [x,y,z] = size(con); if x~=1; x=2; end
                                end
                                tfce_score = limo_tfce(x,squeeze(con(:,:,2)),handles.LIMO.LIMO.data.neighbouring_matrix);
                                cd TFCE; filename2 = sprintf('tfce_%s',FileName); save ([filename2], 'tfce_score'); clear con tfce_score
                                cd ..; cd H0; filename = sprintf('H0_%s',FileName); load(filename);
                                tfce_H0_score = limo_tfce(x,squeeze(H0_ess(:,:,2,:)),handles.LIMO.LIMO.data.neighbouring_matrix);
                                filename2 = sprintf('tfce_%s',filename); save ([filename2], 'tfce_H0_score'); clear H0_con tfce_score
                            elseif strncmp(FileName,'ess',3)
                                load(FileName);
                                if strcmp(handles.LIMO.LIMO.Analysis,'Time-Frequency'); x = 3;
                                else [x,y,z] = size(ess); if x~=1; x=2; end
                                end
                                tfce_score = limo_tfce(x,squeeze(ess(:,:,2)),handles.LIMO.LIMO.data.neighbouring_matrix);
                                cd TFCE; filename2 = sprintf('tfce_%s',FileName); save ([filename2], 'tfce_score'); clear ess tfce_score
                                cd ..; cd H0; filename = sprintf('H0_%s',FileName); load(filename);
                                tfce_H0_score = limo_tfce(x,squeeze(H0_ess(:,:,2,:)),handles.LIMO.LIMO.data.neighbouring_matrix);
                                filename2 = sprintf('tfce_%s',filename); save ([filename2], 'tfce_H0_score'); clear H0_ess tfce_score
                            end
                        end
                        
                        % 2nd level
                        nboot = 1000;
                        if handles.LIMO.LIMO.Level == 2;
                            if handles.bootstrap == 1 && ~exist(sprintf('H0%sH0_%s', filesep, FileName), 'file')
                                if strncmp(FileName,'one_sample',10)
                                    load Yr; limo_random_robust(1,Yr,eval(FileName(28:end-4)),nboot,handles.tfce); clear Yr;
                                    LIMO.design.bootstrap = 1; save LIMO LIMO
                                elseif strncmp(FileName,'two_samples',11)
                                    load Y1r; load Y2r; limo_random_robust(2,Y1r,Y2r,eval(FileName(29:end-4)),nboot,handles.tfce); clear Y1r Y2r;
                                    LIMO.design.bootstrap = 1; save LIMO LIMO
                                elseif strncmp(FileName,'paired_samples',14)
                                    load Y1r; load Y2r; limo_random_robust(3,Y1r,Y2r,eval(FileName(32:end-4)),nboot,handles.tfce); clear Y1r Y2r;
                                    LIMO.design.bootstrap = 1; save LIMO LIMO
                                elseif strncmp(FileName,'Repeated_measures',17)
                                    warndlg2('repeated measure ANOVA bootstrap is not availbale at this stage, please use the random effect GUI','action not performed')
                                else
                                    if strcmp(handles.LIMO.LIMO.Analysis,'Time-Frequency')
                                        limo_eeg_tf(4);
                                    else
                                        limo_eeg(4);
                                    end
                                end
                            end
                            
                            if handles.tfce == 1 && ~exist(sprintf('TFCE%stfce_%s', filesep, FileName), 'file') ...
                                    && exist(sprintf('H0%sH0_%s', filesep, FileName), 'file')
                                mkdir tfce; load(FileName); load(sprintf('H0%sH0_%s', filesep, FileName));
                                if strncmp(FileName,'one_sample',10)
                                    parameter = eval(FileName(28:end-4));
                                    tfce_name = sprintf('tfce_one_sample_ttest_parameter_%g',parameter);
                                    tfce_H0_name = sprintf('tfce_H0_one_sample_ttest_parameter_%g',parameter);
                                    if strcmp(handles.LIMO.LIMO.Analysis,'Time-Frequency');
                                        x = size(one_sample,1);
                                        if x==1
                                            x=2; LIMO.LIMO.data.neighbouring_matrix = [];
                                        else
                                            x=3;
                                        end
                                        tfce_one_sample = limo_tfce(x,squeeze(one_sample(:,:,:,4)),handles.LIMO.LIMO.data.neighbouring_matrix);
                                        save(['tfce', filesep, tfce_name], 'tfce_one_sample'); clear tfce_one_sample;
                                        tfce_H0_one_sample = limo_tfce(x,squeeze(H0_one_sample(:,:,:,1,:)),handles.LIMO.LIMO.data.neighbouring_matrix);
                                        save(['H0', filesep, tfce_H0_name],'tfce_H0_one_sample'); clear tfce_H0_one_sample;
                                    else
                                        x = size(one_sample,1); if x~=1; x=2; end
                                        tfce_one_sample = limo_tfce(x,squeeze(one_sample(:,:,4)),handles.LIMO.LIMO.data.neighbouring_matrix);
                                        save(['tfce', filesep, tfce_name], 'tfce_one_sample'); clear tfce_one_sample;
                                        tfce_H0_one_sample = limo_tfce(x,squeeze(H0_one_sample(:,:,1,:)),handles.LIMO.LIMO.data.neighbouring_matrix);
                                        save(['H0', filesep, tfce_H0_name],'tfce_H0_one_sample'); clear tfce_H0_one_sample;
                                    end
                                elseif strncmp(FileName,'two_samples',11)
                                    parameter = eval(FileName(29:end-4));
                                    tfce_name = sprintf('tfce_two_samples_ttest_parameter_%g',parameter);
                                    tfce_H0_name = sprintf('tfce_H0_two_samples_ttest_parameter_%g',parameter);
                                    if strcmp(handles.LIMO.LIMO.Analysis,'Time-Frequency');
                                        x = size(one_sample,1);
                                        if x==1
                                            x=2; LIMO.LIMO.data.neighbouring_matrix = [];
                                        else
                                            x=3;
                                        end
                                        tfce_two_samples = limo_tfce(x,squeeze(two_samples(:,:,:,4)),handles.LIMO.LIMO.data.neighbouring_matrix);
                                        save(['tfce', filesep, tfce_name], 'tfce_two_samples'); clear tfce_two_samples;
                                        tfce_H0_two_samples = limo_tfce(x,squeeze(H0_two_samples(:,:,:,1,:)),handles.LIMO.LIMO.data.neighbouring_matrix);
                                        save(['H0', filesep, tfce_H0_name],'tfce_H0_two_samples'); clear tfce_H0_two_samples;
                                    else
                                        x = size(two_samples,1); if x~=1; x=2; end
                                        tfce_two_samples = limo_tfce(x,squeeze(two_samples(:,:,4)),handles.LIMO.LIMO.data.neighbouring_matrix);
                                        save(['tfce', filesep, tfce_name], 'tfce_two_samples'); clear tfce_two_samples;
                                        tfce_H0_two_samples = limo_tfce(x,squeeze(H0_two_samples(:,:,1,:)),handles.LIMO.LIMO.data.neighbouring_matrix);
                                        save(['H0', filesep, tfce_H0_name],'tfce_H0_two_samples'); clear tfce_H0_two_samples;
                                    end
                                elseif strncmp(FileName,'paired_samples',14)
                                    parameter = eval(FileName(32:end-4));
                                    tfce_name = sprintf('tfce_paired_samples_ttest_parameter_%g',parameter);
                                    tfce_H0_name = sprintf('tfce_H0_paired_samples_ttest_parameter_%g',parameter);
                                    if strcmp(handles.LIMO.LIMO.Analysis,'Time-Frequency')
                                        x = size(one_sample,1);
                                        if x==1
                                            x=2; LIMO.LIMO.data.neighbouring_matrix = [];
                                        else
                                            x=3;
                                        end
                                        tfce_paired_samples = limo_tfce(x,squeeze(paired_samples(:,:,:,4)),handles.LIMO.LIMO.data.neighbouring_matrix);
                                        save(['tfce', filesep, tfce_name], 'tfce_paired_samples'); clear tfce_paired_samples;
                                        tfce_H0_paired_samples = limo_tfce(x,squeeze(H0_paired_samples(:,:,:,1,:)),handles.LIMO.LIMO.data.neighbouring_matrix);
                                        save(['H0', filesep, tfce_H0_name],'tfce_H0_paired_samples'); clear tfce_H0_paired_samples;
                                    else
                                        x = size(paired_samples,1); if x~=1; x=2; end
                                        tfce_paired_samples = limo_tfce(x,squeeze(paired_samples(:,:,4)),handles.LIMO.LIMO.data.neighbouring_matrix);
                                        save(['tfce', filesep, tfce_name], 'tfce_paired_samples'); clear tfce_paired_samples;
                                        tfce_H0_paired_samples = limo_tfce(x,squeeze(H0_paired_samples(:,:,1,:)),handles.LIMO.LIMO.data.neighbouring_matrix);
                                        save(['H0', filesep, tfce_H0_name],'tfce_H0_paired_samples'); clear tfce_H0_paired_samples;
                                    end
                                elseif strncmp(FileName,'Repeated_measures',17)
                                    msgbox('repeated measure ANOVA tfce is not availbale at this stage, please use the random effect GUI','action not performed','warn')
                                else
                                    if strcmp(handles.LIMO.LIMO.Analysis,'Time-Frequency')
                                        limo_eeg_tf(4);
                                    else
                                        limo_eeg(4);
                                    end
                                end
                            end
                        end
                        
                        % do the figure
                        % -------------
                        limo_display_results(1,FileName,PathName,handles.p,handles.MCC,handles.LIMO.LIMO);
                        cd(handles.dir);
                    end
                    % --------------------------------------------------------
                    %             TOPOPLOT (SCALP MAPS)
                    % --------------------------------------------------------
                case 2
                    %[FileName,PathName,FilterIndex]=uigetfile('*.mat','Select Result to plot');
                    if FilterIndex == 1
                        handles.LIMO = load(fullfile(PathName,'LIMO.mat')); % eeglab mod
                        limo_display_results(2,FileName,PathName,handles.p,handles.MCC,handles.LIMO.LIMO);
                        cd(handles.dir);
                    end
                    
                % --------------------------------------------------------
                %             COURSE PLOT (TIME COURSE)
                % --------------------------------------------------------    
                case 3
                    % --- EEGLAB ADDED ---
                    FileName = 'LIMO.mat';
                    
                    % Getting Channels
                    selected_chan = get(obj.gui_h.listbox_elect2plot,'Value')-1;
                    
                    % Getting Regressor
                    selected_reg = get(obj.gui_h.popupmenu_modelvar2plot,'Value');
                    %---------------------
                    if FilterIndex == 1
                        cd(PathName);
                        try
                            handles.LIMO = load('LIMO.mat');
                            if strncmp(handles.LIMO.LIMO.design.name,'one sample',9)
                                files = dir('*.mat');
                                for i=1:length(files)
                                    if strncmp(files(i).name,'one_sample',10)
                                        FileName = files(i).name;
                                    end
                                end
                            elseif strncmp(handles.LIMO.LIMO.design.name,'two samples',9)
                                files = dir('*.mat');
                                for i=1:length(files)
                                    if strncmp(files(i).name,'two_samples',11)
                                        FileName = files(i).name;
                                    end
                                end
                            elseif strncmp(handles.LIMO.LIMO.design.name,'paired t-test',12)
                                files = dir('*.mat');
                                for i=1:length(files)
                                    if strncmp(files(i).name,'paired_sample',13)
                                        FileName = files(i).name;
                                    end
                                end
                            elseif strncmp(handles.LIMO.LIMO.design.name,'Repeated measures ANOVA',22)
                                if handles.LIMO.LIMO.design.nb_conditions == 1 &&  length(handles.LIMO.LIMO.design.repeated_measure) == 1
                                    FileName = 'Rep_ANOVA_Factor_1.mat';
                                else
                                    [FileName,PathName,FilterIndex]=uigetfile('*.mat','Which Effect to plot?');
                                end
                            end
                            
                        catch
                            LIMO = []; handles.LIMO = LIMO;
                        end
                        limo_display_results(3,FileName,PathName,handles.p,handles.MCC,handles.LIMO.LIMO,0,'channels', {num2str(selected_chan)}, 'regressor',obj.reg_indx(selected_reg));
                        cd(handles.dir);
                    end
                % --------------------------------------------------------
                %             FrequenciesTime/Freq Plane
                % --------------------------------------------------------     
                case 4
                    fprintf(2,'Still working on this......\n');
            end
        end
    end
end
% =================================================================
% =================================================================
function handles = guibuilder(STUDY,analysis)

% Create list "Level of Analysis" (Using First Subject)
%--------------------------------------------------------------------------
subjnames = STUDY.design(STUDY.design_index).cases.value;

listlevel = {'Group (Level 2)'};
for i=1:length(subjnames)
    listlevel{i+1,1} = [subjnames{i} ' (Level 1)'];
end

% Create list for "Measures to Plot" (Using First Subject and first measure available)
%--------------------------------------------------------------------------
measure2plot_list     = {'ERP','Power Spectrum'};
datorica_list         = {'dat', 'ica'};
measure2plot_internal = {'erp','spec'};
datorica              =  find(strcmp(analysis,datorica_list));

for i = 1:length(measure2plot_internal)
    measure_flag{i} = any(cellfun(@(x) ~isempty(x),strfind({STUDY.design(STUDY.currentdesign).limo.datatype},measure2plot_internal{i})));
end
measure_flag = cell2mat(measure_flag);

% Create list "Variable to Plot" (Using First Subject)
%--------------------------------------------------------------------------
if ~logical(sum(measure_flag))
    var2plot_list = {'No Variable has been found'};
    var2plot_indx = 1;
else
    tmp = find(measure_flag);
    measure2plot_indx = tmp(1); % With this we are forcing to pick the first measure available
    [var2plot_list] = getmeasures2plot(STUDY,1,measure2plot_indx,datorica);
    
    % Getting index of  'Condition_effect_1.mat' for default value
    var2plot_indx = find(strcmp(var2plot_list,'Condition_effect_1.mat'), 1);
    if isempty(var2plot_indx)
        var2plot_indx = find(strcmp(var2plot_list,'Covariate_effect_1.mat'),1);
    end
end

% Create list elect2plot_list
%--------------------------------------------------------------------------
tmpchan         =  struct2cell(STUDY.changrp(1,:));
electoplot_list = ['All Channels';{tmpchan{1,1,:}}'];
electoplot_indx = 1;

% TOOLTIPS texts
%--------------------------------------------------------------------------

ttiptext_level     = ['Group or individual analysis. Individual analysis' 10 ... 
                      'is the same a single subject analysis.' 10 ...
                      'Group analysis is the analysis of a group of subjects.'];
ttiptext_m2plot    = ['EEG measure to plot.' 10 '(ERP/Log Spectrum)'];
ttiptext_mvar2plot = ['Selection of experimental conditions, subject groups, other' 10 ...
                      'regressors and contrast. Some options might be shaded out' 10 ... 
                      'based on the ?Plot Type?selected below.'];
ttiptext_plottype   = 'Select type of plot';
ttiptext_pval       = 'Select p-value to compute statistic';
ttiptext_chan       = 'Select the channels to be used in the plot';
%--------------------------------------------------------------------------
plottype_list   = {'Combined Plot','Scalp Maps','Time Course'};
compmethod_list = {''};                                                    % UPDATE THIS
pval_default    = num2str(0.5);
icadefs;
color           = BACKEEGLABCOLOR;  

handles.fig = figure('MenuBar','none',...
    'Name','pop_limoresults',...
    'NumberTitle','off',...
    'Units', 'Points',...
    'Color', color,...
    'Position',[549.134 201.682 273.569 450.29],...
    'Resize', 'off');
%Pannel 1 (Plot Settings)
%--------------------------------------------------------------
% Panel
handles.panel_1 = uipanel('parent',handles.fig,...
    'units','normalized',...
    'position',[0.0182 0.361 0.956 0.627],...
    'title','Plot Settings',...
    'FontSize',12,...
    'backgroundcolor',color,...
    'FontWeight','bold',...
    'Tag', 'panel_1');

% Texts
handles.text_level        = uicontrol('parent',handles.panel_1,...
    'style','text',...
    'units','normalized',...
    'position',[0.0271 0.929 0.333 0.0485],...
    'string','Level of Analysis',...
    'FontSize',11,...
    'backgroundcolor',color,...
    'tooltipString',ttiptext_level,...
    'Tag','text_level');

handles.text_measure2plot  = uicontrol('parent',handles.panel_1,...
    'style','text',...
    'units','normalized',...
    'position',[0.0271 0.81 0.306 0.0672],...
    'string','Measure to Plot',...
    'FontSize',11,...
    'backgroundcolor',color,...
    'tooltipString',ttiptext_m2plot,...
    'Tag','text_measure2plot');

handles.text_modelvar2plot = uicontrol('parent',handles.panel_1,...
    'style','text',...
    'units','normalized',...
    'position',[0.0271 0.716 0.422 0.0597],...
    'string','Model Variable to Plot',...
    'FontSize',11,...
    'backgroundcolor',color,...
    'tooltipString',ttiptext_mvar2plot,...
    'Tag','text_modelvar2plot');

handles.text_plottype = uicontrol('parent',handles.panel_1,...
    'style','text',...
    'units','normalized',...
    'position',[0.0271 0.604 0.186 0.0709],...
    'string','Plot Type',...
    'FontSize',11,...
    'backgroundcolor',color,...
    'tooltipString',ttiptext_mvar2plot,...
    'Tag','text_plottype');

handles.text_pval = uicontrol('parent',handles.panel_1,...
    'style','text',...
    'units','normalized',...
    'position',[0.0271 0.53 0.399 0.049],...
    'string','Threshold (p-value )',...
    'Fontsize',11,...
    'backgroundcolor',color,...
    'tooltipString',ttiptext_pval,...
    'Tag','text_pval');

handles.text_elect2plot = uicontrol('parent',handles.panel_1,...
    'style','text',...
    'units','normalized',...
    'position',[0.0271 0.425 0.345 0.0485],...
    'string','Electrodes to Plot',...
    'FontSize',11,...
    'backgroundcolor',color,...
    'tooltipString',ttiptext_chan,...
    'Tag','text_elect2plot');

% Popups Menus
handles.popupmenu_level         = uicontrol('parent',handles.panel_1,...
    'style','popupmenu',...
    'units','normalized',...
    'position',[0.481 0.907 0.519 0.0709],...
    'string',listlevel,...
    'value',2,...
    'Fontsize',10,...
    'tooltipString',ttiptext_level,...
    'Tag', 'popupmenu_level');

handles.popupmenu_measure2plot  = uicontrol('parent',handles.panel_1,...
    'style','popupmenu',...
    'units','normalized',...
    'position',[0.481 0.806 0.519 0.0709],...
    'string',measure2plot_list,...
    'value', measure2plot_indx, ...
    'Fontsize',10,...
    'tooltipString',ttiptext_m2plot,...
    'Tag', 'popupmenu_measure2plot');

handles.popupmenu_modelvar2plot = uicontrol('parent',handles.panel_1,...
    'style','popupmenu',...
    'units','normalized',...
    'position',[0.481 0.705 0.519 0.0709],...
    'string',var2plot_list,...
    'value', var2plot_indx,...
    'Fontsize',10,...
    'tooltipString',ttiptext_mvar2plot,...
    'Tag','popupmenu_modelvar2plot');

handles.popupmenu_plottype      = uicontrol('parent',handles.panel_1,...
    'style','popupmenu',...
    'units','normalized',...
    'position',[0.481 0.604 0.519 0.0709],...
    'string',plottype_list,...
    'Fontsize',10,...
    'tooltipString',ttiptext_plottype,...
    'Tag','popupmenu_plottype');

% Edit
handles.edit_pval = uicontrol('parent',handles.panel_1,...
    'style','edit',...
    'units','normalized',...
    'position',[0.492 0.496 0.271 0.082],...
    'string',pval_default,...
    'backgroundcolor',[0.929 0.929 0.929],...
    'tooltipString',ttiptext_pval,...
    'Tag','edit_pval');

% List
handles.listbox_elect2plot = uicontrol('parent',handles.panel_1,...
    'style','listbox',...
    'units','normalized',...
    'position',[0.492 0.0261 0.488 0.448],...
    'string',electoplot_list,...
    'Value',electoplot_indx,...
    'backgroundcolor',[0.929 0.929 0.929],...
    'FontSize',11,...
    'enable','off',...
    'tooltipString',ttiptext_chan,...
    'Tag','listbox_elect2plot');

%Pannel 2 (Compute Statistics)
%--------------------------------------------------------------------------
handles.panel_2 = uipanel('parent',handles.fig,...
    'units','normalized',...
    'position',[0.0182 0.163 0.957 0.18],...
    'title','Correct for Multiple Comparisons',... % 'Compute Statistics (optional)'
    'backgroundcolor',color,...
    'FontSize',12,...
    'FontWeight','bold',...
    'Tag','panel_2');

% Texts
handles.text_mcc = uicontrol('parent',handles.panel_2,...
    'style','text',...
    'units','normalized',...
    'position',[0.0426 0.455 0.705 0.47],...
    'string',{'Correct for Multiple';'Comparisons'},...
    'HorizontalAlignment','left',...
    'Fontsize',11,...
    'backgroundcolor',color,...
    'Tag','text_mcc');

handles.text_compmethod = uicontrol('parent',handles.panel_2,...
    'style','text',...
    'units','normalized',...
    'position',[0.0426 0.227 0.438 0.197],...
    'string','Computational Method',...
    'Fontsize',11,...
    'backgroundcolor',color,...
    'Tag','text_compmethod');

% Popups Menus
handles.popupmenu_mcc = uicontrol('parent',handles.panel_2,...
    'style','popupmenu',...
    'units','normalized',...
    'position',[0.535 0.545 0.465 0.273],...
    'string',{'None';'Clustering';'TFCE';'Max'},...
    'Fontsize',11,...
    'Value',1,...
    'backgroundcolor',[0.94 0.94 0.94],...
    'enable','off',...
    'Tag','popupmenu_mcc');

handles.popupmenu_compmethod = uicontrol('parent',handles.panel_2,...
    'style','popupmenu',...
    'units','normalized',...
    'position',[0.535 0.152 0.465 0.273],...
    'string',{'None','Default Method'},...
    'Fontsize',10,...
    'backgroundcolor',[0.94 0.94 0.94],...
    'enable','off',...
    'Tag','popupmenu_compmethod');

% Checkbox
handles.checkbox_stats = uicontrol('parent',handles.fig,...
    'style','checkbox',...
    'units','normalized',...
    'position',[0.76 0.305 0.106 0.051],...
    'string','',...
    'backgroundcolor',color,...
    'Tag','checkbox_stats');

% Button
handles.pushbutton_plot = uicontrol('parent',handles.fig,...
    'style','pushbutton',...
    'units','normalized',...
    'position',[0.113 0.051 0.766 0.0687],...
    'string','PLOT',...
    'backgroundcolor',[0.929 0.929 0.929],...
    'Tag', 'pushbutton_plot');
end

% =================================================================
function [var2plot_list,filespath,limoindx] = getmeasures2plot(STUDY,subjN,measureindx,datoricaindx)

% measureindx  is the index to {'erp','spec'}
% datoricaindx is the index to {'dat', 'ica'}

% Init
filespath = '';
var2plot_list = {'No Variables Computed'};
measure_list  = {'erp','spec'};
datorica_list = {'dat', 'ica'};

requested_datatype = [datorica_list{datoricaindx} measure_list{measureindx}];
measures_computed     = {STUDY.design(STUDY.currentdesign).limo.datatype};

for i = 1:length(measure_list)
    measure_flag{i} = any(cellfun(@(x) ~isempty(x),strfind(measures_computed,measure_list{i})));
end
measure_flag = cell2mat(measure_flag);

% Warning for gui
if ~logical(sum(measure_flag))
    eeglab_warning('No results have been found for this subject');
end
if ismember(requested_datatype,measures_computed)
    limoindx  = find(strcmp(requested_datatype,measures_computed));
    filespath = STUDY.design(STUDY.currentdesign).limo(limoindx).foldername{subjN};
    
    % Generating list
    %----------------
    tmp = dir(filespath);
    if ~isempty({tmp.name})
        var2plot_list        = {tmp.name}';
        var2plot_list{end+1} = 'Add New Variable';
        
        % Cleaning out the list
        %----------------------
        list2cleanout = {'Betas.mat','LIMO.mat','Yr.mat','Yhat.mat','Res.mat','.','..'};
        for i = 1:length(list2cleanout)
            ind2delete{i} = find(strcmp(var2plot_list,list2cleanout{i}));
        end
        nonemptyvals              = find(cellfun(@(x) ~isempty(x), ind2delete));
        ind2delete                = [ind2delete{nonemptyvals}];
        var2plot_list(ind2delete) = [];
    end
end
end

% =================================================================
function eeglab_warning(text)

icadefs;
eeglabcolor = BACKEEGLABCOLOR;
h           = warndlg(text,'EEGLAB Warning');
set(h,'color', eeglabcolor);

uiwait(h);
end