classdef pop_limoresults < handle
    
    %class properties - access is private so nothing else can access these
    properties (Access = private)
        gui_h;
        study;
        limofiles_path;
        limofiles_filename;
        regnames;
        catvarnames;
        datorica_indx;
        reg_indx;
        cat_indx;
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
             
            level_tmp     = get(obj.gui_h.popupmenu_level        ,'Value');
            m2plot_tmp    = get(obj.gui_h.popupmenu_measure2plot ,'Value');
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
                % Getting the Names of the variables
                tmpregnames  = {};
                vartype_tmp  = {};
                reg_indx     = {};
                for i = 1: length(STUDY.design(STUDY.currentdesign).variable)
                    if strcmp(STUDY.design(STUDY.currentdesign).variable(i).vartype,'categorical')
                        if isempty(tmpregnames)
                            for j = 1: length(STUDY.design(STUDY.currentdesign).variable(i).value)
                                tmpregnames{j} = [STUDY.design(STUDY.currentdesign).variable(i).label '_' STUDY.design(STUDY.currentdesign).variable(i).value{j}];
                            end
                            [vartype_tmp{1:length(tmpregnames)}] = deal(STUDY.design(STUDY.currentdesign).variable(i).vartype);
                        else
                            for j = 1:length(STUDY.design(STUDY.currentdesign).variable(i).value)
                                tmpregnames{end+1} = [STUDY.design(STUDY.currentdesign).variable(i).label '_' STUDY.design(STUDY.currentdesign).variable(i).value{j}];
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
                
                %---------------------------------------------------------
                cat_comb = [];
                if ~isempty(cat_indx) && size(cat_indx,2) > 1
                    cat_comb = nchoosek(1:length(cat_indx),2);
                end
                
                % Generating cat index
                %---------------------------------------------------------
                if ~isempty(cat_indx)
                    allcat_indx  = [1: length(cat_indx),(length(tmpregnames) + 1):(length(tmpregnames) + size(cat_comb,1))];
                end
                
                % Adding Regressor to names
                %---------------------------------------------------------
                for i = 1:length(tmpregnames)
                    %tmpregnames{i} = ['Regressor_' tmpregnames{i}];
                    reg_indx{i}    = ['[' num2str(i) ']'];
                end
 
                % Adding Combinations of cat regressors
                %----------------------------------------------------------
                if ~isempty(cat_comb)
                    for j =1:size(cat_comb,1)
                        strtmp = tmpregnames(cat_comb(j,:));
                        tmpregnames{end+1} = [strtmp{1} sprintf(' - %s', strtmp{2:end})];
                        clear strtmp;
                        reg_indx{end+1}    = ['[' num2str(cat_comb(j,:) ) ']'];
                    end
                end
                
                % Retreiving list of independent variables
                if ~isempty(cat_indx)
                    tmpregnames2 = tmpregnames(allcat_indx);
                    cat_indx     =  reg_indx(allcat_indx);
                else
                    tmpregnames2 = [];
                end
                
                % Adding 'Regressor' to names
                %---------------------------------------------------------
                for i = 1:length(tmpregnames)
                    tmpregnames{i} = ['Regressor_' tmpregnames{i}];
                    plusindx = strfind(tmpregnames{i},' - ');
                    if ~isempty(plusindx)
                        tmpregnames{i} = [tmpregnames{i}(1:plusindx+2) 'Regressor_' tmpregnames{i}(plusindx+3:end)];
                    end    
                end
                
                % Updating fields
                %----------------------------------------------------------
                obj.regnames    = tmpregnames;
                obj.catvarnames = tmpregnames2;
                obj.reg_indx    = reg_indx;
                obj.cat_indx    = cat_indx;
                
            catch
                % Case where model is not computed
                %----------------------------------------------------------
                display(['File LIMO.mat not founded in :'  fullfile(obj.limofiles_path,'LIMO.mat')]);
                set(obj.gui_h.popupmenu_modelvar2plot,'Enable','off');
                set(obj.gui_h.popupmenu_plottype     ,'Enable','off');
                set(obj.gui_h.pushbutton_plot        ,'Enable','off');
                eeglab_warning('Make sure to compute the model for this measure');
            end
           
            % CALLBACKS
            %--------------------------------------------------------------
            %set the callback function for button_plot
            set(obj.gui_h.pushbutton_plot, 'Callback', @obj.callback_plot);
            
            %set the callback function for popupmenu_level
            set(obj.gui_h.popupmenu_level,'Callback', @obj.callback_popupmenu_level);
            
            %set the callback function for popupmenu_measure2plot
            set(obj.gui_h.popupmenu_measure2plot,'Callback', @obj.callback_popupmenu_level);
            
            %set the callback function for popupmenu_plottype
            set(obj.gui_h.popupmenu_plottype,'Callback', @obj.callback_popupmenu_plottype);
            
            %set the callback function for popupmenu_dataorresult
            set(obj.gui_h.popupmenu_dataorresult,'Callback', @obj.callback_popupmenu_dataorresult);
            
            %set the callback function for popupmenu_modelvar2plot
            set(obj.gui_h.popupmenu_modelvar2plot,'Callback', @obj.callback_popupmenu_modelvar2plot);
            
            %set the callback function for checkbox_stats
            set(obj.gui_h.checkbox_stats,'Callback', @obj.callback_checkbox_stats);
            
            %set the callback function for checkbox_stats
            set(obj.gui_h.listbox_elect2plot,'Callback', @obj.callback_listbox_elect2plot); 
            
        end
        % =================================================================
        function obj = callback_popupmenu_level(obj,~,~)
            % Getting value (subject from level)
            val_level = get( obj.gui_h.popupmenu_level,'Value');
            if val_level ~= 1
                % Getting value from popupmenu_measure2plot
                % =========================================================
                val_mplot = get(obj.gui_h.popupmenu_measure2plot,'Value');
                
                % Getting value from popupmenu_dataorresult
                % =========================================================
                val_dor = get(obj.gui_h.popupmenu_dataorresult,'Value');
                
                % Getting values popupmenu_modelvar2plot
                % =========================================================
                string_vplot  = get(obj.gui_h.popupmenu_modelvar2plot,'String');
                val_vplot     = get(obj.gui_h.popupmenu_modelvar2plot,'Value');
                val_typeplot  = get(obj.gui_h.popupmenu_plottype     ,'Value'); 
                string_elec   = get(obj.gui_h.listbox_elect2plot     ,'String'); 
                val_elec      = get(obj.gui_h.listbox_elect2plot     ,'Value');
                
                [var2plot_list,filespath] = getmeasures2plot(obj.study,val_level-1,val_mplot,obj.datorica_indx);
                
                % Updating Electrode list
                % =========================================================
                
                if ~strcmp('No Variables Computed',var2plot_list{1})
                    
                    load(fullfile(filespath,'LIMO.mat'),'LIMO');
                    electoplot_list = ['All Channels';{LIMO.data.chanlocs.labels}'];
                    set(obj.gui_h.listbox_elect2plot,'String',electoplot_list);
                    
                    % Getting index for var2plot_list and setting electrode
                    % index
                    if val_typeplot == 1 || val_typeplot == 2
                        
                            % Updating electrode list val
                            %----------------------------------------------
                            set(obj.gui_h.listbox_elect2plot,'Value',1);
                            
                         % Just Checking if the selected variable exist in
                         % the sub (case of Results only)
                         %-------------------------------------------------
                         if val_dor == 2 
                            if ~ismember(string_vplot(val_vplot),var2plot_list)
                                % Getting index of  'Condition_effect_1.mat' for default value
                                % 1 st pass
                                var2plot_indx = find(strcmp(var2plot_list,'Condition_effect_1.mat'), 1);
                                % 2nd pass
                                if isempty(var2plot_indx)
                                    var2plot_indx = find(strcmp(var2plot_list,'Covariate_effect_1.mat'),1);
                                end
                                % Setting to 1 if all fails
                                if isempty(var2plot_indx)
                                    var2plot_indx = 1;
                                end
                                set(obj.gui_h.popupmenu_modelvar2plot,'String',var2plot_list);
                                set(obj.gui_h.popupmenu_modelvar2plot,'Value',var2plot_indx);
                            else
                                var2plot_indx = find(strcmp(string_vplot(val_vplot),var2plot_list));
                                set(obj.gui_h.popupmenu_modelvar2plot,'String',var2plot_list);
                                set(obj.gui_h.popupmenu_modelvar2plot,'Value',var2plot_indx);
                            end
                         end
                        %--------------------------------------------------
                    elseif val_typeplot == 3
                        var2plot_list = string_vplot;

                        % Updating electrode list val
                        %--------------------------------------------------
                        newval_elec = find(strcmp(string_elec(val_elec),electoplot_list));
                        if isempty(newval_elec)
                            newval_elec = 2;
                        end
                        set(obj.gui_h.listbox_elect2plot,'Value',newval_elec);
                    end
                    
                    % Enabling on GUI features
                    %------------------------------------------------------
                    set(obj.gui_h.popupmenu_dataorresult  ,'Enable','on');
                    set(obj.gui_h.popupmenu_modelvar2plot ,'Enable','on');
                    set(obj.gui_h.popupmenu_plottype      ,'Enable','on');
                    set(obj.gui_h.pushbutton_plot         ,'Enable','on');
                else
                    % Enabling off GUI features
                    %------------------------------------------------------
                    set(obj.gui_h.popupmenu_dataorresult  ,'Enable','off');
                    set(obj.gui_h.popupmenu_modelvar2plot ,'Enable','off');
                    set(obj.gui_h.popupmenu_plottype      ,'Enable','off');
                    set(obj.gui_h.pushbutton_plot         ,'Enable','off');
                    
                    eeglab_warning('Make sure to compute the model for this measure');
                end

                % Updating limofiles_path and limofiles_filename
                %----------------------------------------------------------
                obj.limofiles_path     = filespath;
                if strcmp(get(obj.gui_h.popupmenu_modelvar2plot,'Enable'),'off')
                    obj.limofiles_filename = [];
                else
                    obj.limofiles_filename = var2plot_list{get(obj.gui_h.popupmenu_modelvar2plot,'Value')};
                end
                
            else
                % Just to show individual results
                %----------------------------------------------------------
                eeglab_warning('Invalid selection for individual results. Selecting 1st subject instead');
                set( obj.gui_h.popupmenu_level,'Value',2);
                obj = callback_popupmenu_level(obj);
            end
        end
        % =================================================================
        function obj = callback_popupmenu_modelvar2plot(obj,~,~)
            val_level  = get( obj.gui_h.popupmenu_level       ,'Value');
            val_mplot  = get(obj.gui_h.popupmenu_measure2plot ,'Value');
            plottype   = get(obj.gui_h.popupmenu_plottype     ,'Value');
            stringtmp  = get(obj.gui_h.popupmenu_modelvar2plot,'String');
            valtmp     = get(obj.gui_h.popupmenu_modelvar2plot,'Value');
            
            % Just if last line if selected
            %--------------------------------------------------------------
            if valtmp == length(stringtmp) && (plottype ~= 3)
                limo_contrast_manager(fullfile(obj.limofiles_path,'LIMO.mat'));
                htmp = findall(0,'Type','Figure','Tag','figure_limo_contrast_manager');
                if ~isempty(htmp)
                    waitfor(htmp); clear htmp;
                end
                % Get new generated file and add it to the list as current
                % ---------------------------------------------------------
                currentfiles = getmeasures2plot(obj.study,val_level-1,val_mplot,obj.datorica_indx);
                newfile = setdiff(currentfiles,stringtmp);
                if ~isempty(newfile)
                    stringtmp(end+1) = newfile;
                    set(obj.gui_h.popupmenu_modelvar2plot,'String',currentfiles);
                    set(obj.gui_h.popupmenu_modelvar2plot,'Value',length(currentfiles)-1);
                    obj.limofiles_filename = newfile{1};
                else
                    set(obj.gui_h.popupmenu_modelvar2plot,'Value',1);
                end
                
            else
                obj.limofiles_filename = stringtmp{valtmp};
            end
        end
        % =================================================================
        function obj = callback_popupmenu_plottype(obj,~,~)
            val_ptype = get( obj.gui_h.popupmenu_plottype   ,'Value');
            val_dor   = get(obj.gui_h.popupmenu_dataorresult,'Value');
            if val_ptype == 1 || val_ptype == 2
                % Call callback_popupmenu_level
                obj = callback_popupmenu_level(obj);
                
                % Updating listbox_elect2plot
                set(obj.gui_h.listbox_elect2plot,'Value',1);
                set(obj.gui_h.listbox_elect2plot,'Enable','off');
                
            elseif val_ptype == 3
                if val_dor == 2
                    listtmp = obj.regnames;
                elseif val_dor == 1
                    listtmp = obj.catvarnames;
                end
                set(obj.gui_h.popupmenu_modelvar2plot,'String',listtmp);
                set(obj.gui_h.popupmenu_modelvar2plot,'Value',1);
                set(obj.gui_h.listbox_elect2plot     ,'Enable','on');
                set(obj.gui_h.listbox_elect2plot     ,'Value',2);
            end
            
        end
        % =================================================================
        function obj = callback_popupmenu_dataorresult(obj,~,~)     
            val_ptype = get( obj.gui_h.popupmenu_plottype   ,'Value');
            val_dor   = get(obj.gui_h.popupmenu_dataorresult,'Value');
            
            if val_dor == 1
                
                if val_ptype == 1 || val_ptype == 2
                    set(obj.gui_h.popupmenu_modelvar2plot,'String',obj.catvarnames);
                    set(obj.gui_h.popupmenu_modelvar2plot,'Value',1);
                elseif val_ptype == 3
                    set(obj.gui_h.popupmenu_modelvar2plot,'String',obj.catvarnames); % THIS MUST INCLUDE CONT!!!!!!!!!!!!!?????
                    set(obj.gui_h.popupmenu_modelvar2plot,'Value',1);
                end
            elseif val_dor == 2
                if val_ptype == 1 || val_ptype == 2
                    % Call callback_popupmenu_level
                    obj = callback_popupmenu_level(obj);
                    
                elseif val_ptype == 3
                    set(obj.gui_h.popupmenu_modelvar2plot,'String',obj.regnames);
                    set(obj.gui_h.popupmenu_modelvar2plot,'Value',1);
                    
                    % Call callback_popupmenu_level
                    obj = callback_popupmenu_level(obj);
                end 
            end 
        end
        % =================================================================
        function obj = callback_checkbox_stats(obj,~,~)
            
            if get(obj.gui_h.checkbox_stats       ,'Value');
                set(obj.gui_h.popupmenu_mcc       ,'Enable','on');
                set(obj.gui_h.popupmenu_compmethod,'Enable','on');
            else
                set(obj.gui_h.popupmenu_mcc       ,'Enable','off');
                set(obj.gui_h.popupmenu_mcc       ,'Value',1);
                set(obj.gui_h.popupmenu_compmethod,'Enable','off');
                set(obj.gui_h.popupmenu_compmethod,'Value',1);
            end
        end
        % =================================================================
        function obj = callback_listbox_elect2plot(obj,~,~)

            valtmp_plottype = get( obj.gui_h.popupmenu_plottype,'Value');
            valtmp_chan     = get( obj.gui_h.listbox_elect2plot,'Value');            
            
            if (valtmp_plottype == 3) && (valtmp_chan == 1)
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
                
                indxtmp      = get(obj.gui_h.popupmenu_mcc,'Value');
                %valtmp      = get(obj.gui_h.popupmenu_mcc,'string');
                handles.MCC  = indxtmp;
                handles.tfce = 0;
                handles.dir  = pwd;
                indxtmp           = get(obj.gui_h.popupmenu_compmethod,'Value');
                %valtmp            = get(obj.gui_h.popupmenu_mcc,'string');
                if indxtmp == 1
                    handles.bootstrap = 0;
                else
                    handles.bootstrap = 1;
                end
            else
                handles.p         = 0.05;
                handles.MCC       = 1;
                handles.dir       = pwd;
                handles.bootstrap = 0;
                handles.tfce      = 0;
            end
            
            % Loading files
            PathName  = obj.limofiles_path;
            FileName  = obj.limofiles_filename;
            
            % Determine if 'Original Data' or 'Results'
            val_dor   = get(obj.gui_h.popupmenu_dataorresult,'Value');
            
            % Getting the plot type
            ptype_val   = get(obj.gui_h.popupmenu_plottype,'Value');
            % ---------------------

            switch ptype_val
                % --------------------------------------------------------
                %             IMAGE ALL (COMBINED PLOT)
                % --------------------------------------------------------
                case 1                    
                    
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
                    
                % --------------------------------------------------------
                %             TOPOPLOT (SCALP MAPS)
                % --------------------------------------------------------
                case 2
                    handles.LIMO = load(fullfile(PathName,'LIMO.mat')); % eeglab mod
                    limo_display_results(2,FileName,PathName,handles.p,handles.MCC,handles.LIMO.LIMO);
                    cd(handles.dir);

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
                    if val_dor == 1
                        limo_display_results(3,FileName,PathName,handles.p,handles.MCC,handles.LIMO.LIMO,0,'channels', {num2str(selected_chan)}, 'regressor',obj.cat_indx(selected_reg),'plot3type','Original');
                    elseif val_dor == 2
                        plot3type_val = questdlg('Plotting ERP','ERP Options','Modelled','Adjusted','Adjusted');
                        limo_display_results(3,FileName,PathName,handles.p,handles.MCC,handles.LIMO.LIMO,0,'channels', {num2str(selected_chan)}, 'regressor',obj.reg_indx(selected_reg),'plot3type',plot3type_val);
                    end
                    cd(handles.dir);
                    
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

ttiptext{1} = ['Group or individual analysis. Individual analysis' 10 ...          % ttiptext_level
    'is the same a single subject analysis.' 10 ...
    'Group analysis is the analysis of a group of subjects.'];
ttiptext{2} = ['EEG measure to plot.' 10 '(ERP/Log Spectrum)'];                    % ttiptext_m2plot
ttiptext{3} = 'Select to plot Original Data or Results from the Linear Model';     % ttiptext_dorresult
ttiptext{4} = ['Selection of experimental conditions, subject groups, other' 10 ...% ttiptext_mvar2plot
    'regressors and contrast. Some options might be shaded out' 10 ...
    'based on the ?Plot Type?selected below.'];
ttiptext{5} = 'Select type of plot';                                                % ttiptext_plottype
ttiptext{6} = 'Select p-value to compute statistic';                                % ttiptext_pval
ttiptext{7} = 'Select the channels to be used in the plot';                         % ttiptext_chan

ttiptext{8} = 'Select Multiple Comparisons Method';
ttiptext{9} = 'Select Computational Approach to Perform Multiple Comparisons';
% Stuff for texts
%--------------------------------------------------------------------------
textnames = {'text_level','text_measure2plot','text_plottype','text_dataorresult','text_modelvar2plot','text_pval','text_elect2plot','text_mcc','text_compmethod'};
textpos   = {[0.0349 0.941 0.333 0.0428],...
    [0.0349 0.832 0.306 0.0592],...
    [0.0349 0.737 0.186 0.0625],...
    [0.0349 0.655 0.422 0.0526],...
    [0.0349 0.562 0.422 0.0526],...
    [0.0349 0.48 0.399 0.0428 ],...
    [0.0349 0.378 0.345 0.0428],...
    [0.0426 0.591 0.705 0.242 ],...
    [0.0426 0.227 0.481 0.197]};
textstring = {'Level of Analysis','Data Measure','Plot Type','Data or Results','Model Variable to Plot','Threshold (p-value )','Electrodes to Plot','Method','Computational Approach'};
% Stuff for popupmenus
%--------------------------------------------------------------------------
plottype_list    = {'Combined Plot','Scalp Maps','Time Course'};
datorresult_list = {'Original Data';'Linear Model Results'};
mcc_list         = {'None';'Clustering';'TFCE';'Max'};
compmethod_list  = {'None','Default Method'};

popupnames = {'popupmenu_level','popupmenu_measure2plot','popupmenu_plottype','popupmenu_dataorresult','popupmenu_modelvar2plot','popupmenu_mcc','popupmenu_compmethod'};
popuppos   = {[0.481 0.921 0.519 0.0625],... % level
    [0.481 0.829 0.519 0.0625],... % measure2plot
    [0.481 0.737 0.519 0.0625],... % plottype
    [0.481 0.645 0.519 0.0625],... % dataorresult
    [0.481 0.553 0.519 0.0625],... % modelvar2plot
    [0.535 0.561 0.465 0.273 ],... % mcc
    [0.535 0.152 0.465 0.273]};    % compmethod
popupstring = {listlevel,measure2plot_list,plottype_list,datorresult_list,var2plot_list,mcc_list,compmethod_list};
popupindex  = {2,measure2plot_indx,1,2,var2plot_indx,1,1};
%-------------------------------------------------------------------------
compmethod_list = {''};                                                    % UPDATE THIS
pval_default    = num2str(0.5);
icadefs;
color           = BACKEEGLABCOLOR;

% Figure
%--------------------------------------------------------------------------
handles.fig = figure('MenuBar','none',...
    'Name','pop_limoresults',...
    'NumberTitle','off',...
    'Units', 'Points',...
    'Color', color,...
    'Position',[549.072 185.686 273.538 466.212],...
    'Resize', 'off');

%Pannel 1 (Plot Settings)
%--------------------------------------------------------------
% Panel
handles.panel_1 = uipanel('parent',handles.fig,...
    'units','normalized',...
    'position',[0.0182 0.306 0.956 0.683],...
    'title','Plot Settings',...
    'FontSize',12,...
    'backgroundcolor',color,...
    'FontWeight','bold',...
    'Tag', 'panel_1');

%Pannel 2 (Compute Statistics)
%--------------------------------------------------------------------------
handles.panel_2 = uipanel('parent',handles.fig,...
    'units','normalized',...
    'position',[0.0182 0.12 0.956 0.173],...
    'title','Correct for Multiple Comparisons',... % 'Compute Statistics (optional)'
    'backgroundcolor',color,...
    'FontSize',12,...
    'FontWeight','bold',...
    'Tag','panel_2');

% Texts
%--------------------------------------------------------------------------
for i = 1:length(textnames)
    if i > 7 parent = handles.panel_2 ; else parent = handles.panel_1; end
    
    handles.(textnames{i}) = uicontrol('parent'              ,parent ,...
        'style'               ,'text',...
        'units'               ,'normalized',...
        'position'            ,textpos{i},...
        'string'              ,textstring{i},...
        'FontSize'            ,11,...
        'HorizontalAlignment' ,'left',...
        'backgroundcolor'     ,color,...
        'tooltipString'       ,ttiptext{i},...
        'Tag'                 ,textnames{i});
end

% Popups Menus
%--------------------------------------------------------------------------
ttip_popup = {ttiptext{[1:5,8,9]}};
for i = 1:length(popupnames)
    if i > 5 parent = handles.panel_2; else parent = handles.panel_1; end
    
    handles.(popupnames{i}) = uicontrol('parent'        ,parent,...
        'style'         ,'popupmenu',...
        'units'         ,'normalized',...
        'position'      ,popuppos{i},...
        'string'        ,popupstring{i},...
        'value'         ,popupindex{i},...
        'Fontsize'      ,10,...
        'tooltipString' ,ttip_popup{i},...
        'Tag'           ,popupnames{i});
end
set(handles.popupmenu_mcc,       'Enable','off');
set(handles.popupmenu_compmethod,'Enable','off');

% Edit
handles.edit_pval = uicontrol('parent',handles.panel_1,...
    'style','edit',...
    'units','normalized',...
    'position',[0.488 0.451 0.271 0.0724],...
    'string',pval_default,...
    'backgroundcolor',[0.929 0.929 0.929],...
    'tooltipString',ttiptext{6},...
    'Tag','edit_pval');

% List
handles.listbox_elect2plot = uicontrol('parent',handles.panel_1,...
    'style','listbox',...
    'units','normalized',...
    'position',[0.488 0.0263 0.488 0.395],...
    'string',electoplot_list,...
    'Value',electoplot_indx,...
    'backgroundcolor',[0.929 0.929 0.929],...
    'FontSize',11,...
    'enable','off',...
    'tooltipString',ttiptext{7},...
    'Tag','listbox_elect2plot');

%Pannel 2 (Compute Statistics)
%--------------------------------------------------------------------------
% Checkbox
handles.checkbox_stats = uicontrol('parent',handles.fig,...
    'style','checkbox',...
    'units','normalized',...
    'position',[0.755 0.253 0.106 0.0493],...
    'string','',...
    'backgroundcolor',color,...
    'Tag','checkbox_stats');

% Button
handles.pushbutton_plot = uicontrol('parent',handles.fig,...
    'style','pushbutton',...
    'units','normalized',...
    'position',[0.113 0.0321 0.766 0.0664],...
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
measures_computed  = {STUDY.design(STUDY.currentdesign).limo.datatype};

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
        var2plot_list{end+1} = 'Add New Contrast';
        
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