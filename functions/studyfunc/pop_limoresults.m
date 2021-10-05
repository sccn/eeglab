% function to show limo results - deprecated

% Copyright (C) 2016 Ramon Martinez
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

classdef pop_limoresults < handle
    
    %class properties - access is private so nothing else can access these
    properties (Access = public)
        gui_h;
        study;
        limofiles_path;
        limofiles_level1;
        limofiles_level2;
        regnames;
        catvarnames;
        datorica_indx;
        reg_indx;
        cat_indx;
        string_mplot;
        limostruct_indx
    end
    
    methods
        % =================================================================
        function obj = pop_limoresults(STUDY,analysis)
            
            if ~ismember(analysis,{'dat','ica'}), return; end
            
            % Check if STUDY.design.limo
            %-------------------------------------------------------------- 
            datorica_list = {'dat','ica'};
            obj.datorica_indx = find(strcmp(analysis,datorica_list));
           
            [obj.gui_h,STUDY] = guibuilder(STUDY,analysis);
            obj.string_mplot  = get(obj.gui_h.popupmenu_measure2plot,'String'); 
            obj.gui_h         = guihandles(obj.gui_h.fig);
            obj.study         = STUDY;
             
            level_tmp       = get(obj.gui_h.popupmenu_level         ,'Value');
            m2plot_tmp      = get(obj.gui_h.popupmenu_measure2plot  ,'Value');
            var2plot_indx   = get(obj.gui_h.popupmenu_modelvar2plot ,'Value');
            [var2plot,filespath] = getmeasures2plot(obj.study,level_tmp-1,m2plot_tmp,obj.datorica_indx);
           
            % Detecting Measure and updating limofiles_path and
            %--------------------------------------------------------------
            obj.limofiles_path            = filespath;
            obj.limofiles_level1.index    = var2plot_indx;
            obj.limofiles_level1.guiname  = var2plot.guiname;
            obj.limofiles_level1.pathname = var2plot.pathname;
            % Updating Electrode list
            %--------------------------------------------------------------
            try
                load(fullfile(filespath,'LIMO.mat'),'LIMO');
                electoplot_list = ['All Channels';{LIMO.data.chanlocs.labels}'];
                electoplot_indx = 1;
                set(obj.gui_h.listbox_elect2plot,'String',electoplot_list);
                set(obj.gui_h.listbox_elect2plot,'Value',electoplot_indx);
            catch
            end
            
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
                                if  ~iscell(STUDY.design(STUDY.currentdesign).variable(i).value{1})
                                    tmpregnames{j} = [STUDY.design(STUDY.currentdesign).variable(i).label '_' STUDY.design(STUDY.currentdesign).variable(i).value{j}];
                                else
                                    for jvars = 1:length(STUDY.design(STUDY.currentdesign).variable(i).value{j})
                                        if jvars == 1
                                            jointtmp = STUDY.design(STUDY.currentdesign).variable(i).value{j}{jvars};
                                        else
                                            jointtmp = strcat (jointtmp, '&',  STUDY.design(STUDY.currentdesign).variable(i).value{j}{jvars});
                                        end
                                    end
                                    tmpregnames{j} = [STUDY.design(STUDY.currentdesign).variable(i).label '_' jointtmp{:}];
                                end
                            end
                            [vartype_tmp{1:length(tmpregnames)}] = deal(STUDY.design(STUDY.currentdesign).variable(i).vartype);
                        else
                            for j = 1:length(STUDY.design(STUDY.currentdesign).variable(i).value)
                                if  ~iscell(STUDY.design(STUDY.currentdesign).variable(i).value{1})
                                    tmpregnames{end+1} = [STUDY.design(STUDY.currentdesign).variable(i).label '_' STUDY.design(STUDY.currentdesign).variable(i).value{j}];
                                else
                                    for jvars = 1:length(STUDY.design(STUDY.currentdesign).variable(i).value{j})
                                        if jvars == 1
                                            jointtmp = STUDY.design(STUDY.currentdesign).variable(i).value{j}{jvars};
                                        else
                                            jointtmp = strcat (jointtmp, '&',  STUDY.design(STUDY.currentdesign).variable(i).value{j}{jvars});
                                        end
                                    end
                                    tmpregnames{end+1} = [STUDY.design(STUDY.currentdesign).variable(i).label '_' jointtmp];
                                end
                                
                                %
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
                
                % Set regressors order as LIMO assumes it (first all categorical, then all continuous)
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
                
                 % Rgressors names in case of split continuous ones
                 %---------------------------------------------------------
                if ~isempty(cont_indx) && (Nreg ~= length(tmpregnames))
                    display('Split regressor detected. Giving generic names to continuous regressors ...' );
                    if ~isempty(cat_indx), inintcont = length(cat_indx)+1; else inintcont = 1 ; end
                    
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
                
                % Retrieving list of independent variables
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
            
            %set the callback function for buttom_designadd
            set(obj.gui_h.pushbutton_designadd, 'Callback', @obj.callback_pushbutton_designadd);
            
        end
        % =================================================================
        function obj = callback_popupmenu_level(obj,~,~)
            
            val_level = get(obj.gui_h.popupmenu_level,'Value');        % Getting value (subject from level)
            val_mplot = get(obj.gui_h.popupmenu_measure2plot,'Value'); % Getting value from popupmenu_measure2plot
            val_dor = get(obj.gui_h.popupmenu_dataorresult,'Value');   % Getting value from popupmenu_dataorresult
            
            % Getting values popupmenu_modelvar2plot
            string_vplot  = get(obj.gui_h.popupmenu_modelvar2plot ,'String');
            val_vplot     = get(obj.gui_h.popupmenu_modelvar2plot ,'Value' );
            val_typeplot  = get(obj.gui_h.popupmenu_plottype      ,'Value' );
            string_elec   = get(obj.gui_h.listbox_elect2plot      ,'String');
            val_elec      = get(obj.gui_h.listbox_elect2plot      ,'Value' );
             
            % Single subject level
            if val_level ~= 1
   
                [var2plot,filespath,tmp,obj.study] = getmeasures2plot(obj.study,val_level-1,val_mplot,obj.datorica_indx,1);
                
                if strcmp(get(obj.gui_h.popupmenu_measure2plot,'String'),{' '})
                    set(obj.gui_h.popupmenu_measure2plot  ,'String', obj.string_mplot);
                    set(obj.gui_h.popupmenu_measure2plot  ,'Value',1);
                end
                
                % Updating Electrode list
                % =========================================================
                
                if ~strcmp('No Variables Computed',var2plot.guiname{1})
                    
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
                            if ~ismember(string_vplot(val_vplot),var2plot.guiname)
                                % Getting index of  'Condition_effect_1.mat' for default value
                                % 1 st pass
                                var2plot_indx = find(strcmp(var2plot.guiname,'Condition_effect_1'), 1);
                                % 2nd pass
                                if isempty(var2plot_indx)
                                    var2plot_indx = find(strcmp(var2plot.guiname,'Covariate_effect_1'),1);
                                end
                                % Setting to 1 if all fails
                                if isempty(var2plot_indx)
                                    var2plot_indx = 1;
                                end
                                set(obj.gui_h.popupmenu_modelvar2plot,'String',var2plot.guiname);
                                set(obj.gui_h.popupmenu_modelvar2plot,'Value',var2plot_indx);
                            else
                                var2plot_indx = find(strcmp(string_vplot(val_vplot),var2plot.guiname));
                                set(obj.gui_h.popupmenu_modelvar2plot,'String',var2plot.guiname);
                                set(obj.gui_h.popupmenu_modelvar2plot,'Value',var2plot_indx);
                            end
                         end
                        %--------------------------------------------------
                    elseif val_typeplot == 3
                        %var2plot_list = string_vplot;

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
                    set(obj.gui_h.popupmenu_measure2plot  ,'Enable','on');
                    set(obj.gui_h.popupmenu_dataorresult  ,'Enable','on');
                    set(obj.gui_h.popupmenu_modelvar2plot ,'Enable','on');
                    set(obj.gui_h.popupmenu_plottype      ,'Enable','on');
                    set(obj.gui_h.pushbutton_plot         ,'Enable','on');
                    set(obj.gui_h.checkbox_stats	      ,'Enable','on');
                    set(obj.gui_h.edit_pval               ,'Enable','on');
                else
                    % Enabling off GUI features
                    %------------------------------------------------------
                    set(obj.gui_h.popupmenu_dataorresult  ,'Enable','off');
                    set(obj.gui_h.popupmenu_modelvar2plot ,'Enable','off');
                    set(obj.gui_h.popupmenu_plottype      ,'Enable','off');
                    set(obj.gui_h.pushbutton_plot         ,'Enable','off');
                    
                    eeglab_warning('Make sure to compute the model for this measure');
                end

                % Updating limofiles_level1 structure
                %----------------------------------------------------------
                if strcmp(get(obj.gui_h.popupmenu_modelvar2plot,'Enable'),'off')
                    obj.limofiles_level1.guiname = [];
                else
                    obj.limofiles_path          = filespath;
                    obj.limofiles_level1.guiname  = var2plot.guiname;
                    obj.limofiles_level1.pathname = var2plot.pathname;
                    obj.limofiles_level1.index    = get(obj.gui_h.popupmenu_modelvar2plot,'Value');
                end
                
            % Group level   
            else 
                measure_list  = {'erp','spec'};
                datorica_list = {'dat', 'ica'};
                datatype = [datorica_list{obj.datorica_indx} measure_list{get(obj.gui_h.popupmenu_measure2plot,'Value')}];
                
                obj.limostruct_indx = find(~cellfun(@isempty,strfind({obj.study.design(obj.study.currentdesign).limo.datatype},datatype)));  
                if ~isempty(obj.limostruct_indx) && isfield(obj.study.design(obj.study.currentdesign).limo, 'groupmodel') && ~isempty(obj.study.design(obj.study.currentdesign).limo(obj.limostruct_indx).groupmodel)
                    obj.limofiles_level2.guiname  = {obj.study.design(obj.study.currentdesign).limo(obj.limostruct_indx).groupmodel.guiname}';
                    obj.limofiles_level2.pathname = {obj.study.design(obj.study.currentdesign).limo(obj.limostruct_indx).groupmodel.filename}';
                    obj.limofiles_level2.index = 1;
                    set(obj.gui_h.popupmenu_modelvar2plot,'String',obj.limofiles_level2.guiname);
                    set(obj.gui_h.popupmenu_modelvar2plot,'Value',1);
                    
                    % Loading expected electrodes
                    %----------------------------------------------------------
                    try
                        load(obj.study.design(obj.study.currentdesign).limo(obj.limostruct_indx).chanloc);
                        electoplot_list = ['All Channels';{expected_chanlocs.labels}'];
                        set(obj.gui_h.listbox_elect2plot,'String',electoplot_list);
                    catch
                        % implement this
                    end
                    
                    if val_typeplot == 1 || val_typeplot == 2
                        
                        % Updating electrode list val
                        %----------------------------------------------
                        set(obj.gui_h.listbox_elect2plot,'Value',1);
                        
                    elseif val_typeplot == 3
                        
                        % Updating electrode list val
                        %--------------------------------------------------
                        newval_elec = find(strcmp(string_elec(val_elec),electoplot_list));
                        if isempty(newval_elec)
                            newval_elec = 2;
                        end
                        set(obj.gui_h.listbox_elect2plot,'Value',newval_elec);
                    end
                    
                    % Enabling on/off GUI features
                    %------------------------------------------------------
                    set(obj.gui_h.popupmenu_dataorresult  ,'Value', 2);
                    set(obj.gui_h.popupmenu_dataorresult  ,'Enable','off');
                    set(obj.gui_h.popupmenu_modelvar2plot ,'Enable','on');
                    set(obj.gui_h.popupmenu_plottype      ,'Enable','on');
                    set(obj.gui_h.pushbutton_plot         ,'Enable','on');
                    set(obj.gui_h.popupmenu_mcc	          ,'Enable','off');
                    set(obj.gui_h.popupmenu_compmethod    ,'Enable','off');
                    set(obj.gui_h.checkbox_stats	      ,'Enable','off');
                    set(obj.gui_h.edit_pval               ,'Enable','on');
                    set(obj.gui_h.pushbutton_designadd    ,'Enable','off');% Check this
                else
                    % Disabling GUI features
                    %------------------------------------------------------
                    set(obj.gui_h.popupmenu_dataorresult  ,'Enable','off');
                    set(obj.gui_h.popupmenu_modelvar2plot ,'Enable','off');
                    set(obj.gui_h.popupmenu_plottype      ,'Enable','off');
                    set(obj.gui_h.pushbutton_plot         ,'Enable','off');
                    set(obj.gui_h.pushbutton_designadd    ,'Enable','off');
                    set(obj.gui_h.edit_pval               ,'Enable','off');
                    eeglab_warning('No group level results founded');
                end
            end
        end
        % =================================================================
        function obj = callback_popupmenu_modelvar2plot(obj,~,~)
            
            var2plot_indx = get(obj.gui_h.popupmenu_modelvar2plot ,'Value');
            if get(obj.gui_h.popupmenu_level,'Value') ~= 1
                obj.limofiles_level1.index = var2plot_indx;
            else
                obj.limofiles_level2.index = var2plot_indx;
            end 
        end
        % =================================================================
        function obj = callback_popupmenu_plottype(obj,~,~)
            
            val_typeplot = get( obj.gui_h.popupmenu_plottype    ,'Value');
            val_dor      = get(obj.gui_h.popupmenu_dataorresult ,'Value');
            val_level    = get( obj.gui_h.popupmenu_level       ,'Value');
            
            if (val_typeplot == 1 || val_typeplot == 2) 
                if val_level ~= 1
                    % Call callback_popupmenu_level
                    obj = callback_popupmenu_level(obj);
                end
                
                % Updating listbox_elect2plot and disabling it
                set(obj.gui_h.listbox_elect2plot   ,'Value',1);
                set(obj.gui_h.listbox_elect2plot   ,'Enable','off');
                
                % Disabling 'Add design' button
                if  val_dor == 1
                    set(obj.gui_h.pushbutton_designadd ,'Enable','off');
                else
                    set(obj.gui_h.pushbutton_designadd ,'Enable','on');
                end
                
            elseif val_typeplot == 3
                
                if val_level ~= 1
                    if val_dor == 2
                        listtmp = obj.regnames;
                    elseif val_dor == 1
                        listtmp = obj.catvarnames;
                    end
                    set(obj.gui_h.popupmenu_modelvar2plot ,'String',listtmp);
                    set(obj.gui_h.popupmenu_modelvar2plot ,'Value',1);
                end
                set(obj.gui_h.listbox_elect2plot   ,'Enable','on');
                set(obj.gui_h.listbox_elect2plot   ,'Value' ,2);
                set(obj.gui_h.pushbutton_designadd ,'Enable','off');
            end
            
        end
        % =================================================================
        function obj = callback_popupmenu_dataorresult(obj,~,~)     
            val_ptype = get( obj.gui_h.popupmenu_plottype   ,'Value');
            val_dor   = get(obj.gui_h.popupmenu_dataorresult,'Value');
            
            % Case for data
            if val_dor == 1
                
                if val_ptype == 1 || val_ptype == 2
                    set(obj.gui_h.popupmenu_modelvar2plot,'String',obj.catvarnames);
                    set(obj.gui_h.popupmenu_modelvar2plot,'Value',1);
                elseif val_ptype == 3
                    set(obj.gui_h.popupmenu_modelvar2plot,'String',obj.catvarnames); % THIS MUST INCLUDE CONT!!!!!!!!!!!!!?????
                    set(obj.gui_h.popupmenu_modelvar2plot,'Value',1);
                end
                
                % 'Add design' disable
                set(obj.gui_h.pushbutton_designadd ,'Enable','off');
                
            % Case for results    
            elseif val_dor == 2
                if val_ptype == 1 || val_ptype == 2
                    % Call callback_popupmenu_level
                    obj = callback_popupmenu_level(obj);
                    
                % 'Add design' enable
                set(obj.gui_h.pushbutton_designadd ,'Enable','on');
                    
                elseif val_ptype == 3
                    set(obj.gui_h.popupmenu_modelvar2plot,'String',obj.regnames);
                    set(obj.gui_h.popupmenu_modelvar2plot,'Value',1);
                    
                    % Call callback_popupmenu_level
                    obj = callback_popupmenu_level(obj);
                    
                    % 'Add design' disable
                    set(obj.gui_h.pushbutton_designadd ,'Enable','off');
                    
                end
            end
        end
        % =================================================================
        function obj = callback_checkbox_stats(obj,~,~)
            
            if get(obj.gui_h.checkbox_stats        ,'Value');
                set(obj.gui_h.popupmenu_mcc        ,'Enable','on');
                set(obj.gui_h.popupmenu_compmethod ,'Enable','on');
            else
                set(obj.gui_h.popupmenu_mcc        ,'Enable','off');
                set(obj.gui_h.popupmenu_mcc        ,'Value',1);
                set(obj.gui_h.popupmenu_compmethod ,'Enable','off');
                set(obj.gui_h.popupmenu_compmethod ,'Value',1);
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
            
             % Prepare inputs for std_limoresults
            level       = get(obj.gui_h.popupmenu_level,'Value');
            if level == 1        
                levelval = 2;
                testindxval = obj.limofiles_level2.index;
            else
                subjindxval = level - 1;
                levelval    = 1;
                testindxval = obj.limofiles_level1.index;
            end
            measure_list  = {'erp','spec'};
            datorica_list = {'dat', 'ica'};
            measure_val = [datorica_list{obj.datorica_indx} measure_list{get(obj.gui_h.popupmenu_measure2plot,'Value')}];
            
            % DEFS
            checkbox_stat_status = get(obj.gui_h.checkbox_stats,'Value');
            if checkbox_stat_status
                handles.p = str2double(get(obj.gui_h.edit_pval,'string'));
                indxtmp      = get(obj.gui_h.popupmenu_mcc,'Value');
                handles.MCC  = indxtmp;
                handles.tfce = 0;
                handles.dir  = pwd;
                indxtmp           = get(obj.gui_h.popupmenu_compmethod,'Value');
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
                     
            % Determine if 'Original Data' or 'Results'            
            if get(obj.gui_h.popupmenu_dataorresult,'Value') == 1
                flagdataval = 1;
            else
                flagdataval = 0 ;
            end
                
            % Getting the plot type
            plottypeval   = get(obj.gui_h.popupmenu_plottype,'Value');
            if plottypeval == 3
                chanindxval = get( obj.gui_h.listbox_elect2plot,'Value')-1;
                if chanindxval == 0,
                    eeglab_warning('Must select a valid electrode');
                    return;
                end
            end
            % ---------------------
            % --- Data ---
            if flagdataval == 1
                % Level 1
                if levelval == 1
                    if plottypeval ~= 3
                        std_limoresults(obj.study,'plottype',plottypeval,'flagdata',flagdataval,'measure',measure_val,'level',1,'subjindx',subjindxval,'regressor',obj.cat_indx(get(obj.gui_h.popupmenu_modelvar2plot,'value')));
                    else
                        std_limoresults(obj.study,'plottype',3,'flagdata',flagdataval,'measure',measure_val,'level',1,'subjindx',subjindxval,'regressor',obj.cat_indx(get(obj.gui_h.popupmenu_modelvar2plot,'value')),'chanindx',chanindxval);
                    end
                    % Level 2
                else
                    fprintf('Under construction');
                    
                end
            % --- Results ---
            else
                % Level 1
                if levelval == 1
                    if plottypeval ~= 3
                    std_limoresults(obj.study,'plottype',plottypeval,'flagdata',0,'measure',measure_val,'level',1,'subjindx',subjindxval,'testindx',testindxval);
                    else
                        std_limoresults(obj.study,'plottype',3,'flagdata',0,'measure',measure_val,'level',1,'subjindx',subjindxval,'regressor',obj.reg_indx(get(obj.gui_h.popupmenu_modelvar2plot,'value')),'chanindx',chanindxval);
                    end
                        
                % Level 2
                else 
                    if plottypeval ~= 3
                    std_limoresults(obj.study,'plottype',plottypeval,'flagdata',0,'measure',measure_val,'level',2,'testindx',testindxval);   
                    else
                        std_limoresults(obj.study,'plottype',3,'flagdata',0,'measure',measure_val,'level',2,'testindx',testindxval,'chanindx',chanindxval);  
                    end
                end
            end
        end
        % =================================================================
        function obj = callback_pushbutton_designadd(obj,~,~)
            clear -global LIMO;
            
            val_level  = get( obj.gui_h.popupmenu_level        ,'Value');
            val_mplot  = get(obj.gui_h.popupmenu_measure2plot  ,'Value');
            stringtmp  = get(obj.gui_h.popupmenu_modelvar2plot ,'String');
            pwdtmp = pwd; 
            limo_contrast_manager(fullfile(obj.limofiles_path,'LIMO.mat'));
            cd(pwdtmp);
            htmp = findall(0,'Type','Figure','Tag','figure_limo_contrast_manager');
            if ~isempty(htmp)
                waitfor(htmp); clear htmp;
            end
            % Get new generated file and add it to the list as current
            % ---------------------------------------------------------
            [currentfiles,tmp,tmp,STUDY] = getmeasures2plot(obj.study,val_level-1,val_mplot,obj.datorica_indx,1);
            obj.study = STUDY;
            newfile = setdiff(currentfiles.guiname,stringtmp);
            if ~isempty(newfile)
                %stringtmp(end+1) = newfile;
                set(obj.gui_h.popupmenu_modelvar2plot,'String',currentfiles.guiname);
                set(obj.gui_h.popupmenu_modelvar2plot,'Value',length(currentfiles.guiname));
                %obj.limofiles_filename = newfile{1};
            else
                set(obj.gui_h.popupmenu_modelvar2plot,'Value',1);
            end
        end
    end
end
% =================================================================
%% =================================================================
function [handles,STUDY] = guibuilder(STUDY,analysis)

% Create list "Level of Analysis" (Using First Subject)
%--------------------------------------------------------------------------
subjnames = STUDY.design(STUDY.currentdesign).cases.value;

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
    [var2plot,tmp,tmp,STUDY] = getmeasures2plot(STUDY,1,measure2plot_indx,datorica,1);
    if ~iscell(var2plot)
        var2plot_list = var2plot.guiname;
    else
        var2plot_list = var2plot;
    end
    % Getting index of  'Condition_effect_1.mat' for default value
    var2plot_indx = find(strcmp(var2plot.guiname,'Condition_effect_1'), 1);
    if isempty(var2plot_indx)
        var2plot_indx = find(strcmp(var2plot.guiname,'Covariate_effect_1'),1);
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
textstring = {'Level of Analysis','Data Measure','Plot Type','Data or Results','Variable to Plot','Threshold (p-value )','Electrodes to Plot','Method','Computational Approach'};
% Stuff for popupmenus
%--------------------------------------------------------------------------
plottype_list    = {'Combined Plot','Scalp Maps Series','Time Course'};
datorresult_list = {'Data';'Linear Model Results'};
mcc_list         = {'None';'Clustering';'TFCE';'Max'};
compmethod_list  = {'None','Default Method'};
xalign1          = 0.450; %0.481;
popupnames = {'popupmenu_level','popupmenu_measure2plot','popupmenu_plottype','popupmenu_dataorresult','popupmenu_modelvar2plot','popupmenu_mcc','popupmenu_compmethod'};
popuppos   = {[xalign1 0.921 0.550 0.0625],...           % level
              [xalign1 0.829 0.550 0.0625],...           % measure2plot
              [xalign1 0.737 0.550 0.0625],...           % plottype
              [xalign1 0.645 0.550 0.0625],...           % dataorresult
              [xalign1+0.08 0.553 0.550-0.08 0.0625],... % modelvar2plot
              [0.535 0.561 0.465 0.273 ],...             % mcc
              [0.535 0.152 0.465 0.273]};                % compmethod
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
    'Position',[549.072 185.686 283.538 466.212],...
    'Resize', 'off');

%Panel 1 (Plot Settings)
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

%Panel 2 (Compute Statistics)
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
    'position',[xalign1+0.02 0.451 0.271 0.0724],...
    'string',pval_default,...
    'backgroundcolor',[0.929 0.929 0.929],...
    'tooltipString',ttiptext{6},...
    'Tag','edit_pval');

% List
handles.listbox_elect2plot = uicontrol('parent',handles.panel_1,...
    'style','listbox',...
    'units','normalized',...
    'position',[xalign1+0.02 0.0263 0.488 0.395],...
    'string',electoplot_list,...
    'Value',electoplot_indx,...
    'backgroundcolor',[0.929 0.929 0.929],...
    'FontSize',11,...
    'enable','off',...
    'tooltipString',ttiptext{7},...
    'Tag','listbox_elect2plot');

%Panel 2 (Compute Statistics)
%--------------------------------------------------------------------------
% Checkbox
handles.checkbox_stats = uicontrol('parent',handles.fig,...
    'style','checkbox',...
    'units','normalized',...
    'position',[0.73 0.253 0.106 0.0493],...
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

% Button add design
handles.pushbutton_designadd = uicontrol('parent',handles.panel_1,...
    'style','pushbutton',...
    'units','normalized',...
    'position',[xalign1+0.02 0.553+0.005 0.07 0.0625 - 0.005],...
    'string','+',...
    'FontWeight','bold',...
    'FontSize', 12,...
    'tooltipString','Add contrast',...
    'backgroundcolor',[0.929 0.929 0.929],...
    'Tag', 'pushbutton_designadd');
end

%% =================================================================
function [var2plot,filespath,limoindx,STUDY] = getmeasures2plot(STUDY,subjN,measureindx,datoricaindx,checknewfile_flag)

% measureindx  is the index to {'erp','spec'}
% datoricaindx is the index to {'dat', 'ica'}

% Init
filespath = [];
limoindx = [];
var2plot.guiname{1}  = 'No Variables Computed';
measure_list         = {'erp','spec'};
datorica_list        = {'dat', 'ica'};
if ~exist('checknewfile_flag','var'), checknewfile_flag = 0; end

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
    % Generating list
    %----------------
    filespath         = STUDY.design(STUDY.currentdesign).limo(limoindx).basefolder_level1{subjN};
    
    if checknewfile_flag
        tmp = dir(filespath);
        if ~isempty({tmp.name})
            var2plot_list        = {tmp.name}';
            
            % Cleaning out the list
            %----------------------
            list2cleanout = {'Betas.mat','LIMO.mat','Yr.mat','Yhat.mat','Res.mat','.','..'};
            for i = 1:length(list2cleanout)
                ind2delete{i} = find(strcmp(var2plot_list,list2cleanout{i}));
            end
            nonemptyvals              = find(cellfun(@(x) ~isempty(x), ind2delete));
            ind2delete                = [ind2delete{nonemptyvals}];
            var2plot_list(ind2delete) = [];
            
            %Update STUDY
           
            [newfiles, newfilesindx] = setdiff(cellfun(@(x) fullfile(filespath,x),var2plot_list,'UniformOutput', false),STUDY.design(STUDY.currentdesign).limo(limoindx).model(subjN).filename);
            if ~isempty(newfiles)
                STUDY.design(STUDY.currentdesign).limo(limoindx).model(subjN).guiname(newfilesindx) = cellfun(@(x) x(1:end-4),var2plot_list(newfilesindx),'UniformOutput', false);
                STUDY.design(STUDY.currentdesign).limo(limoindx).model(subjN).filename(newfilesindx) = cellfun(@(x) fullfile(filespath,x),var2plot_list(newfilesindx),'UniformOutput', false);
                STUDY = pop_savestudy( STUDY, [],'filepath', STUDY.filepath,'savemode','resave');
                display('pop_limoresults: New contrast files added to STUDY ');
            end
        end
    end
    
    var2plot.guiname  = (STUDY.design(STUDY.currentdesign).limo(limoindx).model(subjN).guiname);
    var2plot.pathname = (STUDY.design(STUDY.currentdesign).limo(limoindx).model(subjN).filename);
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
