%  STUDY        - studyset structure containing some or all files in ALLEEG
%  measure      - Measure to display. One from: {'daterp','datspec','icaerp','icaspec'}
%  dataflag     - [1|0] Flag for 'data' or 'limoresults' respectively
%  level        - [1|2] Level of analysis to plot. [1] Single subject
%                 analysis results. [2] Group level analysis results.
%  stats        - Cell array of  {'property','property_value',...} of stats settings.
%                 Properties can be: pvalue ('p'), ('MCC, Number of bootstrap
%                 ('bootstrap'),('tfce'). If not provided, defaults are:
%                 {'p',0.05,'MCC',1,'dir',pwd,'bootstrap',0,'tfce',0};
% plottype      - Type of plot [1|2|3] 
% If Plot type 1 ('plottype', 1), 
% subjindx      - Index on STUDY.datasetinfo.subject of the subject to
%                 plot.
% testindx      - Index (if level 1of the test to plot.
% I.E.
% std_limoresults(STUDY,'flagdata',1,'plottype',1,'measure','daterp','level',1,'subjindx',1,'testindx',1); % To do
% std_limoresults(STUDY,'flagdata',0,'plottype',1,'measure','daterp','level',1,'subjindx',1,'testindx',1);
%
% If Plot type 2 ('plottype', 2)
%
%
% Optional inputs:
%  'measure' - ['daterp'|'icaerp'|'datspec'|'icaspec'|'datersp'|'icaersp']
% Example
%  %Plot 1 (imagesc + topo)
% std_limoresults(STUDY,'plottype',1,'flagdata',1,'measure','daterp','level',1,'subjindx',1,'regressor',{'[1]'});% To do (data single subj level)
% % To do (data group level)
% std_limoresults(STUDY,'plottype',1,'flagdata',0,'measure','daterp','level',1,'subjindx',1,'testindx',1);       % (results single subj level)
% std_limoresults(STUDY,'plottype',1,'flagdata',0,'measure','daterp','level',2,'testindx',1);                    % (results group level)
% 
%  %Plot 2 (topo)
% std_limoresults(STUDY,'plottype',2,'flagdata',1,'measure','daterp','level',1,'subjindx',1,'regressor',{'[1]'});% (data single subj level)
% % To do (data group level)
% std_limoresults(STUDY,'plottype',2,'flagdata',0,'measure','daterp','level',1,'subjindx',1,'testindx',1);       % (results single subj level)
% std_limoresults(STUDY,'plottype',2,'flagdata',0,'measure','daterp','level',2,'testindx',1);                    % (results group level)
% 
% % Plot 3 (Time serie + significance)
% std_limoresults(STUDY,'plottype',3,'flagdata',1,'measure','daterp','level',1,'subjindx',1,'regressor',{'[1]'},'chanindx',6);                % (data single subj level)
% % To do (data group level)
% std_limoresults(STUDY,'plottype',3,'flagdata',0,'measure','daterp','level',1,'subjindx',1,'testindx',1,'regressor',{'[1 2]'},'chanindx',6); % (results single subj level)
% std_limoresults(STUDY,'plottype',3,'flagdata',0,'measure','daterp','level',2,'testindx',1,'chanindx',6);     

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

function std_limoresults(STUDY,varargin)

if nargin < 1
    help std_limoresults;
    return;
end
measure_list = {'daterp', 'datspec', 'icaerp', 'icaspec'};
opt = finputcheck(varargin, { 'plottype'       'integer'   [1 2 3]              []     ;
                              'flagdata'       'integer'   [1,0]                ''     ;
                              'measure'        'string'    measure_list         ''     ;
                              'subjindx'       'integer'   []                   []      ;
                              'level'          'integer'   [1 2]                []      ;
                              'testindx'       'integer'   []                   []      ;
                              'chanindx'       'integer'   []                   []      ;
                              'regressor'      'cell'      []                   {};
                              'stats'          'cell'      { }                  {}}, 'std_limoresults');
if ischar(opt), error(opt); end

% Stats
for i = 1:2:numel(opt.stats)
    handles.(opt.stats{i}) = opt.stats{i+1};
end

try handles.p         ; catch, handles.p         = 0.05; end
try handles.MCC       ; catch, handles.MCC       = 1;    end
try handles.dir       ; catch, handles.dir       = pwd;  end
try handles.bootstrap ; catch, handles.bootstrap = 0;    end
try handles.tfce      ; catch, handles.tfce      = 0;    end

% Limo indx
limostruct_indx = find(~cellfun(@isempty,strfind({STUDY.design(STUDY.currentdesign).limo.datatype},opt.measure)));
 if isempty(limostruct_indx)
     fprintf(2,['std_limoresults error: Can not find measeure ' opt.measure '\n']);
     exit;
 end
 
% --- EEGLAB ADDED ---
% Loading files
if opt.level == 2  
    files                   = {STUDY.design(STUDY.currentdesign).limo(limostruct_indx).groupmodel.filename};
    [PathName,  tmp1, tmp2] = fileparts(files{opt.testindx});
    FileName                = [tmp1 tmp2]; clear tmp1 tmp2;
elseif opt.level == 1
    if opt.flagdata 
        PathName = fileparts(STUDY.design(STUDY.currentdesign).limo(limostruct_indx).model(opt.subjindx).filename{1});
        FileName                = ''; 
    else
        if opt.plottype == 3, opt.testindx = 1; end
            
        [PathName,  tmp1, tmp2] = fileparts(STUDY.design(STUDY.currentdesign).limo(limostruct_indx).model(opt.subjindx).filename{opt.testindx});
        FileName                = [tmp1 tmp2]; clear tmp1 tmp2;
    end
else
    fprintf(2,['std_limoresults error: Not a valid input ''' num2str(opt.level) ''' \n']);
    return;
end
% ---------------------

switch opt.plottype
    % --------------------------------------------------------
    %             IMAGE ALL (COMBINED PLOT)
    % --------------------------------------------------------
    case 1
        if ~(opt.flagdata)
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
                        msgbox('repeated measure ANOVA tfce is not available at this stage, please use the random effect GUI','action not performed','warn')
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
        else
            fprintf(2,'std_limoresults error: Under construction \n');
        end
        
        % --------------------------------------------------------
        %             TOPOPLOT (SCALP MAPS)
        % --------------------------------------------------------
    case 2
        handles.LIMO = load(fullfile(PathName,'LIMO.mat'));
        if opt.flagdata
            limo_display_results(2,'Yr.mat',PathName,handles.p,handles.MCC,handles.LIMO.LIMO,1,'regressor',opt.regressor);
        else
            limo_display_results(2,FileName,PathName,handles.p,handles.MCC,handles.LIMO.LIMO);
        end
        cd(handles.dir);
        
        % --------------------------------------------------------
        %             COURSE PLOT (TIME COURSE)
        % --------------------------------------------------------
    case 3
        % --- EEGLAB ADDED ---
        FileName = 'LIMO.mat';
        
        % Getting Channels
        selected_chan = opt.chanindx; %get(obj.gui_h.listbox_elect2plot,'Value')-1;
        
        % Getting Regressor
        %selected_reg = get(obj.gui_h.popupmenu_modelvar2plot,'Value');
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
        if opt.flagdata && opt.level ==1
            limo_display_results(3,FileName,PathName,handles.p,handles.MCC,handles.LIMO.LIMO,0,'channels', {num2str(selected_chan)}, 'regressor',opt.regressor,'plot3type','Original');
        elseif ~(opt.flagdata) && opt.level == 1
            plot3type_val = questdlg('Plotting ERP','ERP Options','Modelled','Adjusted','Adjusted');
            limo_display_results(3,FileName,PathName,handles.p,handles.MCC,handles.LIMO.LIMO,0,'channels', {num2str(selected_chan)}, 'regressor',opt.regressor,'plot3type',plot3type_val);
        elseif opt.level == 2
            limo_display_results(3,FileName,PathName,handles.p,handles.MCC,handles.LIMO.LIMO,0,'channels', {num2str(selected_chan)});
        end
        cd(handles.dir);
        
        % --------------------------------------------------------
        %             FrequenciesTime/Freq Plane
        % --------------------------------------------------------
    case 4
        fprintf(2,'std_limoresults error: Under construction \n');
end
end
