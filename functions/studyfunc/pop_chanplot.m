% pop_chanplot() - graphic user interface (GUI)-based function with plotting 
%                options for visualizing. Only channel measures (e.g., spectra, 
%                ERPs, ERSPs, ITCs) that have been computed and saved in the study EEG 
%                datasets can be visualized. These can be computed using the GUI-based 
%                pop_precomp().
% Usage:    
%                >> STUDY = pop_chanplot(STUDY, ALLEEG);   
% Inputs:
%   ALLEEG     - Top-level EEGLAB vector of loaded EEG structures for the dataset(s) 
%                in the STUDY. ALLEEG for a STUDY set is typically loaded using 
%                pop_loadstudy(), or in creating a new STUDY, using pop_createstudy().  
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG
%   datasets in ALLEEG.
%
% Outputs:
%   STUDY      - The input STUDY set structure modified according to specified user edits,
%                if any. Plotted channel measure means (maps, ERSPs, etc.) are added to 
%                the STUDY structure after they are first plotted to allow quick replotting.  
%
% Graphic interface buttons:
%  "Select channel to plot" - [list box] Displays available channels to plot (format is
%                'channel name (number of channels)'). The presented channels depend s
%                on the optional input variable 'channels'. Selecting (clicking on) a 
%                channel from the list will display the selected channel channels in the 
%                "Select channel(s) to plot" list box. Use the plotting buttons below 
%                to plot selected measures of the selected channel.
%  "Select channel(s) to plot" - [list box] Displays the ICA channels of the currently 
%                selected channel (in the "Select channel to plot" list box). Each channel 
%                has the format: 'subject name, channel index'.
%  "Plot channel properties" - [button] Displays in one figure all the mean channel measures
%                (e.g., dipole locations, scalp maps, spectra, etc.) that were calculated
%                and saved in the EEG datsets. If there is more than one condition, the ERP 
%                and the spectrum will have different colors for each condition. The ERSP 
%                and ITC plots will show only the first condition; clicking on the subplot 
%                will open a new figure with the different conditions displayed together. 
%                Uses the command line function std_propplot().
%  "Plot ERSPs" - [button] Displays the channel channel ERSPs. 
%                If applied to a channel, channel ERSPs are plotted in one figure  
%                (per condition) with the channel mean ERSP. If "All # channel centroids" 
%                option is selected, plots all average ERSPs of the channels in one figure 
%                per condition. If applied to channels, display the ERSP images of specified 
%                channel channels in separate figures, using one figure for all conditions.
%                Uses the command line functions std_erspplot().
%  "Plot ITCs" - [button] Same as  "Plot ERSPs" but with ITC.
%                Uses the command line functions std_itcplot().
%  "Plot spectra" - [button] Displays the channel channel spectra.   
%                If applied to a channel, displays channel spectra plus the average channel 
%                spectrum in bold. For a specific channel, displays the channel channel 
%                spectra plus the average channel spectrum (in bold) in one figure per condition.
%                If the "All # channel centroids" option is selected, displays the average 
%                spectrum of all channels in the same figure, with spectrum for different 
%                conditions (if any) plotted in different colors.  
%                If applied to channels, displays the spectrum of specified channel 
%                channels in separate figures using one figure for all conditions.  
%                Uses the command line functions std_specplot().
%  "Plot ERPs" - [button] Same as "Plot spectra" but for ERPs.
%                Uses the command line functions std_erpplot().
%  "Plot ERPimage" - [button] Same as "Plot ERP" but for ERPimave.
%                Uses the command line functions std_erpimplot().
%
% Authors: Arnaud Delorme, Scott Makeig, SCCN/INC/UCSD, October 11, 2004

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, October 11, 2004, arno@sccn.ucsd.edu
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

% Coding notes: Useful information on functions and global variables used.

function [STUDY, com] = pop_chanplot(varargin)

icadefs;
com = [];
if ~ischar(varargin{1})
    if nargin < 2
        error('pop_chanplot(): You must provide ALLEEG and STUDY structures');
    end
    STUDY  = varargin{1};
    STUDY.etc.erpparams.topotime    = []; % [] for channels and NaN for components
    STUDY.etc.specparams.topofreq   = []; % NaN -> GUI disabled
    STUDY.etc.erspparams.topotime   = [];
    STUDY.etc.erspparams.topofreq   = [];
    STUDY.etc.erpimparams.topotime  = [];
    STUDY.etc.erpimparams.topotrial = [];
    
    % test path
    % ---------
    pathwarn = 'off';
    if ~isempty(STUDY.filename)
        if ~strcmpi(pwd, STUDY.filepath) && ~strcmpi(pwd, STUDY.filepath(1:end-1))
            if length(STUDY.datasetinfo(1).filepath) < 1
                pathwarn = 'on';
            elseif STUDY.datasetinfo(1).filepath(1) == '.'
                pathwarn = 'on';
            end
        end
        if isempty(STUDY.filepath) && exist(STUDY.datasetinfo(1).filename) == 2
            pathwarn = 'off';
        end
        if strcmpi(pathwarn, 'on')
            warndlg2(strvcat('You have changed your working path and data files are', ...
                             'no longer available; Cancel, and go back to your STUDY folder'), 'warning');
        end
    end
    
    STUDY.tmphist = '';
    ALLEEG = varargin{2};
    if ~isfield(STUDY, 'changrp')
        STUDY = std_changroup(STUDY, ALLEEG);
        disp('Warning: history not saved for group creation');
    elseif isempty(STUDY.changrp)
        STUDY = std_changroup(STUDY, ALLEEG);
        disp('Warning: history not saved for group creation');
    end
    
    show_chan          = ['pop_chanplot(''showchan'',gcf);'];
    show_onechan       = ['pop_chanplot(''showchanlist'',gcf);'];
	plot_chan_maps     = ['pop_chanplot(''topoplot'',gcf); ']; 
    plot_onechan_maps  = ['pop_chanplot(''plotchantopo'',gcf); ']; 
    plot_chan_ersps    = ['pop_chanplot(''erspplot'',gcf); '];
    plot_onechan_ersps = ['pop_chanplot(''plotchanersp'',gcf); '];
    plot_chan_itcs     = ['pop_chanplot(''itcplot'',gcf); '];
    plot_onechan_itcs  = ['pop_chanplot(''plotchanitc'',gcf); '];
    plot_chan_erpim    = ['pop_chanplot(''erpimageplot'',gcf); '];
    plot_onechan_erpim = ['pop_chanplot(''plotchanerpimage'',gcf); '];
    plot_chan_spectra  = ['pop_chanplot(''specplot'',gcf); '];
    plot_onechan_spectra = ['pop_chanplot(''plotchanspec'',gcf); '];
    plot_chan_erp      = ['pop_chanplot(''erpplot'',gcf); '];
    plot_onechan_erp   = ['pop_chanplot(''plotchanerp'',gcf); '];
    plot_chan_dip      = ['pop_chanplot(''dipplot'',gcf); '];
    plot_onechan_dip   = ['pop_chanplot(''plotchandip'',gcf); '];
    plot_chan_sum      = ['pop_chanplot(''plotsum'',gcf); '];
    plot_onechan_sum   = ['pop_chanplot(''plotonechanum'',gcf); '];
    rename_chan        = ['pop_chanplot(''renamechan'',gcf);']; 
    move_onechan       = ['pop_chanplot(''movecomp'',gcf);'];
    move_outlier       = ['pop_chanplot(''moveoutlier'',gcf);'];
    create_chan        = ['pop_chanplot(''createchan'',gcf);'];
    reject_outliers    = ['pop_chanplot(''rejectoutliers'',gcf);'];
    merge_channels     = ['pop_chanplot(''mergechannels'',gcf);'];
    erp_opt            = ['pop_chanplot(''erp_opt'',gcf);'];
    spec_opt           = ['pop_chanplot(''spec_opt'',gcf);'];
    erpim_opt          = ['pop_chanplot(''erpim_opt'',gcf);'];
    ersp_opt           = ['pop_chanplot(''ersp_opt'',gcf);'];
    stat_opt           = ['pop_chanplot(''stat_opt'',gcf);'];
    create_group       = ['pop_chanplot(''create_group'',gcf);'];
    edit_group         = ['pop_chanplot(''edit_group'',gcf);'];
    delete_group       = ['pop_chanplot(''delete_group'',gcf);'];
    saveSTUDY          = [ 'set(findobj(''parent'', gcbf, ''userdata'', ''save''), ''enable'', fastif(get(gcbo, ''value'')==1, ''on'', ''off''));' ];
    browsesave         = [ '[filename, filepath] = uiputfile2(''*.study'', ''Save STUDY with .study extension -- pop_chan()''); ' ... 
                           'set(faindobj(''parent'', gcbf, ''tag'', ''studyfile''), ''string'', [filepath filename]);' ];
    sel_all_chans      = ['pop_chanplot(''sel_all_chans'',gcf);'];
                       
    % list of channel groups
    % ----------------------
    show_options = {};
    for index = 1:length(STUDY.changrp)
        show_options{end+1} = [ 'All ' STUDY.changrp(index).name ];
    end
    
    % enable buttons
    % --------------
    filename = fullfile(STUDY.datasetinfo(1).filepath, STUDY.datasetinfo(1).subject);
    if exist([filename '.datspec'] ) || exist([filename '_ses-01.datspec']), spec_enable = 'on'; else  spec_enable  = 'off'; end
    if exist([filename '.daterp']  ) || exist([filename '_ses-01.daterp'])  , erp_enable = 'on'; else   erp_enable  = 'off'; end
    if exist([filename '.dattimef']) || exist([filename '_ses-01.dattimef']) ,ersp_enable = 'on'; else  ersp_enable  = 'off'; end
    if exist([filename '.dattimef']) || exist([filename '_ses-01.dattimef'])  ,itc_enable = 'on'; else   itc_enable  = 'off'; end
    if exist([filename '.daterpim']) || exist([filename '_ses-01.daterpim']),erpim_enable = 'on'; else erpim_enable  = 'off'; end
    
    if isfield(ALLEEG(1).dipfit, 'model'), dip_enable   = 'on'; else dip_enable   = 'off'; end
    
    % userdata below
    % --------------
    fig_arg{1}{1} = ALLEEG;
    fig_arg{1}{2} = STUDY;
    fig_arg{1}{3} = STUDY.changrp;
    fig_arg{1}{4} = { STUDY.changrp.name };
    fig_arg{2}    = length(STUDY.changrp);
        
    std_line = [0.9 0.35 0.9];
    geometry = { [0.8 3] [1] [0.6 0.35 0.1 0.1 0.9] std_line std_line std_line std_line std_line std_line  };
    str_name       = sprintf('STUDY name ''%s'' - ''%s''', STUDY.name, STUDY.design(STUDY.currentdesign).name);
    if length(str_name) > 80, str_name = [ str_name(1:80) '...''' ]; end
             
    % list of designs
    uilist   = { ...
        {'style' 'text'       'string' 'Select design:' 'FontWeight' 'Bold' 'HorizontalAlignment' 'center'} ...
        {'style' 'popupmenu'  'string' { STUDY.design.name } 'FontWeight' 'Bold' 'tag' 'design' 'value' STUDY.currentdesign } ...
        { } ...
        {'style' 'text'       'string' 'Select channel to plot' 'FontWeight' 'Bold' } ...
        {'style' 'pushbutton' 'string' 'Sel. all' 'callback' sel_all_chans } {} {} ...
        {'style' 'text'       'string' 'Select subject(s) to plot' 'FontWeight' 'Bold'} ...
        {'style' 'listbox'    'string' show_options 'value' 1 'max' 2 'tag' 'chan_list' 'Callback' show_chan } ...
        {'style' 'pushbutton' 'enable' 'on'         'string'  [ 'STATS' 10 'params' ]  'callback' stat_opt } ...
        {'style' 'listbox'    'string' '' 'tag' 'chan_onechan' 'max' 2 'min' 1  'callback' show_onechan } ... 
        {'style' 'pushbutton' 'enable'   erp_enable 'string' 'Plot ERPs'        'Callback' plot_chan_erp} ...
        {'style' 'pushbutton' 'enable'   erp_enable 'string' 'Params'           'Callback' erp_opt }  ...
        {'style' 'pushbutton' 'enable'   erp_enable 'string' 'Plot ERP(s)'      'Callback' plot_onechan_erp} ...
        {'style' 'pushbutton' 'enable'  spec_enable 'string' 'Plot spectra'     'Callback' plot_chan_spectra} ...
        {'style' 'pushbutton' 'enable'  spec_enable 'string' 'Params'           'Callback' spec_opt }  ...
        {'style' 'pushbutton' 'enable'  spec_enable 'string' 'Plot spectra'     'Callback' plot_onechan_spectra} ...
        {'style' 'pushbutton' 'enable' erpim_enable 'string' 'Plot ERPimage'    'Callback' plot_chan_erpim } ...
        {'style' 'pushbutton' 'enable' erpim_enable 'string' 'Params'           'Callback' erpim_opt } ... 
        {'style' 'pushbutton' 'enable' erpim_enable 'string' 'Plot ERPimage(s)' 'Callback' plot_onechan_erpim } ...        
        {'style' 'pushbutton' 'enable'  ersp_enable 'string' 'Plot ERSPs'       'Callback' plot_chan_ersps} ...
        {'vertexpand' 2.15 'style' 'pushbutton' 'enable'  ersp_enable 'string' 'Params' 'Callback' ersp_opt }  ...
        {'style' 'pushbutton' 'enable'  ersp_enable 'string' 'Plot ERSP(s)'     'Callback' plot_onechan_ersps}...
        {'style' 'pushbutton' 'enable'   itc_enable 'string' 'Plot ITCs'        'Callback' plot_chan_itcs} { }  ...
        {'style' 'pushbutton' 'enable'   itc_enable 'string' 'Plot ITC(s)'      'Callback' plot_onechan_itcs}...
        };
    %    {'style' 'pushbutton' 'string' 'Plot channel properties' 'Callback' plot_chan_sum} {} ... 
    %{'style' 'pushbutton' 'string' 'Plot channel properties (soon)' 'Callback' plot_onechan_sum 'enable' 'off'}
    
    % additional UI given on the command line
    % ---------------------------------------
    geomvert = [ 1 0.5 1 5 1 1 1 1 1];
    if nargin > 2
        addui = varargin{3};
        if ~isfield(addui, 'uilist')
            error('Additional GUI definition (argument 4) requires the field "uilist"');
        end
        if ~isfield(addui, 'geometry')
            addui.geometry = mat2cell(ones(1,length(addui.uilist)));
        end
        uilist = { uilist{:}, addui.uilist{:} };
        geometry = { geometry{:} addui.geometry{:} };
        geomvert = [ geomvert ones(1,length(addui.geometry)) ];
    end
    [out_param userdat] = inputgui( 'geometry' , geometry, 'uilist', uilist, ...
        'helpcom', 'pophelp(''pop_chanplot'')', ...
        'title', 'View and edit current channels -- pop_chanplot()' , 'userdata', fig_arg, ...
        'geomvert', geomvert, 'eval', show_chan );
	
   if ~isempty(userdat)
       ALLEEG = userdat{1}{1};
       STUDY  = userdat{1}{2};
   end

   % history
   % -------
   com = STUDY.tmphist;
   STUDY = rmfield(STUDY, 'tmphist');
   
else
    hdl = varargin{2};  %figure handle
    userdat  = get(varargin{2}, 'userdat');    
    ALLEEG   = userdat{1}{1};
    STUDY    = userdat{1}{2};
    cls      = userdat{1}{3};
    allchans = userdat{1}{4};
    
    design  = get(findobj('parent', hdl, 'tag', 'design')      , 'value');
    changrp = get(findobj('parent', hdl, 'tag', 'chan_list')   , 'value');
    onechan = get(findobj('parent', hdl, 'tag', 'chan_onechan'), 'value');
	if ~std_checkdesign(STUDY, design)
        return;
    end
   
    try
        switch  varargin{1}

            case {'topoplot', 'erspplot','itcplot','specplot', 'erpplot', 'erpimageplot' }
                changrpstr = allchans(changrp);
                plotting_option = varargin{1};
                plotting_option = [ plotting_option(1:end-4) 'plot' ];
                a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''channels'','  vararg2str({changrpstr}) ', ''design'', ' int2str(design) ');' ];
                 % update Study history
                eval(a); STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, a);  
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat);

            case {'plotchanersp','plotchanitc','plotchanspec', 'plotchanerp','plotchanerpimage' }
                changrpstr    = allchans(changrp);

                %if length(changrp) > 1
                %    subject = STUDY.subject{onechan-1};
                %else
                %    changrpstruct = STUDY.changrp(changrp);
                %    allsubjects   = unique_bc({ STUDY.datasetinfo([ changrpstruct.setinds{:} ]).subject });
                %    subject = allsubjects{onechan-1};
                %end

                plotting_option = varargin{1};
                plotting_option = [ plotting_option(9:end) 'plot' ];
                if onechan(1) ~= 1  % check that not all onechan in channel are requested
                     subject = STUDY.design(STUDY.currentdesign).cases.value{onechan-1};
                     a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''channels'','  vararg2str({changrpstr}) ', ''subject'', ''' subject ''', ''design'', ' int2str(design) ' );' ];
                     eval(a); STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, a);  
                 else
                    a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''channels'','  vararg2str({changrpstr}) ', ''plotsubjects'', ''on'', ''design'', ' int2str(design) ' );' ];
                    eval(a); STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, a);
                 end
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat); 

            case 'stat_opt' % save the list of selected channels
                [STUDY com] = pop_statparams(STUDY);
                if ~isempty(com)
                    STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, com);
                end
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat); %update information (STUDY)

            case 'erp_opt' % save the list of selected channels
                [STUDY com] = pop_erpparams(STUDY);
                if ~isempty(com)
                    STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, com);
                end
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat); %update information (STUDY)     

            case 'spec_opt' % save the list of selected channels
                [STUDY com] = pop_specparams(STUDY);
                if ~isempty(com)
                    STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, com);
                end
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat); %update information (STUDY)     

            case 'ersp_opt' % save the list of selected channels
                [STUDY com] = pop_erspparams(STUDY);
                if ~isempty(com)
                    STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, com);
                end
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat); %update information (STUDY)     

            case 'erpim_opt' % save the list of selected channels
                [STUDY com] = pop_erpimparams(STUDY);
                if ~isempty(com)
                    STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, com);
                end
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat); %update information (STUDY)     

            case 'showchanlist' % save the list of selected channels
                if length(changrp) == 1
                    STUDY.changrp(changrp).selected = onechan;
                end
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat); %update information (STUDY)     

           case 'showchan'
                cind     = get(findobj('parent', hdl, 'tag', 'chan_list')   , 'value');
                changrp  = STUDY.changrp(cind);

                % Find datasets availaible
                % ------------------------
                %setind = STUDY.setind .* (changrp.chaninds > 0); % set to 0 the cell not
                %%                                       % containing any electrode
                %allchansets = unique_bc( setind(find(setind(:))) );

                % Generate channel list
                % ---------------------
                chanid{1} = 'All subjects';
                if length(changrp) == 1
                    allsubjects = STUDY.design(STUDY.currentdesign).cases.value;
                    for l = 1:length(allsubjects)
                        chanid{end+1} = [ allsubjects{l} ' ' changrp.name ];
                    end
                else
                    for l = 1:length(STUDY.design(STUDY.currentdesign).cases.value)
                        chanid{end+1} = [ STUDY.design(STUDY.currentdesign).cases.value{l} ];
                    end
                end

                selected = 1;
                if isfield(changrp, 'selected') && length(cind) == 1
                    if ~isempty(STUDY.changrp(cind).selected)
                        selected = min(STUDY.changrp(cind).selected, 1+length(chanid));
                        STUDY.changrp(cind).selected = selected;
                    end
                end

                set(findobj('parent', hdl, 'tag', 'chan_onechan'), 'value', selected, 'String', chanid);

            case 'sel_all_chans'
                set(findobj('parent', hdl, 'tag', 'chan_list'), 'value', [1:length(STUDY.changrp)]);

                % Generate channel list
                % ---------------------
                chanid{1} = 'All subjects';
                for l = 1:length(STUDY.design(STUDY.currentdesign).cases.value)
                    chanid{end+1} = [ STUDY.design(STUDY.currentdesign).cases.value{l} ' All' ];
                end
                selected = 1;
                set(findobj('parent', hdl, 'tag', 'chan_onechan'), 'value', selected, 'String', chanid);

            case 'plotsum'
                changrpstr = allchans(changrp);
                [STUDY] = std_propplot(STUDY, ALLEEG, allchans(changrp));
                a = ['STUDY = std_propplot(STUDY, ALLEEG, ' vararg2str({ allchans(changrp) }) ' );'  ];
                STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, a);  
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat);    

            case 'create_group'
                channames = { STUDY.changrp(changrp).name };
                for i=1:length(channames), channames{i} = [ ' ' channames{i} ]; end
                channamestr = strcat(channames{:});
                res = inputdlg2({ 'Name of channel group', 'Channels to group' }, 'Create channel group', 1, { '' channamestr(2:end) });
                if isempty(res), return; end
                STUDY.changrp(end+1).name = res{1};
                allchans(end+1)         = { res{1} };
                chanlabels = parsetxt(res{2});
                if length(chanlabels) == 1
                    warndlg2('Cannot create a channel group with a single channel');
                    return;
                end
                STUDY.changrp(end).channels = chanlabels;
                tmp = std_chanlookup( STUDY, ALLEEG, STUDY.changrp(end));
                STUDY.changrp(end).chaninds = tmp.chaninds;
                userdat{1}{2} = STUDY;
                userdat{1}{4} = allchans;
                set(hdl, 'userdat',userdat);    

                % list of channel groups
                % ----------------------
                tmpobj  = findobj('parent', hdl, 'tag', 'chan_list');
                tmptext = get(tmpobj, 'string');
                tmptext{end+1} = [ 'All ' STUDY.changrp(end).name ];
                set(tmpobj, 'string', tmptext, 'value', length(tmptext));

            case 'edit_group'
                if length(changrp) > 1, return; end
                if length(STUDY.changrp(changrp).channels) < 2, return; end
                channames = STUDY.changrp(changrp).channels;
                for i=1:length(channames), channames{i} = [ ' ' channames{i} ]; end
                channamestr = strcat(channames{:});
                res = inputdlg2({ 'Name of channel group', 'Channels to group' }, 'Create channel group', ...
                                1, { STUDY.changrp(changrp).name channamestr(2:end) });
                if isempty(res), return; end
                STUDY.changrp(end+1).name = '';
                STUDY.changrp(changrp)    = STUDY.changrp(end);
                STUDY.changrp(end)        = [];
                STUDY.changrp(changrp).name = res{1};
                allchans(changrp)         = { res{1} };
                chanlabels = parsetxt(res{2});
                STUDY.changrp(changrp).channels = chanlabels;
                tmp = std_chanlookup( STUDY, ALLEEG, STUDY.changrp(end));
                STUDY.changrp(changrp).chaninds = tmp.chaninds;
                userdat{1}{2} = STUDY;
                userdat{1}{4} = allchans;
                set(hdl, 'userdat',userdat);    

                % list of channel groups
                % ----------------------
                show_options = {};
                for index = 1:length(STUDY.changrp)
                    show_options{end+1} = [ 'All ' STUDY.changrp(index).name ];
                end
                tmpobj  = findobj('parent', hdl, 'tag', 'chan_list');
                set(tmpobj, 'string', show_options, 'value', changrp);

            case 'delete_group'
                if length(changrp) > 1, return; end
                if length(STUDY.changrp(changrp).channels) < 2, return; end
                STUDY.changrp(changrp)    = [];

                % list of channel groups
                % ----------------------
                show_options = {};
                for index = 1:length(STUDY.changrp)
                    show_options{end+1} = [ 'All ' STUDY.changrp(index).name ];
                end
                tmpobj  = findobj('parent', hdl, 'tag', 'chan_list');
                set(tmpobj, 'string', show_options, 'value', changrp-1);

            case 'renamechan'
                STUDY.saved = 'no';
                chan_name_list = get(findobj('parent', hdl, 'tag', 'chan_list'), 'String');
                chan_num = get(findobj('parent', hdl, 'tag', 'chan_list'), 'Value') -1;
                if chan_num == 0  % 'all subjects' option 
                    return;
                end
                % Don't rename 'Notchan' and 'Outliers'  channels.
                if strncmpi('Notchan',STUDY.channel(cls(chan_num)).name,8) || strncmpi('Outliers',STUDY.channel(cls(chan_num)).name,8) || ...
                        strncmpi('Parentchannel',STUDY.channel(cls(chan_num)).name,13)
                    warndlg2('The Parentchannel, Outliers, and Notchan channels cannot be renamed');
                    return;
                end
                old_name = STUDY.channel(cls(chan_num)).name;
                rename_param  = inputgui( { [1] [1] [1]}, ...
                    { {'style' 'text' 'string' ['Rename ' old_name] 'FontWeight' 'Bold'} {'style' 'edit' 'string' '' 'tag' 'chan_rename' } {} }, ...
                '', 'Rename channel - from pop_chanplot()' );
                if ~isempty(rename_param) %if not canceled
                    new_name = rename_param{1};
                    STUDY = std_renamechan(STUDY, ALLEEG, cls(chan_num), new_name);
                    % update Study history
                    a = ['STUDY = std_renamechan(STUDY, ALLEEG, ' num2str(cls(chan_num)) ', ' STUDY.channel(cls(chan_num)).name  ');'];
                    STUDY.tmphist =  sprintf('%s\n%s',  STUDY.tmphist, a);  

                    new_name = [ STUDY.channel(cls(chan_num)).name ' (' num2str(length(STUDY.channel(cls(chan_num)).onechan))  ' ICs)'];
                    chan_name_list{chan_num+1} = renamechan( chan_name_list{chan_num+1}, new_name);
                    set(findobj('parent', hdl, 'tag', 'chan_list'), 'String', chan_name_list);
                    set(findobj('parent', hdl, 'tag', 'chan_rename'), 'String', '');
                    userdat{1}{2} = STUDY;
                    set(hdl, 'userdat',userdat); %update STUDY
                end            
        end
    catch
        eeglab_error;
    end
end

function newname = renamechan(oldname, newname);
    
    tmpname = deblank(oldname(end:-1:1));
    strpos  = strfind(oldname, tmpname(end:-1:1));
    
    newname = [ oldname(1:strpos-1) newname ];
