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
%                to plot selected measures of the selected channel. Additional editing 
%                options (renaming the channel, rejecting outliers, moving channels to 
%                another channel) are also available. The option 'All N channel centroids' 
%                at the top of the list displays all the channels in the list except the 
%                'Notchannel', 'Outlier' and 'Parentchannel' channels. Selecting this option 
%                will plot the channel centroids (i.e. ERP, ERSP, ...) in a single figure.
%  "Select channel(s) to plot" - [list box] Displays the ICA channels of the currently 
%                selected channel (in the "Select channel to plot" list box). Each channel 
%                has the format: 'subject name, channel index'. Multiple channels can be 
%                selected from the list. Use the plotting buttons below to plot different 
%                measures of the selected channels on different figures. Selecting the 
%                "all subjects" option is  equivalent to using the channel plotting buttons. 
%                Additional editing options are reassigning the selected channels to 
%                another channel or moving them to the outlier channel.
%  "Plot channel properties" - [button] Displays in one figure all the mean channel measures
%                (e.g., dipole locations, scalp maps, spectra, etc.) that were calculated
%                and saved in the EEG datsets. If there is more than one condition, the ERP 
%                and the spectrum will have different colors for each condition. The ERSP 
%                and ITC plots will show only the first condition; clicking on the subplot 
%                will open a new figure with the different conditions displayed together. 
%                Uses the command line function std_propplot().
%  "Plot scalp maps"  - [button] Displays the scalp maps of channel channels.
%                If applied to a channel, scalp maps of the channel channels
%                are plotted along with the channel mean scalp map in one figure. 
%                If "All # channel centroids" option is selected, all channel scalp map
%                means are plotted in the same figure. If applied to channels, displays
%                the scalp maps of the specified channel channels in separate figures.
%                Uses the command line functions std_plotmap() and std_plotchanmap().
%  "Plot ERSPs" - [button] Displays the channel channel ERSPs. 
%                If applied to a channel, channel ERSPs are plotted in one figure  
%                (per condition) with the channel mean ERSP. If "All # channel centroids" 
%                option is selected, plots all average ERSPs of the channels in one figure 
%                per condition. If applied to channels, display the ERSP images of specified 
%                channel channels in separate figures, using one figure for all conditions.
%                Uses the command line functions std_plotersp() and std_plotchannelsp().
%  "Plot ITCs" - [button] Same as  "Plot ERSPs" but with ITC.
%                Uses the command line functions std_plotitc() and std_plotchanitc().
%  "Plot dipoles" - [button] Displays the dipoles of the channel channels.
%                If applied to a channel, plots the channel channel dipoles (in blue) 
%                plus the average channel dipole (in red). If "All # channel centroids" option 
%                is selected, all channel plots are displayed in one figure each channel in 
%                a separate subplot. If applied to channels, displays the ERSP images of the
%                specified channel. For specific channels displays channels dipole (in blue) 
%                plus the average channel dipole (in Red) in separate figures. 
%                Uses the command line functions std_dipplot() and std_plotchandip().
%  "Plot spectra" - [button] Displays the channel channel spectra.   
%                If applied to a channel, displays channel spectra plus the average channel 
%                spectrum in bold. For a specific channel, displays the channel channel 
%                spectra plus the average channel spectrum (in bold) in one figure per condition.
%                If the "All # channel centroids" option is selected, displays the average 
%                spectrum of all channels in the same figure, with spectrum for different 
%                conditions (if any) plotted in different colors.  
%                If applied to channels, displays the spectrum of specified channel 
%                channels in separate figures using one figure for all conditions.  
%                Uses the command line functions std_plotspec() and std_plotonechanpec().
%  "Plot ERPs" - [button] Same as "Plot spectra" but for ERPs.
%                Uses the command line functions std_ploterp() and std_plotchannelp().
%  "Create new channel" - [button] Creates a new empty channel.
%                Opens a popup window in which a name for the new channel can be entered.
%                If no name is given the default name is 'Cls #', where '#' is the next
%                available channel number. For changes to take place, press the popup 
%                window 'OK' button, else press the 'Cancel' button. After the empty 
%                channel is created, channels can be moved into it using, 
%                'Reassign selected channel(s)' (see below). Uses the command line 
%                function std_createchant().
%  "Rename selected channel" - [button] Renames a channel using the selected (mnemonic) name. 
%                Opens a popup window in which a new name for the selected channel can be 
%                entered. For changes to take place, press the popup window 'OK' button, 
%                else press the 'Cancel' button. Uses the command line function std_renamechant().
%  "Reject outlier channels" - [button] rejects outlier channels to an outlier channel.
%                Opens a popup window to specify the outlier threshold. Move outlier 
%                channels that are more than x standard deviations devs from the 
%                channel centroid to an outlier channel. For changes to take place, 
%                press the popup window 'OK' button, else press the 'Cancel' button. 
%                Uses the command line function std_rejectoutliers().
%  "Merge channels" - [button] Merges several channels into one channel.
%                Opens a popup window in which the channels to merge may be specified 
%                An optional name can be given to the merged channel. If no name is given, 
%                the default name is 'Cls #', where '#' is the next available channel number.   
%                For changes to take place, press the popup window 'OK' button, else press
%                the 'Cancel' button. Uses the command line function std_mergechant().
%  "Remove selected outlier channel(s)" - [button] Moves selected channel(s) to the 
%                outlier channel. The channels that will be moved are the ones selected 
%                in the "Select channel(s) to plot" list box. Opens a popup window in which 
%                a list of the selected channel(s) is presented. For changes to take place,
%                press the popup window 'OK' button, else press the 'Cancel' button. 
%                Uses the command line function std_moveoutlier().
%  "Reassign selected channel(s)" - [button] Moves selected channel(s) from one channel 
%                to another. The channels that will reassign are the ones selected in the
%                "Select channel(s) to plot" list box. Opens a popup window in which 
%                a list of possible channels to which to move the selected channel(s) is 
%                presented. For changes to take place, press the popup window 'OK' button, 
%                else press the 'Cancel' button. Uses the command line function std_movecomp().
%  "Save STUDY set to disk" - [check box] Saves the STUDY set structure modified according 
%                to specified user edits to the disk. If no file name is entered will
%                overwrite the current STUDY set file. 
%
% See also:  pop_prechant(), pop_chant().         
%
% Authors: Arnaud Delorme, Hilit Serby, Scott Makeig, SCCN/INC/UCSD, October 11, 2004

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, October 11, 2004, hilit@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% Coding notes: Useful information on functions and global variables used.

function [STUDY, com] = pop_chanplot(varargin)

icadefs;
com = [];
if ~isstr(varargin{1})
    if nargin < 2
        error('pop_chanplot(): You must provide ALLEEG and STUDY structures');
    end
    STUDY  = varargin{1};
    STUDY.etc.erpparams.topotime  = []; % [] for channels and NaN for components
    STUDY.etc.specparams.topofreq = []; % NaN -> GUI disabled
    STUDY.etc.erspparams.topotime = [];
    STUDY.etc.erspparams.topofreq = [];
    
    % test path
    % ---------
    pathwarn = 'off';
    if ~strcmpi(pwd, STUDY.filepath)
        if length(STUDY.datasetinfo(1).filepath) < 1
            pathwarn = 'on';
        elseif STUDY.datasetinfo(1).filepath(1) == '.'
            pathwarn = 'on';
        end;
    end;
    if strcmpi(pathwarn, 'on')
        warndlg2(strvcat('You have changed your working path and data files are', ...
                         'no longer available; Cancel, and go back to your STUDY folder'), 'warning');
    end;
        
    oldhistory = STUDY.history;
    STUDY.history = '';
    ALLEEG = varargin{2};
    if ~isfield(STUDY, 'changrp')
        STUDY = std_changroup(STUDY, ALLEEG);
        disp('Warning: history not saved for group creation');
    elseif isempty(STUDY.changrp)
        STUDY = std_changroup(STUDY, ALLEEG);
        disp('Warning: history not saved for group creation');
    end;
    
    show_chan          = [ 'pop_chanplot(''showchan'',gcf);'];
    show_onechan       = [ 'pop_chanplot(''showchanlist'',gcf);'];
	plot_chan_maps     = [ 'pop_chanplot(''topoplot'',gcf); ']; 
    plot_onechan_maps  = [ 'pop_chanplot(''plotchantopo'',gcf); ']; 
    plot_chan_ersps    = ['pop_chanplot(''erspplot'',gcf); '];
    plot_onechan_ersps = ['pop_chanplot(''plotchanersp'',gcf); '];
    plot_chan_itcs     = ['pop_chanplot(''itcplot'',gcf); '];
    plot_onechan_itcs  = ['pop_chanplot(''plotchanitc'',gcf); '];
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
    ersp_opt           = ['pop_chanplot(''ersp_opt'',gcf);'];
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
    end;
                           
    % Create default ERSP / ITC time/freq. paramters 
    % ----------------------------------------------
    % Find the first entry in STUDY.setind that is not NaN
    ref_ind = 0;
    found_ind1 = 0;
    for ri = 1:size(STUDY.setind,1)
        for ci = 1:size(STUDY.setind,2)
            if ~isnan(STUDY.setind(ri,ci))
                ref_ind = STUDY.setind(ri,ci);
                found_ind1 = 1;
                break
            end
        end
        if found_ind1 == 1, break; end
    end
    if ref_ind == 0     % If STUDY.setind contains only NaNs, or is empty.
        error('STUDY contains no datasets');
    end

    % enable buttons
    % --------------
    filename = fullfile( ALLEEG(ref_ind).filepath, ALLEEG(ref_ind).filename(1:end-3));
    if exist([filename 'datspec']) , spec_enable = 'on'; else spec_enable  = 'off'; end;
    if exist([filename 'daterp'] )  , erp_enable = 'on'; else erp_enable   = 'off'; end;
    if exist([filename 'datersp']) , ersp_enable = 'on'; else ersp_enable  = 'off'; end;
    
    if isfield(ALLEEG(1).dipfit, 'model'), dip_enable   = 'on'; else dip_enable   = 'off'; end;
    
    % userdata below
    % --------------
    fig_arg{1}{1} = ALLEEG;
    fig_arg{1}{2} = STUDY;
    fig_arg{1}{3} = STUDY.changrp;
    fig_arg{1}{4} = { STUDY.changrp.name };
    fig_arg{2}    = length(STUDY.changrp);
        
    geometry = { [4] [1] [0.7 0.3 0.3 0.1 0.9] [1 0.3 1] [1 0.3 1] [1 0.3 1] [1 0.3 1] ...
                 [1 0.3 1] [1 0.3 1] [1] [1 0.3 1] [1 0.3 1] [1] [0.2 1 1.5 0.3] [1] };
    uilist   = { ...
        {'style' 'text' 'string' ['Study ''' STUDY.name '''' ] ...
            'FontWeight' 'Bold' 'HorizontalAlignment' 'center'} {} ...
        {'style' 'text'       'string' 'Select channel to plot' 'FontWeight' 'Bold' } ...
        {'style' 'pushbutton' 'string' 'Sel. all' 'callback' sel_all_chans } {} {} ...
        {'style' 'text'       'string' 'Select subject(s) to plot' 'FontWeight' 'Bold'} ...
        {'style' 'listbox'    'string' show_options 'value' 1 'max' 2 'tag' 'chan_list' 'Callback' show_chan } {} ...
        {'style' 'listbox'    'string' '' 'tag' 'chan_onechan' 'max' 2 'min' 1 'callback'    show_onechan } ... 
        {'style' 'pushbutton' 'enable'   erp_enable 'string' 'Plot ERPs' 'Callback' plot_chan_erp} ...
        {'style' 'pushbutton' 'enable'   erp_enable 'string' 'Params' 'Callback' erp_opt }  ...
        {'style' 'pushbutton' 'enable'   erp_enable 'string' 'Plot ERP(s)' 'Callback' plot_onechan_erp} ...
        {'style' 'pushbutton' 'enable'  spec_enable 'string' 'Plot spectra' 'Callback' plot_chan_spectra} ...
        {'style' 'pushbutton' 'enable'  spec_enable 'string' 'Params' 'Callback' spec_opt }  ...
        {'style' 'pushbutton' 'enable'  spec_enable 'string' 'Plot spectra' 'Callback' plot_onechan_spectra} ...
        {'style' 'pushbutton' 'enable'  ersp_enable 'string' 'Plot ERSPs' 'Callback' plot_chan_ersps} ...
        {'vertshift' 'style' 'pushbutton' 'enable'  ersp_enable 'string' 'Params' 'Callback' ersp_opt }  ...
        {'style' 'pushbutton' 'enable'  ersp_enable 'string' 'Plot ERSP(s)' 'Callback' plot_onechan_ersps}...
        {'style' 'pushbutton' 'enable'  ersp_enable 'string' 'Plot ITCs' 'Callback' plot_chan_itcs} { }  ...
        {'style' 'pushbutton' 'enable'  ersp_enable 'string' 'Plot ITC(s)' 'Callback' plot_onechan_itcs}...
        {'style' 'pushbutton' 'string' 'Plot channel properties' 'Callback' plot_chan_sum} {} ... 
        {'style' 'pushbutton' 'string' 'Plot channel properties (soon)' 'Callback' plot_onechan_sum 'enable' 'off'} {} ...
        {'style' 'pushbutton' 'string' 'Create channel group (soon)' 'Callback' create_group 'enable' 'off'} {} ...
        {'style' 'pushbutton' 'string' 'Delete channel group (soon)' 'Callback' delete_group 'enable' 'off'} ...
        {'style' 'pushbutton' 'string' 'Edit channel group (soon)' 'Callback' edit_group 'enable' 'off'} {} {} {} ...
        {'style' 'checkbox'   'string' '' 'tag' 'saveSTUDY' 'Callback' saveSTUDY 'value' 0} ...
        {'style' 'text'       'string' 'Save STUDY set to disk'} ...
        {'style' 'edit'       'string' fullfile(STUDY.filepath, STUDY.filename) 'enable' 'off' 'tag' 'studyfile' 'userdata' 'save'} ...
        {'style' 'pushbutton' 'string' '...' 'tag' 'browsesave' 'Callback' browsesave 'enable' 'off' 'userdata' 'save'} {} };
    
   [out_param userdat] = inputgui( 'geometry' , geometry, 'uilist', uilist, ...
                                   'helpcom', 'pophelp(''pop_chanoutput'')', ...
                                   'title', 'View and edit current channels -- pop_chanplot()' , 'userdata', fig_arg, ...
                                   'geomvert', [ 1 1 1 5 1 1 1 1 1 1 1 1 1 1 1], 'eval', show_chan );
	
   if ~isempty(userdat)
       ALLEEG = userdat{1}{1};
       STUDY = userdat{1}{2};
       % If save updated STUDY to disk
        if out_param{3}
            if ~isempty(out_param{4})
                [filepath filename ext] = fileparts(out_param{4});
                a = sprintf('%s%s%s%s%s%s', 'STUDY = pop_savestudy(STUDY, ALLEEG' , '''filename'', ''', [filename ext], ''', ''filepath'', ''', filepath, ''');' );
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);              
                STUDY = pop_savestudy(STUDY, ALLEEG, 'filename', [filename ext], 'filepath', filepath);
            else
                if (~isempty(STUDY.filename)) & (~isempty(STUDY.filepath))
                    a = sprintf('%s%s%s%s%s%s', 'STUDY = pop_savestudy(STUDY, ALLEEG, ' , '''filename'', ''', STUDY.filename, ''', ''filepath'', ''', STUDY.filepath, ''');' );
                    STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);              
                    STUDY = pop_savestudy(STUDY, ALLEEG, 'filename', STUDY.filename, 'filepath', STUDY.filepath);
                else
                    a = sprintf('%s%s%s%s%s%s', 'STUDY = pop_savestudy(STUDY, ALLEEG);' );
                    STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);              
                    STUDY = pop_savestudy(STUDY, ALLEEG);
                end
           end
        end
   end

   % history
   % -------
   com = STUDY.history;
   STUDY.history =  sprintf('%s%s', oldhistory, com);              
   
else
    hdl = varargin{2};  %figure handle
    userdat  = get(varargin{2}, 'userdat');    
    ALLEEG   = userdat{1}{1};
    STUDY    = userdat{1}{2};
    cls      = userdat{1}{3};
    allchans = userdat{1}{4};
    
    changrp = get(findobj('parent', hdl, 'tag', 'chan_list')   , 'value');
    onechan = get(findobj('parent', hdl, 'tag', 'chan_onechan'), 'value');
   
    switch  varargin{1}
        
        case {'topoplot', 'erspplot','itcplot','specplot', 'erpplot'}
            changrpstr = allchans(changrp);
            plotting_option = varargin{1};
            plotting_option = [ plotting_option(1:end-4) 'plot' ];
            a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''channels'','  vararg2str({changrpstr}) ');' ];
             % update Study history
            eval(a); STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);  
            userdat{1}{2} = STUDY;
            set(hdl, 'userdat',userdat); 

        case {'plotchantopo', 'plotchanersp','plotchanitc','plotchanspec', 'plotchanerp','plotchandip'}
            changrpstr    = allchans(changrp);
            %if length(changrp) > 1
            %    subject = STUDY.subject{onechan-1};
            %else
            %    changrpstruct = STUDY.changrp(changrp);
            %    allsubjects   = unique({ STUDY.datasetinfo([ changrpstruct.setinds{:} ]).subject });
            %    subject = allsubjects{onechan-1};
            %end;
            plotting_option = varargin{1};
            plotting_option = [ plotting_option(9:end) 'plot' ];
            if onechan(1) ~= 1  % check that not all onechan in channel are requested
                 subject = STUDY.subject{onechan-1};
                 a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''channels'','  vararg2str({changrpstr}) ', ''subject'', ''' subject ''' );' ];
                 eval(a); STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);  
             else
                a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''channels'','  vararg2str({changrpstr}) ', ''plotsubjects'', ''on'' );' ];
                eval(a); STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);
             end;
            userdat{1}{2} = STUDY;
            set(hdl, 'userdat',userdat); 
            
        case 'erp_opt' % save the list of selected channels
            [STUDY com] = pop_erpparams(STUDY);
            if ~isempty(com)
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, com);
            end;
            userdat{1}{2} = STUDY;
            set(hdl, 'userdat',userdat); %update information (STUDY)     

        case 'spec_opt' % save the list of selected channels
            [STUDY com] = pop_specparams(STUDY);
            if ~isempty(com)
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, com);
            end;
            userdat{1}{2} = STUDY;
            set(hdl, 'userdat',userdat); %update information (STUDY)     
         
        case 'ersp_opt' % save the list of selected channels
            [STUDY com] = pop_erspparams(STUDY);
            if ~isempty(com)
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, com);
            end;
            userdat{1}{2} = STUDY;
            set(hdl, 'userdat',userdat); %update information (STUDY)     
            
        case 'showchanlist' % save the list of selected channels
            if length(changrp) == 1
                STUDY.changrp(changrp).selected = onechan;
            end;
            userdat{1}{2} = STUDY;
            set(hdl, 'userdat',userdat); %update information (STUDY)     
               
       case 'showchan'
            cind     = get(findobj('parent', hdl, 'tag', 'chan_list')   , 'value');
            changrp  = STUDY.changrp(cind);
            
            % Find datasets availaible
            % ------------------------
            %setind = STUDY.setind .* (changrp.chaninds > 0); % set to 0 the cell not
            %%                                       % containing any electrode
            %allchansets = unique( setind(find(setind(:))) );
            
            % Generate channel list
            % ---------------------
            chanid{1} = 'All subjects';
            if length(changrp) == 1
                allsubjects = unique({ STUDY.datasetinfo([ changrp.setinds{:} ]).subject });
                for l = 1:length(allsubjects)
                    chanid{end+1} = [ allsubjects{l} ' ' changrp.name ];
                end;
            else
                for l = 1:length(STUDY.subject)
                    chanid{end+1} = [ STUDY.subject{l} ];
                end;
            end;
                
            selected = 1;
            if isfield(changrp, 'selected') & length(cind) == 1
                if ~isempty(STUDY.changrp(cind).selected)
                    selected = min(STUDY.changrp(cind).selected, 1+length(chanid));
                    STUDY.changrp(cind).selected = selected;
                end;
            end;

            set(findobj('parent', hdl, 'tag', 'chan_onechan'), 'value', selected, 'String', chanid);
        
        case 'sel_all_chans'
            set(findobj('parent', hdl, 'tag', 'chan_list'), 'value', [1:length(STUDY.changrp)]);
            
            % Generate channel list
            % ---------------------
            chanid{1} = 'All subjects';
            for l = 1:length(STUDY.subject)
                chanid{end+1} = [ STUDY.subject{l} ' All' ];
            end;
            selected = 1;
            set(findobj('parent', hdl, 'tag', 'chan_onechan'), 'value', selected, 'String', chanid);

        case 'plotsum'
            changrpstr = allchans(changrp);
            [STUDY] = std_propplot(STUDY, ALLEEG, allchans(changrp));
            a = ['STUDY = std_propplot(STUDY, ALLEEG, ' vararg2str({ allchans(changrp) }) ' );'  ];
            STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);  
            userdat{1}{2} = STUDY;
            set(hdl, 'userdat',userdat);    
                       
        case 'create_group'
            channames = { STUDY.changrp(changrp).name };
            for i=1:length(channames), channames{i} = [ ' ' channames{i} ]; end;
            channamestr = strcat(channames{:});
            res = inputdlg2({ 'Name of channel group', 'Channels to group' }, 'Create channel group', 1, { '' channamestr(2:end) });
            if isempty(res), return; end;
            STUDY.changrp(end+1).name = res{1};
            allchans(end+1)         = { res{1} };
            chanlabels = parsetxt(res{2});
            if length(chanlabels) == 1
                warndlg2('Cannot create a channel group with a single channel');
                return;
            end;
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
            if length(changrp) > 1, return; end;
            if length(STUDY.changrp(changrp).channels) < 2, return; end;
            channames = STUDY.changrp(changrp).channels;
            for i=1:length(channames), channames{i} = [ ' ' channames{i} ]; end;
            channamestr = strcat(channames{:});
            res = inputdlg2({ 'Name of channel group', 'Channels to group' }, 'Create channel group', ...
                            1, { STUDY.changrp(changrp).name channamestr(2:end) });
            if isempty(res), return; end;
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
            end;
            tmpobj  = findobj('parent', hdl, 'tag', 'chan_list');
            set(tmpobj, 'string', show_options, 'value', changrp);
            
        case 'delete_group'
            if length(changrp) > 1, return; end;
            if length(STUDY.changrp(changrp).channels) < 2, return; end;
            STUDY.changrp(changrp)    = [];
            
            % list of channel groups
            % ----------------------
            show_options = {};
            for index = 1:length(STUDY.changrp)
                show_options{end+1} = [ 'All ' STUDY.changrp(index).name ];
            end;
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
            if strncmpi('Notchan',STUDY.channel(cls(chan_num)).name,8) | strncmpi('Outliers',STUDY.channel(cls(chan_num)).name,8) | ...
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
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);  
                
                new_name = [ STUDY.channel(cls(chan_num)).name ' (' num2str(length(STUDY.channel(cls(chan_num)).onechan))  ' ICs)'];
                chan_name_list{chan_num+1} = renamechan( chan_name_list{chan_num+1}, new_name);
                set(findobj('parent', hdl, 'tag', 'chan_list'), 'String', chan_name_list);
                set(findobj('parent', hdl, 'tag', 'chan_rename'), 'String', '');
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat); %update STUDY
            end            
    end
end

function newname = renamechan(oldname, newname);
    
    tmpname = deblank(oldname(end:-1:1));
    strpos  = strfind(oldname, tmpname(end:-1:1));
    
    newname = [ oldname(1:strpos-1) newname ];
