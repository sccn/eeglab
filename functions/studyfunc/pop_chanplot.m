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
%                "All channels" option is  equivalent to using the channel plotting buttons. 
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
    oldhistory = STUDY.history;
    STUDY.history = '';
    ALLEEG = varargin{2};
    if ~isfield(STUDY, 'changrp')
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
    merge_channels      = ['pop_chanplot(''mergechannels'',gcf);'];
    erp_opt            = ['pop_chanplot(''erp_opt'',gcf);'];
    spec_opt           = ['pop_chanplot(''spec_opt'',gcf);'];
    ersp_opt           = ['pop_chanplot(''ersp_opt'',gcf);'];
    saveSTUDY          = [ 'set(findobj(''parent'', gcbf, ''userdata'', ''save''), ''enable'', fastif(get(gcbo, ''value'')==1, ''on'', ''off''));' ];
    browsesave         = [ '[filename, filepath] = uiputfile2(''*.study'', ''Save STUDY with .study extension -- pop_chan()''); ' ... 
                           'set(faindobj(''parent'', gcbf, ''tag'', ''studyfile''), ''string'', [filepath filename]);' ];
    
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
    if exist([filename 'icatopo']), scalp_enable = 'on'; else scalp_enable = 'off'; end;
    if exist([filename 'icaersp']) , ersp_enable = 'on'; else ersp_enable  = 'off'; end;
    
    if isfield(ALLEEG(1).dipfit, 'model'), dip_enable   = 'on'; else dip_enable   = 'off'; end;
    
    % userdata below
    % --------------
    fig_arg{1}{1} = ALLEEG;
    fig_arg{1}{2} = STUDY;
    fig_arg{1}{3} = STUDY.changrp;
    fig_arg{1}{4} = { STUDY.changrp.name };
    fig_arg{2}    = length(STUDY.changrp);
        
    geometry = { [4] [1] [1 0.3 1] [1 0.3 1] [1 0.3 1] [1 0.3 1] [1 0.3 1] ...
                 [1 0.3 1] [1 0.3 1] [1] [0.2 1 1.5 0.3] [1] };
    uilist   = { ...
        {'style' 'text' 'string' ['Study ''' STUDY.name '''' ] ...
            'FontWeight' 'Bold' 'HorizontalAlignment' 'center'} {} ...
        {'style' 'text'       'string' 'Select channel to plot' 'FontWeight' 'Bold' } {} ...
        {'style' 'text'       'string' 'Select channel(s) to plot' 'FontWeight' 'Bold'} ...
        {'style' 'listbox'    'string' show_options 'value' 1 'tag' 'chan_list' 'Callback' show_chan } {} ...
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
        {'style' 'pushbutton' 'string' 'Plot channel properties' 'Callback' plot_onechan_sum 'enable' 'off'} {} ...
        {'style' 'checkbox'   'string' '' 'tag' 'saveSTUDY' 'Callback' saveSTUDY 'value' 0} ...
        {'style' 'text'       'string' 'Save STUDY set to disk'} ...
        {'style' 'edit'       'string' fullfile(STUDY.filepath, STUDY.filename) 'enable' 'off' 'tag' 'studyfile' 'userdata' 'save'} ...
        {'style' 'pushbutton' 'string' '...' 'tag' 'browsesave' 'Callback' browsesave 'enable' 'off' 'userdata' 'save'} {} };
    
   [out_param userdat] = inputgui( 'geometry' , geometry, 'uilist', uilist, ...
                                   'helpcom', 'pophelp(''pop_chanoutput'')', ...
                                   'title', 'View and edit current channels -- pop_chanplot()' , 'userdata', fig_arg, ...
                                   'geomvert', [ 1 1 1 3 1 1 1 1 1 1 1 1], 'eval', show_chan );
	
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
            changrpstr = allchans(changrp);
            plotting_option = varargin{1};
            plotting_option = [ plotting_option(9:end) 'plot' ];
            if onechan(1) ~= 1  % check that not all onechan in channel are requested
                a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''channels'','  vararg2str({changrpstr}) ', ''subject'', ''' STUDY.subject{onechan-1} ''' );' ];
                eval(a); STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);  
             else
                a = ['STUDY = std_' plotting_option '(STUDY,ALLEEG,''channels'','  vararg2str({changrpstr}) ', ''plotsubjects'', ''on'' );' ];
                eval(a); STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);
             end;
            userdat{1}{2} = STUDY;
            set(hdl, 'userdat',userdat); 
            
        case 'erp_opt' % save the list of selected channels
            STUDY = pop_erpparams(STUDY);
            userdat{1}{2} = STUDY;
            set(hdl, 'userdat',userdat); %update information (STUDY)     

        case 'spec_opt' % save the list of selected channels
            STUDY = pop_specparams(STUDY);
            userdat{1}{2} = STUDY;
            set(hdl, 'userdat',userdat); %update information (STUDY)     
         
        case 'ersp_opt' % save the list of selected channels
            STUDY = pop_erspparams(STUDY);
            userdat{1}{2} = STUDY;
            set(hdl, 'userdat',userdat); %update information (STUDY)     
            
        case 'showchanlist' % save the list of selected channels
            STUDY.changrp(changrp).selected = onechan;
            userdat{1}{2} = STUDY;
            set(hdl, 'userdat',userdat); %update information (STUDY)     
               
       case 'showchan'
            cind     = get(findobj('parent', hdl, 'tag', 'chan_list')   , 'value');
            changrp  = STUDY.changrp(cind);
            
            len = length(STUDY.changrp(cind).chaninds);
            chanid{1} = 'All channels';

            % Find datasets availaible
            % ------------------------
            %setind = STUDY.setind .* (changrp.chaninds > 0); % set to 0 the cell not
            %%                                       % containing any electrode
            %allchansets = unique( setind(find(setind(:))) );
            
            % Generate channel list
            % ---------------------
            for l = 1:length(STUDY.subject)
                chanid{end+1} = [ STUDY.subject{l} ' ' changrp.name ];
            end;
            selected = 1;
            if isfield(changrp, 'selected')
                if ~isempty(STUDY.changrp(cind).selected)
                    selected = min(STUDY.changrp(cind).selected, 1+length(chanid));
                    STUDY.changrp(cind).selected = selected;
                end;
            end;

            set(findobj('parent', hdl, 'tag', 'chan_onechan'), 'value', selected, 'String', chanid);
           
        case 'plotsum'
            changrpstr = allchans(changrp);
            [STUDY] = std_propplot(STUDY, ALLEEG, allchans(changrp));
            a = ['STUDY = std_propplot(STUDY, ALLEEG, ' vararg2str(allchans(changrp)) ' );'  ];
            STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);  
            userdat{1}{2} = STUDY;
            set(hdl, 'userdat',userdat);    
           
        case 'plotonechanum'
            for ci = 1 : length(comp_ind)
            end
            
        case 'renamechan'
            STUDY.saved = 'no';
            chan_name_list = get(findobj('parent', hdl, 'tag', 'chan_list'), 'String');
            chan_num = get(findobj('parent', hdl, 'tag', 'chan_list'), 'Value') -1;
            if chan_num == 0  % 'all channels' option 
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
            
        case 'movecomp'
            STUDY.saved = 'no';
            old_chan = get(findobj('parent', hdl, 'tag', 'chan_list'), 'value') -1;
            comp_ind = get(findobj('parent', hdl, 'tag', 'chan_onechan'), 'Value'); 
            if old_chan == 0 % 'all channels' option 
                return;
            end
            % Don't reassign channels of 'Notchan' or the 'Parentchannel'.
            if strncmpi('Parentchannel',STUDY.channel(cls(old_chan)).name,13)  
                warndlg2('Cannot reassign channels of ''Parentchannel''.');
                return;
			end
            old_name = STUDY.channel(cls(old_chan)).name;
            ncomp = length(comp_ind); % number of selected channels
            optionalcls =[];
            for k = 1:length(cls)
                if (~strncmpi('Parentchannel',STUDY.channel(cls(k)).name,13))  & (k~= old_chan)
                    optionalcls = [optionalcls cls(k)];
                end
            end                    
            reassign_param  = inputgui( { [1] [1] [1]}, ...
                { {'style' 'text' 'string' strvcat(['Reassign ' fastif(ncomp >1, [num2str(length(comp_ind)) ' currently selected channels'], ...
                                                              'currently selected channel') ], ...
                            [' from ' old_name ' to the channel selected below']) 'FontWeight' 'Bold'} ...
                  {'style' 'listbox' 'string' {STUDY.channel(optionalcls).name} 'tag' 'new_chan'} {} }, ...
                  '', 'Reassign channel - from pop_chanplot()' ,[] , 'normal', [2 3 1] );
            if ~isempty(reassign_param) %if not canceled
                new_chan = reassign_param{1};
                comp_to_disp = get(findobj('parent', hdl, 'tag', 'chan_onechan'), 'String');      
                if strcmp(comp_to_disp{comp_ind(1)},'All channels')
                    warndlg2('Cannot move all the channels of the channel - abort move channels', 'Aborting move channels');
                    return;
                end
                STUDY = std_movecomp(STUDY, ALLEEG,  cls(old_chan), optionalcls(new_chan), comp_ind - 1);                
                % update Study history
                a = ['STUDY = std_movecomp(STUDY, ALLEEG, ' num2str(cls(old_chan)) ', ' num2str(optionalcls(new_chan)) ', [' num2str(comp_ind - 1) ']);'];
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);
                newind = find(cls == optionalcls(new_chan));
                
                % update GUI
                % ----------
                chan_name_list = get(findobj('parent', hdl, 'tag', 'chan_list'), 'String');     
                newname = [STUDY.channel(optionalcls(new_chan)).name ' (' num2str(length(STUDY.channel(optionalcls(new_chan)).onechan))  ' ICs)'];
                chan_name_list{newind+1} = renamechan( chan_name_list{newind+1}, newname);
                newname = [STUDY.channel(cls(old_chan)).name ' (' num2str(length(STUDY.channel(cls(old_chan)).onechan))  ' ICs)'];
                chan_name_list{old_chan+1} = renamechan( chan_name_list{old_chan+1}, newname);
                set( findobj('parent', hdl, 'tag', 'chan_list'), 'String', chan_name_list);
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat); 
                pop_chanplot('showchan',hdl);
            end          
            
        case 'moveoutlier'
            STUDY.saved = 'no';
            old_chan = get(findobj('parent', hdl, 'tag', 'chan_list'), 'value') -1;
            comp_ind = get(findobj('parent', hdl, 'tag', 'chan_onechan'), 'Value'); 
            if ~isempty(find(comp_ind ==1))
                warndlg2('Cannot remove all the channel channels');
                return;
            end
            if old_chan == 0 % 'all channels' option 
                return;
            end
            if strncmpi('Notchan',STUDY.channel(cls(old_chan)).name,8) | strncmpi('Parentchannel',STUDY.channel(cls(old_chan)).name,13)    % There are no outliers to 'Notchan'
                warndlg2('Cannot reassign channels of ''Notchan'' or ''Parentchannel''.');
                return;
            end
            comp_list = get(findobj('parent', hdl, 'tag', 'chan_onechan'), 'String'); 
            ncomp = length(comp_ind);
            old_name = STUDY.channel(cls(old_chan)).name;            
            if strncmpi('Outliers',STUDY.channel(cls(old_chan)).name,8)  % There are no outliers of 'Outliers'
                warndlg2('Cannot use ''Outliers'' channels for this option.');
                return;
			end
            reassign_param  = inputgui( { [1] [1] [1]}, ...
                { {'style' 'text' 'string' ['Remove ' fastif(ncomp >1, [num2str(length(comp_ind)) ' currently selected channels below '], 'currently selected channel below ') ...
                            'from ' old_name ' to its outlier channel?'] 'FontWeight' 'Bold'} ...
                  {'style' 'listbox' 'string' {comp_list{comp_ind}} 'tag' 'new_chan'} {} }, ...
                  '', 'Remove outliers - from pop_chanplot()' ,[] , 'normal', [1 3 1] );
            if ~isempty(reassign_param) %if not canceled
                STUDY = std_moveoutlier(STUDY, ALLEEG,  cls(old_chan), comp_ind - 1);
                chan_name_list = get(findobj('parent', hdl, 'tag', 'chan_list'), 'String');
                outlier_chan = std_findoutlierchan(STUDY,cls(old_chan)); %find the outlier channel for this channel
                oind = find(cls == outlier_chan); % the outlier chan index (if already exist) in the channel list GUI
                if ~isempty(oind) % the outlier chan is already presented in the channel list GUI
                    newname = [STUDY.channel(outlier_chan).name ' (' num2str(length(STUDY.channel(outlier_chan).onechan))  ' ICs)'];
                    chan_name_list{oind+1} = renamechan( chan_name_list{oind+1}, newname);
                elseif outlier_chan == length(STUDY.channel) % update the list with the Outlier channel (if didn't exist before)
                    chan_name_list{end+1} = [STUDY.channel(outlier_chan).name ' (' num2str(length(STUDY.channel(outlier_chan).onechan))  ' ICs)'];
                    userdat{2} = userdat{2} + 1; % update N, number of channels in edit window 
                    cls(end +1) = length(STUDY.channel); % update the GUI channels list with the outlier channel
                    userdat{1}{3} = cls;  % update cls, the channel indices in edit window
                end
                newname = [STUDY.channel(cls(old_chan)).name ' (' num2str(length(STUDY.channel(cls(old_chan)).onechan))  ' ICs)'];
                chan_name_list{old_chan+1} = renamechan(chan_name_list{old_chan+1}, newname);
                set(findobj('parent', hdl, 'tag', 'chan_list'), 'String', chan_name_list);
                % update Study history
                a = ['STUDY = std_moveoutlier(STUDY, ALLEEG, ' num2str(cls(old_chan)) ', [' num2str(comp_ind - 1) ']);'];
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat); 
                pop_chanplot('showchan',hdl);    
            end
            
        case 'rejectoutliers'
            STUDY.saved = 'no';
            chan = get(findobj('parent', hdl, 'tag', 'chan_list'), 'Value') -1;
            if chan
                std_name = STUDY.channel(cls(chan)).name;
                % Cannot reject outliers from 'Notchan', 'Parentchannel' and 'Outlier' channels
                if strncmpi('Notchan',std_name,8) | strncmpi('Parentchannel', std_name,13) | ...
                        strncmpi('Outliers',std_name,8)
                    warndlg2('Cannot reject outliers of ''Notchan'' or ''Outliers'' or ''Parentchannel'' channels.');
                    return;
			    end
                channels = cls(chan);
            else
                std_name = 'All channels';
                channels = [];
                for k = 1:length(cls)
                     if ~strncmpi('Notchan',STUDY.channel(cls(k)).name,8) & ~strncmpi('Outliers',STUDY.channel(cls(k)).name,8) & ...
                             ~strncmpi('Parentchannel',STUDY.channel(cls(k)).name,13)  
                        channels = [ channels cls(k)];
                    end
                end
            end
            reject_param  = inputgui( { [1] [1] [4 1 2] [1]}, ...
                { {'style' 'text' 'string' ['Reject "' std_name  '" outliers ' ] 'FontWeight' 'Bold'} {} ...
                   {'style' 'text' 'string' 'Move outlier channels that are more than'} {'style' 'edit' 'string' '3' 'tag' 'outliers_std' } ...
                  {'style' 'text' 'string' 'standard deviations' } ...
                  {'style' 'text' 'string' [ 'from the "' std_name  '" centroid to an outlier channel.']} }, ...
                  '', 'Reject outliers - from pop_chanplot()' );
            if ~isempty(reject_param) %if not canceled
                ostd = reject_param{1}; % the requested outlier std
                [STUDY] = std_rejectoutliers(STUDY, ALLEEG, channels, str2num(ostd));  
                % update Study history
                a = ['STUDY = std_rejectoutliers(STUDY, ALLEEG, [ ' num2str(channels) ' ], ' ostd ');'];
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);
                chan_name_list = get(findobj('parent', hdl, 'tag', 'chan_list'), 'String');
                for k = 1:length(channels)
                    outlier_chan = std_findoutlierchan(STUDY,channels(k)); %find the outlier channel for this channel
                    oind = find(cls == outlier_chan); % the outlier chan index (if already exist) in the channel list GUI
                    if ~isempty(oind) % the outlier chan is already presented in the channel list GUI
                        newname = [STUDY.channel(outlier_chan).name ' (' num2str(length(STUDY.channel(outlier_chan).onechan))  ' ICs)'];
                        chan_name_list{oind+1} = renamechan( chan_name_list{oind+1}, newname);
                    else % update the list with the outlier channel 
                        chan_name_list{end+1} = [STUDY.channel(outlier_chan).name ' (' num2str(length(STUDY.channel(outlier_chan).onechan))  ' ICs)'];
                        userdat{2} = userdat{2} + 1; % update N, number of channels in edit window 
                        cls(end +1) = outlier_chan; % update the GUI channels list with the outlier channel
                        userdat{1}{3} = cls;  % update cls, the channel indices in edit window
                    end
                    clsind = find(cls == channels(k));
                    newname = [STUDY.channel(channels(k)).name ' (' num2str(length(STUDY.channel(channels(k)).onechan))  ' ICs)'];
                    chan_name_list{clsind+1} = renamechan( chan_name_list{clsind+1}, newname);
                    set(findobj('parent', hdl, 'tag', 'chan_list'), 'String', chan_name_list);
                end
                % If outlier channel doesn't exist in the GUI window add it 
                userdat{1}{2} = STUDY;
                set(hdl, 'userdat',userdat); 
                pop_chanplot('showchan',hdl);
            end
            
        case 'createchan'
            STUDY.saved = 'no';
            create_param  = inputgui( { [1] [1 1] [1]}, ...
                { {'style' 'text' 'string' 'Create new empty channel' 'FontWeight' 'Bold'} ...
                   {'style' 'text' 'string' 'Enter channel name:'} {'style' 'edit' 'string' '' } {} }, ...
                  '', 'Create new empty channel - from pop_chanplot()' );
            if ~isempty(create_param) %if not canceled
                chan_name = create_param{1}; % the name of the new channel
                [STUDY] = std_createchan(STUDY, ALLEEG, chan_name); 
                % Update channel list
                chan_name_list = get(findobj('parent', hdl, 'tag', 'chan_list'), 'String');
                chan_name_list{end+1} = [STUDY.channel(end).name ' (0 ICs)']; %update the channel list with the new channel
                % update the first option on the GUI list : 'All 10 channel centroids'
                % with the new number of channel centroids
                ti = strfind(chan_name_list{1},'channel'); %get the number of channels centroid 
                cent = num2str(str2num(chan_name_list{1}(5:ti-2))+1); % new number of centroids
                chan_name_list{1} = [chan_name_list{1}(1:4) cent chan_name_list{1}(ti-1:end)]; %update list
                set(findobj('parent', hdl, 'tag', 'chan_list'), 'String', chan_name_list);
                % update Study history
                if isempty(chan_name)
                    a = ['STUDY = std_createchan(STUDY, ALLEEG);'];
                else
                    a = ['STUDY = std_createchan(STUDY, ALLEEG, ' chan_name ');'];
                end                
                STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);  
                userdat{1}{2} = STUDY;
                userdat{2} = userdat{2} + 1; % update N, the number of channel options in edit window 
                cls(end +1) = length(STUDY.channel); % update the GUI channels list with the new channel
                userdat{1}{3} = cls;  % update cls, the channel indices in edit window
                set(hdl, 'userdat',userdat); %update STUDY, cls and N
            end
            
        case 'mergechannels'
            STUDY.saved = 'no';
            chan_names = get(findobj('parent', hdl, 'tag', 'chan_list'), 'string') ;
            optionalcls =[];
            for k = 2:length(chan_names)
                if (~strncmpi('Notchan',chan_names{k},8)) & (~strncmpi('Outliers',chan_names{k},8)) & ...
                        (~strncmpi('Parentchannel',chan_names{k},13))
                    optionalcls = [optionalcls k];
                end
            end           
            reassign_param  = inputgui( { [1] [1] [1] [2 1] [1]}, ...
                { {'style' 'text' 'string' 'Select channels to Merge' 'FontWeight' 'Bold'} ...
                  {'style' 'listbox' 'string' chan_names(optionalcls) 'tag' 'new_chan' 'max' 3 'min' 1} {} ...
                  {'style' 'text' 'string' 'Optional, enter a name for the merged channel:' 'FontWeight' 'Bold'} ...
                  {'style' 'edit' 'string' ''} {} }, ...
                  '', 'Merge channels - from pop_chanplot()' ,[] , 'normal', [1 3 1 1 1] );
              if ~isempty(reassign_param)
                  std_mrg = cls(optionalcls(reassign_param{1})-1);
                  name = reassign_param{2};
                  allleaves = 1;
                  N = userdat{2};
                  for k = 1: N %check if all leaves
                      if ~isempty(STUDY.channel(cls(k)).child) 
                          allleaves = 0;
                      end
                  end                     
                  [STUDY] = std_mergechan(STUDY, ALLEEG, std_mrg, name); 
                  % 
                  % update Study history
                  % 
                  if isempty(name)
                      a = ['STUDY = std_mergechan(STUDY, ALLEEG, [' num2str(std_mrg) ']);'];
                  else
                      a = ['STUDY = std_mergechan(STUDY, ALLEEG, [' num2str(std_mrg) '], ' name ');'];
                  end                  
                  STUDY.history =  sprintf('%s\n%s',  STUDY.history, a);
                  userdat{1}{2} = STUDY;
                  %
                  % Replace the merged channels with the one new merged channel 
                  % in the GUI if all channels are leaves
                  %
                  if allleaves                      
                    %
                    % Update channel list
                    %
                    chan_names{end+1} = [STUDY.channel(end).name ' (' num2str(length(STUDY.channel(end).onechan))  ' ICs)']; 
                    %
                    % update the channel list with the new channel
                    %
                    chan_names([optionalcls(reassign_param{1})]) = [];
                    cls = setdiff(cls, std_mrg); % remove from the GUI channels list the merged channels
                    cls(end+1) = length(STUDY.channel); % update the GUI channels list with the new channel
                    N  = length(cls);
                    %
                    % update the first option on the GUI list : 'All 10 channel centroids'
                    % with the new number of channel centroids
                    %
                    ti = strfind(chan_names{1},'channel'); %get the number of channels centroid 
                    cent = num2str(str2num(chan_names{1}(5:ti-2))+1- length(std_mrg)); % new number of centroids
                    chan_names{1} = [chan_names{1}(1:4) cent chan_names{1}(ti-1:end)]; %update list
                    set(findobj('parent', hdl, 'tag', 'chan_list'), 'String', chan_names);
                    %
                    % update Study history
                    %
                    userdat{2} = N; % update N, the number of channel options in edit window 
                    userdat{1}{3} = cls;  % update cls, the channel indices in edit window
                  end
                  set(hdl, 'userdat',userdat); %update information (STUDY)     
                  pop_chanplot('showchan',hdl);
              end
    end
end

function newname = renamechan(oldname, newname);
    
    tmpname = deblank(oldname(end:-1:1));
    strpos  = strfind(oldname, tmpname(end:-1:1));
    
    newname = [ oldname(1:strpos-1) newname ];

% --------------------------
% create groups for channels
% --------------------------
function STUDY = std_changroup(STUDY, ALLEEG);

% union of all channel structures
% -------------------------------
alllocs = ALLEEG(STUDY.datasetinfo(1).index).chanlocs;
alllabs = { alllocs.labels };
for index = 2:length(STUDY.datasetinfo)
   tmplocs = ALLEEG(STUDY.datasetinfo(index).index).chanlocs;
   alllocs = eeg_mergechan(alllocs, tmplocs);
end;

% create group for each electrode
% -------------------------------
for indc = 1:length(alllocs)
    STUDY.changrp(indc).name = [ alllocs(indc).labels ];
    STUDY.changrp(indc).channels = { alllocs(indc).labels };
    tmp = std_chanlookup( STUDY, ALLEEG, STUDY.changrp(indc));
    STUDY.changrp(indc).chaninds = tmp.chaninds;
    STUDY.changrp(indc).centroid = [];
end;
%STUDY.changrp(indc).name = [ 'full montage' ];
%STUDY.changrp(indc).channels = { alllocs.labels };
%tmp = std_chanlookup( STUDY, ALLEEG, STUDY.changrp(indc));
%STUDY.changrp(indc).chaninds = tmp.chaninds;

% ---------------
% channel look-up
% ---------------
function changrp = std_chanlookup( STUDY, ALLEEG, changrp);

    changrp.chaninds = [];
    changrp.chaninds = zeros(1,length(STUDY.datasetinfo));
    for ir = 1:length(STUDY.datasetinfo)
        datind  = STUDY.datasetinfo(ir).index;
        tmplocs = { ALLEEG(STUDY.datasetinfo(datind).index).chanlocs.labels };

        for indc = 1:length(changrp.channels)
            ind = strmatch( changrp.channels(indc), tmplocs, 'exact');
            if ~isempty(ind)
                 changrp.chaninds(ir) = ind;
            else changrp.chaninds(ir) = NaN;
            end;
        end;
    end;
