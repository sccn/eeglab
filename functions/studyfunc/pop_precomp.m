% pop_precomp() - precompute measures (spectrum, ERP, ERSP) for a collection of data
%                 channels.  Calls std_precomp().
% Usage:    
%                >> [STUDY, ALLEEG] = pop_precomp(STUDY, ALLEEG); % pop up interactive window
% Inputs:
%   STUDY        - STUDY set structure containing (loaded) EEG dataset structures
%   ALLEEG       - ALLEEG vector of EEG structures, else a single EEG dataset.
%
% Outputs:
%   STUDY        - the input STUDY set with added pre-clustering data for use by pop_clust() 
%   ALLEEG       - the input ALLEEG vector of EEG dataset structures modified by adding 
%                  pre-clustering data (pointers to .mat files that hold cluster measure information).
%
% Authors: Arnaud Delorme, CERCO, CNRS, 2006-
%
% See also: std_precomp()

% Copyright (C) Arnaud Delorme, CERCO, CNRS, arno@salk.edu
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

% $Log: not supported by cvs2svn $
% Revision 1.11  2007/08/09 18:50:11  arno
% remove test of channel labels
%
% Revision 1.10  2007/08/09 18:44:40  arno
% fix function if no channel labels for datasets
%
% Revision 1.9  2007/04/06 21:46:26  arno
% remove debug msg
%
% Revision 1.8  2007/04/06 21:45:39  arno
% fix channel selection
%
% Revision 1.7  2007/04/06 21:36:01  arno
% fix testing ERSP and ITC
%
% Revision 1.6  2007/02/28 12:02:21  arno
% recompute tag
%
% Revision 1.4  2007/01/26 18:02:55  arno
% warning message
%
% Revision 1.3  2006/11/22 21:05:56  arno
% fixing time limits
%
% Revision 1.2  2006/11/08 01:49:00  arno
% remove newtimef
%
% Revision 1.1  2006/09/12 18:43:23  arno
% Initial revision
%

function [STUDY, ALLEEG, com] = pop_precomp(varargin)

com = '';

if ~isstr(varargin{1}) %intial settings
    if length(varargin) < 2
        error('pop_precomp(): needs both ALLEEG and STUDY structures');
    end
    STUDY  = varargin{1};
    ALLEEG = varargin{2};
    comps = false;
    if nargin > 2
        if strcmpi(varargin{3}, 'components')
            comps = true;
        end;
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
         
    % callbacks
    % ---------
    erspparams_str = [ '''cycles'', [3 0.5], ''padratio'', 1' ];
    specparams_str = '';
    set_ersp       = ['pop_precomp(''setersp'',gcf);']; 
    test_ersp      = ['pop_precomp(''testersp'',gcf);']; 
    set_itc        = ['pop_precomp(''setitc'',gcf);']; 
    set_spec       = ['pop_precomp(''setspec'',gcf);']; 
    test_spec      = ['pop_precomp(''testspec'',gcf);']; 
    str_name       = ['Pre-compute channel measures for STUDY ''' STUDY.name '''' ];
    chanlist       = ['pop_precomp(''chanlist'',gcf);']; 
    chaneditbox    = ['pop_precomp(''chaneditbox'',gcf);']; 
    warninterp     = ['warndlg2(''EEGLAB has not yet been tested if some channels are missing; use at your own risk'');' ];
    
    if comps == true
        str_name       = ['Pre-compute component measures for STUDY ''' STUDY.name '''' ];
        guiadd1 = { {'style' 'checkbox'   'string' '' 'tag' 'compallersp' 'value' 0 }  ...
                    {'style' 'text'       'string' 'Compute ERP/spectrum/ERSP for all components (set) or only those selected by RV (unset)' } };
        guiadd2 = { {'style' 'checkbox'   'string' '' 'tag' 'scalp_on' 'value' 0 }  ...
                    {'style' 'text'       'string' 'Scalp maps' } };
        geomadd1 = { [0.33 6] };
        geomadd2 = { [0.33 6] };
    else
        str_name       = ['Pre-compute channel measures for STUDY ''' STUDY.name '''' ];
        guiadd1 = { {'style' 'text'       'string' 'Channel list (default:all)' 'FontWeight' 'Bold'} ...
            {'Style' 'edit'       'string' '' 'tag' 'chans' 'callback' chaneditbox }, ...
            {'style' 'pushbutton' 'string'  '...', 'enable' fastif(isempty(ALLEEG(1).chanlocs), 'off', 'on') ...
            'callback' chanlist }, ...
            {'style' 'checkbox'   'string' '' 'tag' 'interpolate_on' 'value' 1 'callback' warninterp }  ...
            {'style' 'text'       'string' 'Interpolate missing channels (datasets will be modified on disk)' } };
        guiadd2 = {};
        geomadd1 = { [2 3 0.5] [0.33 6] }; 
        geomadd2 = { };
    end;
            
    gui_spec = { ...
    {'style' 'text'       'string' str_name 'FontWeight' 'Bold' 'horizontalalignment' 'left'} {} ...
    guiadd1{:} ...
    {} {'style' 'text'    'string' 'List of measures to precompute' 'FontWeight' 'Bold' 'horizontalalignment' 'left'} ...
    {'style' 'checkbox'   'string' '' 'tag' 'erp_on' 'value' 0 }  ...
	{'style' 'text'       'string' 'ERPs' } ...
    {'style' 'checkbox'   'string' '' 'tag' 'spectra_on' 'value' 0 'Callback' set_spec } ...
	{'style' 'text'       'string' 'Power spectrum' } {} ...
    {'style' 'text'       'string' 'Parameters' 'tag' 'spec_push' 'value' 1 'enable' 'off'}...
    {'style' 'edit'       'string' specparams_str 'tag' 'spec_params' 'enable' 'off' } ...
    {'style' 'pushbutton' 'string' 'Test' 'tag' 'spec_test' 'enable' 'off' 'callback' test_spec}...
    {'style' 'checkbox'   'string' '' 'tag' 'ersp_on' 'value' 0 'Callback' set_ersp }  ...
	{'style' 'text'       'string' 'ERSPs' 'horizontalalignment' 'center' } {} ...
    {'vertshift' 'style'  'text'       'string' 'Time/freq. parameters' 'tag' 'ersp_push' 'value' 1 'enable' 'off'} ...
    {'vertshift' 'style'  'edit'       'string' erspparams_str 'tag' 'ersp_params' 'enable' 'off'}...
    {'vertshift' 'style'  'pushbutton' 'string' 'Test' 'tag' 'ersp_test' 'enable' 'off' 'callback' test_ersp }...
    {'style' 'checkbox'   'string' '' 'tag' 'itc_on' 'value' 0 'Callback' set_itc } ...
	{'style' 'text'       'string' 'ITCs' 'horizontalalignment' 'center' } {'link2lines' 'style'  'text'   'string' '' } {} {} {} ...
    guiadd2{:} ...
    {} {'style' 'checkbox'   'string' '' 'tag' 'recomp_on' 'value' 0 } ...
	{'style' 'text'       'string' 'Recompute even if present on disk' }
    };
  
	%{'style' 'checkbox'   'string' '' 'tag' 'precomp_PCA'  'Callback' precomp_PCA 'value' 0} ...
	%{'style' 'text'       'string' 'Do not prepare dataset for clustering at this time.' 'FontWeight' 'Bold'  } {} ...

    % find the list of all channels
    % -----------------------------
    allchans  = { };
    keepindex = 0;
    for index = 1:length(ALLEEG)
        tmpchans = { ALLEEG(index).chanlocs.labels };
        allchans = unique({ allchans{:} ALLEEG(index).chanlocs.labels });
        if length(allchans) == length(tmpchans), keepindex = index; end;
    end;
    if keepindex, allchans = { ALLEEG(keepindex).chanlocs.labels }; end;
    
    chanlist = {};
    firsttimeersp = 1;
    fig_arg = { ALLEEG STUDY allchans chanlist firsttimeersp };
    geomline = [0.45 1 0.3 2 3 0.7 ];
    geometry = { [1] [1] geomadd1{:}  [1] [1] [0.33 6] [0.45 1.5 0.3 1.5 3 0.7 ] geomline geomline geomadd2{:} 1 [0.05 1] };
    geomvert = [ 1 0.5 fastif(length(geomadd1) == 1,1,[1 1]) 0.5 1 1 1 1 1 1 fastif(length(geomadd2) == 1,1,[]) 1 1];
	[precomp_param, userdat2, strhalt, os] = inputgui( 'geometry', geometry, 'uilist', gui_spec, 'geomvert', geomvert, ...
                                                      'helpcom', ' pophelp(''std_precomp'')', ...
                                                      'title', 'Select and compute component measures for later clustering -- pop_precomp()', ...
                                                      'userdata', fig_arg);	
	if isempty(precomp_param), return; end;
    
    if comps == 1
        options = { STUDY ALLEEG 'components' };
    else
        options = { STUDY ALLEEG userdat2{4} };
    end
    if ~isfield(os, 'interpolate_on'), os.interpolate_on = 0; end;
    if ~isfield(os, 'scalp_on'),    os.scalp_on = 0; end;
    if ~isfield(os, 'compallersp'), os.compallersp = 0; end;
    
    % interpolate option is on
    % ------------------------
    if os.interpolate_on == 1 
        options = { options{:} 'interpolate' 'on' };
    end

    % compallersp option is on
    % ------------------------
    if os.compallersp == 1 
        options = { options{:} 'allcomps' 'on' };
    end    
    
    % recompute option is on
    % ----------------------
    if os.recomp_on == 1 
        options = { options{:} 'recompute' 'on' };
    end
    
    % ERP option is on
    % ----------------
    if os.erp_on == 1 
        options = { options{:} 'erp' 'on' };
    end
    
    % SCALP option is on
    % ----------------
    if os.scalp_on == 1 
        options = { options{:} 'scalp' 'on' };
    end
    
    % Spectrum option is on
    % --------------------
    if os.spectra_on== 1 
        tmpparams = eval( [ '{' os.spec_params '}' ] );
        options = { options{:} 'spec' 'on' 'specparams' tmpparams };
    end
    
    % ERSP option is on
    % -----------------
    if os.ersp_on  == 1 
        tmpparams = eval( [ '{' os.ersp_params '}' ] );
        options = { options{:} 'ersp' 'on' 'erspparams' tmpparams };
    end
    
    % ITC option is on 
    % ----------------
    if os.itc_on  == 1 
        tmpparams = eval( [ '{' os.ersp_params '}' ] );
        options = { options{:} 'itc' 'on' };
        if os.ersp_on  == 0, options = { options{:} 'erspparams' tmpparams }; end;
    end       
        
    % evaluate command
    % ----------------
    if length(options) == 4
        warndlg2('No measure selected: aborting.'); 
        return; 
    end;
    [STUDY ALLEEG] = std_precomp(options{:});
    com = sprintf('%s\n[STUDY ALLEEG] = std_precomp(STUDY, ALLEEG, %s);', ...
                  STUDY.history, vararg2str(options(3:end)));
    
else
    hdl = varargin{2}; %figure handle
    userdat = get(varargin{2}, 'userdata');    
    ALLEEG  = userdat{1};
    STUDY   = userdat{2};
    allchans = userdat{3};
    chansel  = userdat{4};
    firsttimeersp = userdat{5};

    switch  varargin{1}
               
        case 'chanlist'
            [tmp tmp2 tmp3] = pop_chansel(allchans, 'select', chansel);
            set(findobj('parent', hdl, 'tag', 'chans'), 'string', tmp2);
            userdat{4} = tmp3;
            set(hdl, 'userdata',userdat); 
     
        case 'chaneditbox'
            userdat{4} = parsetxt(get(findobj('parent', hdl, 'tag', 'chans'), 'string'));
            set(hdl, 'userdata',userdat); 
            
        case { 'setitc' 'setersp' }
            if firsttimeersp
                warndlg2(strvcat('Checking both ''ERSP'' and ''ITC'' does not require further', ...
                                 'computing time. However it requires disk space'));
            end;
            set_itc  = get(findobj('parent', hdl, 'tag', 'itc_on'), 'value'); 
            set_ersp = get(findobj('parent', hdl, 'tag', 'ersp_on'), 'value'); 
            if  (~set_ersp & ~set_itc )
                set(findobj('parent', hdl,'tag', 'ersp_push'),   'enable', 'off');
                set(findobj('parent', hdl,'tag', 'ersp_params'), 'enable', 'off');
                set(findobj('parent', hdl,'tag', 'ersp_test'),   'enable', 'off');                
            else
                set(findobj('parent', hdl,'tag', 'ersp_push'),   'enable', 'on');
                set(findobj('parent', hdl,'tag', 'ersp_params'), 'enable', 'on');                
                set(findobj('parent', hdl,'tag', 'ersp_test'),   'enable', 'on');                
            end
            userdat{5} = 0;
            set(hdl, 'userdata',userdat); 
            
        case 'setspec'
            set_spec = get(findobj('parent', hdl, 'tag', 'spectra_on'), 'value'); 
            if set_spec
                 set(findobj('parent', hdl,'tag', 'spec_push'),   'enable', 'on');
                 set(findobj('parent', hdl,'tag', 'spec_params'), 'enable', 'on');
                 set(findobj('parent', hdl,'tag', 'spec_test'),   'enable', 'on');
            else set(findobj('parent', hdl,'tag', 'spec_push'),   'enable', 'off');
                 set(findobj('parent', hdl,'tag', 'spec_params'), 'enable', 'off');
                 set(findobj('parent', hdl,'tag', 'spec_test'),   'enable', 'off');
            end

        case 'testspec'
            try,
                spec_params = eval([ '{' get(findobj('parent', hdl, 'tag', 'spec_params'), 'string') '}' ]); 

                EEG = eeg_checkset(ALLEEG(1), 'loaddata');
                data = EEG.data(1,:,1:max(EEG.trials,10));
                figure; spectopo( data, EEG.pnts, EEG.srate, spec_params{:});
                tmplim = ylim;
                text( EEG.srate/4, mean(tmplim), strvcat('This is a test plot performed on', ...
                                                         'the first 10 trials of the first', ....
                                                         'dataset to determine if you like', ...
                                                         'the frequency resolution.', ...
                                                         'Frequency range may be adjusted', ...
                                                         'after computation'));
            catch, warndlg2('Error while calling function, check parameters'); end;

        case 'testersp'
            %try,
                ersp_params = eval([ '{' get(findobj('parent', hdl, 'tag', 'ersp_params'), 'string') '}' ]); 
                tmpstruct = struct(ersp_params{:});
                [ tmpX tmpt tmpf ersp_params ] = std_ersp(ALLEEG(1), 'channels', 1, 'type', 'ersp', 'recompute', 'on', 'getparams', 'on', ersp_params{:});

                set_itc  = get(findobj('parent', hdl, 'tag', 'itc_on'), 'value'); 
                set_ersp = get(findobj('parent', hdl, 'tag', 'ersp_on'), 'value'); 
                opt = {};
                for ind = length(ersp_params)-1:-2:1, 
                    if strcmpi(ersp_params{ind}, 'plotitc'),  ersp_params(ind:ind+1) = []; end;
                    if strcmpi(ersp_params{ind}, 'plotersp'), ersp_params(ind:ind+1) = []; end;
                end;
                if ~set_itc,  opt = { opt{:} 'plotitc',  'off' }; end;
                if ~set_ersp, opt = { opt{:} 'plotersp', 'off' }; end;

                EEG = eeg_checkset(ALLEEG(1), 'loaddata');
                data = EEG.data(1,:,1:min(EEG.trials,10));
                %f1   = fieldnames( ersp_params);
                %f2   = struct2cell(ersp_params);
                %ersp_params = {f1{:}; f2{:} }; ersp_params = ersp_params(:)';
                
                figure; pos = get(gcf, 'position'); pos(3)=pos(3)*2; set(gcf, 'position', pos);
                subplot(1,2,1); newtimef( data, EEG.pnts, [ EEG.xmin EEG.xmax ]*1000, EEG.srate, tmpstruct.cycles, opt{:}, 'maxfreq', EEG.srate/2, ersp_params{:});
                subplot(1,2,2); 
                text( 0.2, 0.8, strvcat('This is a test plot performed on', ...
                                                         'the first 10 trials of the first', ....
                                                         'dataset to determine if you like', ...
                                                         'the time and frequency resolution.', ...
                                                         ' ', ...
                                                         'Time and frequency range may also be', ...
                                                         'adjusted after computation.'));
                axis off;
            %catch, warndlg2(strvcat('Error while calling function, check parameters', lasterr)); 
            %end;
                                                 
    end;
end
STUDY.saved = 'no';
