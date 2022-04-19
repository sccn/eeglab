% pop_runica() - Run an ICA decomposition of an EEG dataset using runica(), 
%                binica(), or another ICA or other linear decomposition. 
% Usage:
%   >> OUT_EEG = pop_runica( EEG ); % pops-up a data entry window
%   >> OUT_EEG = pop_runica( EEG, 'key', 'val' ); % no pop_up
%
% Graphic interface:
%   "ICA algorithm to use" - [edit box] The ICA algorithm to use for 
%                 ICA decomposition. Command line equivalent: 'icatype'
%   "Commandline options" - [edit box] Command line options to forward
%                 to the ICA algorithm. Command line equivalent: 'options' 
% Inputs:
%   EEG         - input EEG dataset or array of datasets
%
% Optional inputs:
%   'icatype'   - ['runica'|'binica'|'jader'|'fastica'] ICA algorithm 
%                 to use for the ICA decomposition. The nature of any 
%                 differences in the results of these algorithms have 
%                 not been well characterized. {default: binica(), if
%                 found, else runica()}
%   'dataset'   - [integer array] dataset index or indices.
%   'chanind'   - [integer array or cell array] subset of channel indices 
%                 for running the ICA decomposition. Alternatively, you may
%                 also enter channel types here in a cell array.
%   'concatenate' - ['on'|'off'] 'on' concatenate all input datasets 
%                 (assuming there are several). 'off' run ICA independently
%                 on each dataset. Default is 'off'.
%   'concatcond'  - ['on'|'off'] 'on' concatenate conditions for input datasets 
%                 of the same sessions and the same subject. Default is 'off'.
%   'reorder'   - ['on'|'off'] re-order components by variance if that's not
%                 already the case. Default is 'on'.
%   'key','val' - ICA algorithm options (see ICA routine help messages).
% 
% Adding a new algorithm:
% Add the algorithm to the list of algorithms line 366 to 466, for example
%
%    case 'myalgo', [EEG.icaweights] = myalgo( tmpdata, g.options{:} );
%
% where "myalgo" is the name of your algorithm (and Matlab function). 
% tmpdata is the 2-D array containing the EEG data (channels x points) and
% g.options{} contains custom options for your algorithm (there is no
% predetermined format for these options). The output EEG.icaweights is the
% mixing matrix (or inverse of the unmixing matrix).
%
% Note:
% 1) Infomax (runica, binica) is the ICA algorithm we use most. It is based 
%    on Tony Bell's infomax algorithm as implemented for automated use by 
%    Scott Makeig et al. using the natural gradient of Amari et al. It can 
%    also extract sub-Gaussian sources using the (recommended) 'extended' option 
%    of Lee and Girolami. Function runica() is the all-Matlab version; function 
%    binica() calls the (1.5x faster) binary version (a separate download) 
%    translated into C from runica() by Sigurd Enghoff.
% 2) jader() calls the JADE algorithm of Jean-Francois Cardoso. This is 
%    included in the EEGLAB toolbox by his permission. See >> help jader
% 3) To run fastica(), download the fastICA toolbox from its website,
%    http://www.cis.hut.fi/projects/ica/fastica/, and make it available 
%    in your Matlab path. According to its authors, default parameters
%    are not optimal: Try args 'approach', 'sym' to estimate components 
%    in parallel.
%
% Outputs:
%   OUT_EEG = The input EEGLAB dataset with new fields icaweights, icasphere 
%             and icachansind (channel indices). 
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: runica(), binica(), jader(), fastica()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% 01-25-02 reformated help & license -ad 
% 03-07-02 add the eeglab options -ad
% 03-18-02 add other decomposition options -ad
% 03-19-02 text edition -sm

function [ALLEEG, com] = pop_runica( ALLEEG, varargin )

com = '';
if nargin < 1   
    help pop_runica;
    return;
end

% find available algorithms
% -------------------------
allalgs   = { 'runica' 'binica' 'jader' 'jadeop' 'jade_td_p' 'MatlabshibbsR' 'fastica' ...
              'tica' 'erica' 'simbec' 'unica' 'amuse' 'fobi' 'evd' 'evd24' 'sons' 'sobi' 'ng_ol' ...
              'acsobiro' 'acrsobibpf' 'pearson_ica' 'egld_ica' 'eeA' 'tfbss' 'icaML' 'icaMS' 'picard' }; % do not use egld_ica => too slow
selectalg = {};
linenb    = 1;
count     = 1;
for index = length(allalgs):-1:1
    if exist(allalgs{index}) ~= 2 && exist(allalgs{index}) ~= 6
        allalgs(index) = [];
    end
end

% special AMICA
% -------------
selectamica = 0;
defaultopts = [ '''extended'', 1' ] ;
if nargin > 1
    if ischar(varargin{1})
        if strcmpi(varargin{1}, 'selectamica')
            selectamica = 1;
            allalgs = { 'amica' allalgs{:} };
            defaultopts = sprintf('''outdir'', ''%s''', fullfile(pwd, 'amicaout'));
        elseif strcmpi(varargin{1}, 'selectamicaloc')
            selectamica = 1;
            allalgs = { 'amica' allalgs{:} };
            defaultopts = sprintf('''outdir'', ''%s'', ''qsub'', ''off''', fullfile(pwd, 'amicaout'));
        end
    end
end

% popup window parameters
% -----------------------
fig = [];
if nargin < 2 || selectamica
    commandchans = [ 'tmpchans = get(gcbf, ''userdata'');' ...
                     'tmpchans = tmpchans{1};' ...
                     'set(findobj(gcbf, ''tag'', ''chantype''), ''string'', ' ...
                     '       int2str(pop_chansel( tmpchans )));' ...
                     'clear tmpchans;' ];
    commandtype = ['tmptype = get(gcbf, ''userdata'');' ...
                   'tmptype = tmptype{2};' ...
                   'if ~isempty(tmptype),' ...
                   '    [tmps,tmpv, tmpstr] = listdlg2(''PromptString'',''Select type(s)'', ''ListString'', tmptype);' ...
				   '    if tmpv' ...
				   '        set(findobj(''parent'', gcbf, ''tag'', ''chantype''), ''string'', tmpstr);' ...
				   '    end;' ...
                   'else,' ...
                   '    warndlg2(''No channel type'', ''No channel type'');' ...
                   'end;' ...
				   'clear tmps tmpv tmpstr tmptype tmpchans;' ];
    cb_ica = [ 'if get(gcbo, ''value'') < 3, ' ...
               '     set(findobj(gcbf, ''tag'', ''params''), ''string'', ''''''extended'''', 1'');' ...
               'else,' ...
               '    tmpStr =  get(gcbo, ''string'');' ....
               '    tmpAlgo = tmpStr(get(gcbo, ''value''));' ...
               '    if strcmpi(tmpAlgo, ''picard''),' ...
               '        set(findobj(gcbf, ''tag'', ''params''), ''string'', ''''''maxiter'''', 500'');' ...
               '    else,' ...
               '        set(findobj(gcbf, ''tag'', ''params''), ''string'', '''');' ...
               '    end;' ...
               'end;' ];
               
    promptstr    = { { 'style' 'text'       'string' 'ICA algorithm to use (click to select)' } ...
                     { 'style' 'listbox'    'string' allalgs 'callback', cb_ica } ...
                     { 'style' 'text'       'string' 'Commandline options (See help messages)' } ...
                     { 'style' 'edit'       'string' defaultopts 'tag' 'params' } ...
                     { 'style'  'checkbox'  'string' 'Reorder components by variance (if that''s not already the case)' 'value' 1 } ...
                     { 'style' 'text'       'string' 'Channel type(s) or channel indices' } ...
                     { 'style' 'edit'       'string' '' 'tag' 'chantype' }  ...
                     { 'style' 'pushbutton' 'string' '... types' 'callback' commandtype } ...
                     { 'style' 'pushbutton' 'string' '... channels' 'callback' commandchans } };
    geometry = { [2 1.5] [2 1.5] [1] [2 1 1 1] };
    geomvert = [ 1.5 1 1 1];
    if length(ALLEEG) > 1
        cb1 = 'set(findobj(''parent'', gcbf, ''tag'', ''concat2''), ''value'', 0);';
        cb2 = 'set(findobj(''parent'', gcbf, ''tag'', ''concat1''), ''value'', 0);';
        promptstr = { promptstr{:}, ...
                     { 'style' 'text'       'string' 'Concatenate all datasets (check=yes; uncheck=run ICA on each dataset)?' }, ...
                     { 'style' 'checkbox'   'string' '' 'value' 0 'tag' 'concat1' 'callback' cb1 }, ...
                     { 'style' 'text'       'string' 'Concatenate datasets for the same subject and session (check=yes)?' }, ...
                     { 'style' 'checkbox'   'string' '' 'value' 1 'tag' 'concat2' 'callback' cb2 } };
        geometry = { geometry{:} [ 2 0.2 ] [ 2 0.2 ]};
        geomvert = [ geomvert 1 1];
    end
    
    % channel types
    % -------------
    if isfield(ALLEEG(1).chanlocs, 'type')
        tmpchanlocs = ALLEEG(1).chanlocs;
        alltypes = { tmpchanlocs.type };
        indempty = cellfun('isempty', alltypes);
        alltypes(indempty) = '';
        try
            alltypes = unique_bc(alltypes);
        catch
            alltypes = '';
        end
    else
        alltypes = '';
    end
    
    % channel labels
    % --------------
    if ~isempty(ALLEEG(1).chanlocs)
        tmpchanlocs = ALLEEG(1).chanlocs;        
        alllabels = { tmpchanlocs.labels };
    else
        for index = 1:ALLEEG(1).nbchan
            alllabels{index} = int2str(index);
        end
    end
    
    % gui
    % ---
    result       = inputgui( 'geometry', geometry, 'geomvert', geomvert, 'uilist', promptstr, ...
                             'helpcom', 'pophelp(''pop_runica'')', ...
                             'title', 'Run ICA decomposition -- pop_runica()', 'userdata', { alllabels alltypes } );
    if length(result) == 0 return; end      
    options = { 'icatype' allalgs{result{1}} 'dataset' [1:length(ALLEEG)] 'options' eval( [ '{' result{2} '}' ]) 'reorder' fastif(result{3}, 'on', 'off') };
    if ~isempty(result{4})
        if ~isempty(str2num(result{4})), options = { options{:} 'chanind' str2num(result{4}) };
        else                             options = { options{:} 'chanind' parsetxt(result{4}) }; 
        end
    end
    if length(result) > 4
        options = { options{:} 'concatenate' fastif(result{5}, 'on', 'off') };
        options = { options{:} 'concatcond'  fastif(result{6}, 'on', 'off') };
    end
else 
    if mod(length(varargin),2) == 1
        options = { 'icatype' varargin{1:end} };
    else
        options = varargin;
    end
end

% decode input arguments
% ----------------------
[ g, addoptions ] = finputcheck( options, { 'icatype'        'string'  allalgs   'runica'; ...
                            'dataset'        'integer' []        [1:length(ALLEEG)];
                            'options'        'cell'    []        {};
                            'concatenate'    'string'  { 'on','off' }   'off';
                            'concatcond'     'string'  { 'on','off' }   'off';
                            'reorder'        'string'  { 'on','off' }   'on';
                            'chanind'        { 'cell','integer' } { [] [] }        [];}, ...
                            'pop_runica', 'ignore');
if ischar(g), error(g); end
if ~isempty(addoptions), g.options = { g.options{:} addoptions{:}}; end

% select datasets, create new big dataset if necessary
% ----------------------------------------------------
if length(g.dataset) == 1
    EEG = ALLEEG(g.dataset);
    EEG = eeg_checkset(EEG, 'loaddata');
    if isfield(EEG.etc, 'ic_classification')
        EEG.etc = rmfield(EEG.etc, 'ic_classification');
    end
elseif length(ALLEEG) > 1 && ~strcmpi(g.concatenate, 'on') && ~strcmpi(g.concatcond, 'on')
    [ ALLEEG, com ] = eeg_eval( 'pop_runica', ALLEEG, 'warning', 'off', 'params', ...
           { 'icatype' g.icatype 'options' g.options 'chanind' g.chanind } );
    return;
elseif length(ALLEEG) > 1 && strcmpi(g.concatcond, 'on')
    allsubjects = { ALLEEG.subject };
    allsessions = { ALLEEG.session };
    alltags     = zeros(1,length(allsubjects));
    if any(cellfun('isempty', allsubjects))
        errordlg2( [ 'Aborting: Subject names missing from at least one dataset file.' 10 ...
                     'Subject names must be stored within the datasets. To do so,' 10 ...
                     'use the STUDY > Edit STUDY Info menu and check the box' 10 ...
                     '"Dataset info (condition, group, ...) differs from study info..."' ]);
    end
    dats = {};
    for index = 1:length(allsubjects)
        if ~alltags(index)
            allinds = strmatch(allsubjects{index}, allsubjects, 'exact');
            rmind = [];
            % if we have different sessions they will not be concatenated
            for tmpi = setdiff_bc(allinds,index)'
                if ~isequal(allsessions(index), allsessions(tmpi))
                    rmind = [rmind tmpi];
                end
            end
            allinds = setdiff_bc(allinds, rmind);
            fprintf('Found %d datasets for subject ''%s'' session %d\n', length(allinds), allsubjects{index}, allsessions{index});
            dats = { dats{:} allinds };
            alltags(allinds) = 1;
        end
    end
    fprintf('**************************\nNOW RUNNING ALL DECOMPOSITIONS\n****************************\n');
    eeglab_options;
    if option_parallel
        for index = 1:length(dats)
            TMPEEG{index} = ALLEEG(dats{index});
        end
        parfor index = 1:length(dats)
            TMPEEG{index} = pop_runica(TMPEEG{index}, 'icatype', g.icatype, ...
                'options', g.options, 'chanind', g.chanind, 'concatenate', 'on');
        end
        for index = 1:length(dats)
            ALLEEG(dats{index}) = TMPEEG{index};
            for idat = 1:length(dats{index})
                ALLEEG(dats{index}(idat)).saved = 'no';
                pop_saveset(ALLEEG(dats{index}(idat)), 'savemode', 'resave');
                ALLEEG(dats{index}(idat)).saved = 'yes';
            end
        end
    else
        for index = 1:length(dats)
            ALLEEG(dats{index}) = pop_runica(ALLEEG(dats{index}), 'icatype', g.icatype, ...
                'options', g.options, 'chanind', g.chanind, 'concatenate', 'on');
            for idat = 1:length(dats{index})
                ALLEEG(dats{index}(idat)).saved = 'no';
                pop_saveset(ALLEEG(dats{index}(idat)), 'savemode', 'resave');
                ALLEEG(dats{index}(idat)).saved = 'yes';
                ALLEEG(dats{index}(idat)) = update_datafield(ALLEEG(dats{index}(idat)));
            end
        end
    end
    com = sprintf('EEG = pop_runica(EEG, %s);', ...
              vararg2str({ 'icatype' g.icatype 'concatcond' 'on' 'options' g.options }) );
    return;
else
    disp('Concatenating datasets...');
    EEG = ALLEEG(g.dataset(1));
    
    if isfield(EEG.etc, 'ic_classification')
        EEG.etc = rmfield(EEG.etc, 'ic_classification');
    end
    
    % compute total data size
    % -----------------------
    totalpnts = 0;
    for i = g.dataset
        totalpnts = totalpnts+ALLEEG(g.dataset(i)).pnts*ALLEEG(g.dataset(i)).trials;
    end
    EEG.data = zeros(EEG.nbchan, totalpnts);
    
    % copy data
    % ---------
    cpnts = 1;
    for i = g.dataset
        tmplen = ALLEEG(g.dataset(i)).pnts*ALLEEG(g.dataset(i)).trials;
        TMP = eeg_checkset(ALLEEG(g.dataset(i)), 'loaddata');
        EEG.data(:,cpnts:cpnts+tmplen-1) = reshape(TMP.data, size(TMP.data,1), size(TMP.data,2)*size(TMP.data,3));
        cpnts = cpnts+tmplen;
    end
    EEG.icaweights = [];
    EEG.trials = 1;
    EEG.pnts   = size(EEG.data,2);
    EEG.saved  = 'no';
end

% Store and then remove current EEG ICA weights and sphere
% ---------------------------------------------------
fprintf('\n');
if ~isempty(EEG.icaweights)
    fprintf('Saving current ICA decomposition in "EEG.etc.oldicaweights" (etc.).\n');
    if ~isfield(EEG,'etc'), EEG.etc = []; end
    if ~isfield(EEG.etc,'oldicaweights')
        EEG.etc.oldicaweights = {};
        EEG.etc.oldicasphere = {};
        EEG.etc.oldicachansind = {};
    end
    tmpoldicaweights  = EEG.etc.oldicaweights;
    tmpoldicasphere   = EEG.etc.oldicasphere;
    tmpoldicachansind = EEG.etc.oldicachansind;
    EEG.etc.oldicaweights = { EEG.icaweights    tmpoldicaweights{:} };
    EEG.etc.oldicasphere  = { EEG.icasphere     tmpoldicasphere{:}  };
    EEG.etc.oldicachansind  = { EEG.icachansind tmpoldicachansind{:}  };
    fprintf('               Decomposition saved as entry %d.\n',length(EEG.etc.oldicaweights));
end
EEG.icaweights = [];
EEG.icasphere  = [];
EEG.icawinv    = [];
EEG.icaact     = [];

% select sub_channels
% -------------------
if isempty(g.chanind)
    g.chanind = 1:EEG.nbchan;
end
if iscell(g.chanind)
    datatype = {EEG.chanlocs.type};
    tmpChanInd = [];
    for iChan = 1:length(datatype)
        if ~isempty(datatype{iChan}) && ~isempty(strmatch(datatype{iChan}, g.chanind))
            tmpChanInd = [ tmpChanInd iChan ];
        end
    end
    g.chanind = tmpChanInd;
end
EEG.icachansind = g.chanind;

% is pca already an option?
% -------------------------
pca_opt = 0;
for i = 1:length(g.options)
    if ischar(g.options{i})
        if strcmpi(g.options{i}, 'pca')
            pca_opt = 1;
            pca_ind = i;
        end
    end
end
if pca_opt
    fprintf([ 'Warning: you have used PCA to reduce dimensionality so ICA\n' ...
        '         is not modeling the entire data, only the PCA-reduced data.\n' ]);
end

%------------------------------
% compute ICA on a definite set
% -----------------------------
tmpdata = reshape( EEG.data(g.chanind,:,:), length(g.chanind), EEG.pnts*EEG.trials);
tmprank = getrank(double(tmpdata(:,1:min(3000, size(tmpdata,2)))), pca_opt);
tmpdata = tmpdata - repmat(mean(tmpdata,2), [1 size(tmpdata,2)]); % zero mean 
if ~strcmpi(g.icatype, 'binica')
    try
        disp('Attempting to convert data matrix to double precision for more accurate ICA results.')
        tmpdata = double(tmpdata);
        tmpdata = tmpdata - repmat(mean(tmpdata,2), [1 size(tmpdata,2)]); % zero mean (more precise than single precision)
    catch
        disp('*************************************************************')
        disp('Not enough memory to convert data matrix to double precision.')
        disp('All computations will be done in single precision. Matlab 7.x')
        disp('under 64-bit Linux and others is imprecise in this mode.')
        disp('We advise use of "binica" instead of "runica."')
        disp('*************************************************************')
    end
end
switch lower(g.icatype)
    case 'runica' 
        try if ismatlab, g.options = {  g.options{:}, 'interrupt', 'on' }; end; catch, end
        if tmprank == size(tmpdata,1) || pca_opt
            [EEG.icaweights,EEG.icasphere] = runica( tmpdata, 'lrate', 0.001,  g.options{:} );
        else 
            if nargin < 2
                uilist = { { 'style' 'text' 'string' [ 'EEGLAB has detected that the rank of your data matrix' 10 ...
                                                       'is lower the number of input data channels. This might' 10 ...
                                                       'be because you are including a reference channel or' 10 ...
                                                       'because you are running a second ICA decomposition.' 10 ...
                                                       sprintf('The proposed dimension for ICA is %d (out of %d channels).', tmprank, size(tmpdata,1)) 10 ...
                                                       'Rank computation may be inaccurate so you may edit this' 10 ...
                                                       'number below. If you do not understand, simply press OK.' ] } { } ...
                           { 'style' 'text' 'string' 'Proposed rank:' } ...
                           { 'style' 'edit' 'string' num2str(tmprank) } };
                res = inputgui('uilist', uilist, 'geometry', { [1] [1] [1 1] }, 'geomvert', [6 1 1]);
                if isempty(res), return; end
                tmprank = str2num(res{1});
                g.options = [g.options  { 'pca' tmprank }];
            else
                g.options = [g.options  {'pca' tmprank }]; % automatic for STUDY (batch processing)
            end
            disp(['Data rank (' int2str(tmprank) ') is smaller than the number of channels (' int2str(size(tmpdata,1)) ').']);
            [EEG.icaweights,EEG.icasphere] = runica( tmpdata, 'lrate', 0.001, g.options{:} );
        end
     case 'binica'
        icadefs;
        fprintf(['Warning: If the binary ICA function does not work, check that you have added the\n' ...
                 'binary file location (in the EEGLAB directory) to your Unix /bin directory (.cshrc file)\n']);
        if exist(ICABINARY) ~= 2
            error('Pop_runica(): binary ICA executable not found. Edit icadefs.m file to specify the ICABINARY location');
        end
        tmprank = getrank(tmpdata(:,1:min(3000, size(tmpdata,2))));
        if tmprank == size(tmpdata,1) || pca_opt
            [EEG.icaweights,EEG.icasphere] = binica( tmpdata, 'lrate', 0.001, g.options{:} );
        else 
            disp(['Data rank (' int2str(tmprank) ') is smaller than the number of channels (' int2str(size(tmpdata,1)) ').']);
            [EEG.icaweights,EEG.icasphere] = binica( tmpdata, 'lrate', 0.001, 'pca', tmprank, g.options{:} );
        end
    case 'amica' 
        tmprank = getrank(tmpdata(:,1:min(3000, size(tmpdata,2))));
        fprintf('Now Running AMICA\n');
        if length(g.options) > 1
            if ischar(g.options{2})
                fprintf('See folder %s for outputs\n', g.options{2});
            end
        end
        fprintf('To import results, use menu item "Tools > Run AMICA > Load AMICA components\n');
        modres = runamica( tmpdata, [], size(tmpdata,1), size(tmpdata,2), g.options{:} );
        if ~isempty(modres)
            EEG.icaweights = modres.W;
            EEG.icasphere  = modres.S;
        else
            return;
        end
     case 'picard' 
         options2 = g.options;
         if pca_opt
             if g.options{pca_ind+1} < 0
                 fprintf('Decreasing data dimension from %d to %d using PCA...\n', size(tmpdata,1), size(tmpdata,1)+g.options{pca_ind+1});
                 [tmpdata,eigvec] = runpca(tmpdata, size(tmpdata,1)+g.options{pca_ind+1});
             else
                 fprintf('Decreasing data dimension from %d to %d using PCA...\n', size(tmpdata,1), g.options{pca_ind+1});
                 [tmpdata,eigvec] = runpca(tmpdata, g.options{pca_ind+1});
             end
             options2(pca_ind:pca_ind+1) = [];
         end
        [~, EEG.icaweights] = picard( tmpdata, 'verbose', true, options2{:});
         if pca_opt
            EEG.icaweights = EEG.icaweights*pinv(eigvec);
         end
     case 'pearson_ica' 
        if isempty(g.options)
            disp('Warning: EEGLAB default for pearson ICA is 1000 iterations and epsilon=0.0005');
            [~, EEG.icaweights] = pearson_ica( tmpdata, 'maxNumIterations', 1000,'epsilon',0.0005);
        else    
            [~, EEG.icaweights] = pearson_ica( tmpdata, g.options{:});
        end
     case 'egld_ica', disp('Warning: This algorithm is very slow!!!');
        [~, EEG.icaweights] = egld_ica( tmpdata, g.options{:} );
     case 'tfbss' 
        if  isempty(g.options)
             [~, EEG.icaweights] = tfbss( tmpdata, size(tmpdata,1), 8, 512 );
        else    
             [~, EEG.icaweights] = tfbss( tmpdata, g.options{:} );
        end
     case 'jader',         [EEG.icaweights] = jader( tmpdata, g.options{:} );
     case 'matlabshibbsr', [EEG.icaweights] = MatlabshibbsR( tmpdata, g.options{:} );
     case 'eea',           [EEG.icaweights] = eeA( tmpdata, g.options{:} );
     case 'icaml',         [~, EEG.icawinv] = icaML( tmpdata, g.options{:} );
     case 'icams',         [~, EEG.icawinv] = icaMS( tmpdata, g.options{:} );
     case 'fastica',       [~, EEG.icawinv, EEG.icaweights] = fastica( tmpdata, 'displayMode', 'off', g.options{:} );
     case { 'tica' 'erica' 'simbec' 'unica' 'amuse' 'fobi' 'evd' 'sons' ...
            'jadeop' 'jade_td_p' 'evd24' 'sobi' 'ng_ol' 'acsobiro' 'acrsobibpf' } 
        fig = figure('tag', 'alg_is_run', 'visible', 'off');
        
        if isempty(g.options), g.options = { size(tmpdata,1) }; end
        switch lower(g.icatype)
         case 'tica',     EEG.icaweights = tica( tmpdata, g.options{:} );
         case 'erica',    EEG.icaweights = erica( tmpdata, g.options{:} );
         case 'simbec',   EEG.icaweights = simbec( tmpdata, g.options{:} );
         case 'unica',    EEG.icaweights = unica( tmpdata, g.options{:} );
         case 'amuse',    EEG.icaweights = amuse( tmpdata );
         case 'fobi',     [~, EEG.icaweights] = fobi( tmpdata, g.options{:} );
         case 'evd',      EEG.icaweights = evd( tmpdata, g.options{:} );
         case 'sons',     EEG.icaweights = sons( tmpdata, g.options{:} );
         case 'jadeop',   EEG.icaweights = jadeop( tmpdata, g.options{:} );
         case 'jade_td_p',EEG.icaweights = jade_td_p( tmpdata, g.options{:} );
         case 'evd24',    EEG.icaweights = evd24( tmpdata, g.options{:} );
         case 'sobi',     EEG.icawinv    = sobi( tmpdata, g.options{:} );
         case 'ng_ol',    [~, EEG.icaweights] = ng_ol( tmpdata, g.options{:} );
         case 'acsobiro', EEG.icawinv   = acsobiro( tmpdata, g.options{:} );
         case 'acrsobibpf', EEG.icawinv = acrsobibpf( tmpdata, g.options{:} );
        end
        clear tmp;
        close(fig);
     otherwise, error('Pop_runica: unrecognized algorithm');
end

% update weight and inverse matrices etc...
% -----------------------------------------
if ~isempty(fig), try close(fig); catch, end; end
if isempty(EEG.icaweights)
    EEG.icaweights = pinv(EEG.icawinv);
end
if isempty(EEG.icasphere)
    EEG.icasphere  = eye(size(EEG.icaweights,2));
end
if isempty(EEG.icawinv)
    EEG.icawinv    = pinv(EEG.icaweights*EEG.icasphere); % a priori same result as inv
end

% Reorder components by variance
% ------------------------------
meanvar = sum(EEG.icawinv.^2).*sum(transpose((EEG.icaweights *  EEG.icasphere)*EEG.data(EEG.icachansind,:)).^2)/((length(EEG.icachansind)*EEG.pnts)-1);
[~, windex] = sort(meanvar);
windex = windex(end:-1:1); % order large to small
meanvar = meanvar(windex);
EEG.icaweights = EEG.icaweights(windex,:);
EEG.icawinv    = pinv( EEG.icaweights *  EEG.icasphere );
if ~isempty(EEG.icaact)
    EEG.icaact = EEG.icaact(windex,:,:);
end

% copy back data to datasets if necessary
% ---------------------------------------
if length(g.dataset) > 1
    for i = g.dataset
        ALLEEG(i).icaweights = EEG.icaweights;
        ALLEEG(i).icasphere  = EEG.icasphere;
        ALLEEG(i).icawinv    = EEG.icawinv;
        ALLEEG(i).icachansind = g.chanind;
    end          
    ALLEEG = eeg_checkset(ALLEEG);
else
    EEG = eeg_checkset(EEG);
    ALLEEG = eeg_store(ALLEEG, EEG, g.dataset);
end

if nargin < 2 || selectamica
    if ~isempty(g.options)
        com = sprintf('EEG = pop_runica(EEG, ''icatype'', ''%s'', %s);', g.icatype, vararg2str(g.options) ); %vararg2str({ 'icatype' g.icatype 'dataset' g.dataset 'options' g.options }) );
    else
        com = sprintf('EEG = pop_runica(EEG, ''icatype'', ''%s'');',g.icatype );
    end
end

return;

function tmprank2 = getrank(tmpdata, pca_opt)
    
    tmprank = rank(tmpdata);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Here: alternate computation of the rank by Sven Hoffman
    %tmprank = rank(tmpdata(:,1:min(3000, size(tmpdata,2)))); old code
    covarianceMatrix = cov(tmpdata', 1);
    [~, D] = eig (covarianceMatrix);
    rankTolerance = 1e-7;
    tmprank2=sum (diag (D) > rankTolerance);
    if tmprank ~= tmprank2
        if nargin >= 2 && pca_opt ~= 0
            fprintf('Warning: fixing rank computation inconsistency (%d vs %d) most likely because running under Linux 64-bit Matlab\n', tmprank, tmprank2);
        end
        tmprank2 = max(tmprank, tmprank2);
    end
            
function EEG = update_datafield(EEG)
    if ~isfield(EEG, 'datfile'), EEG.datfile = ''; end
    if ~isempty(EEG.datfile)
        EEG.data = EEG.datfile;
    else 
        EEG.data = 'in set file';
    end
    EEG.icaact = [];            
