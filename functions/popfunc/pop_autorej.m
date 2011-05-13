% pop_autorej() - perform automatic artifact epoch detection and rejection 
%
% Usage:
%   >> [EEG, rmepochs] = pop_autorej( EEG, 'key', 'val');
%
% Inputs:
%   EEG       - input EEG structure
%
% Optional inputs:
%   'threshold'  - [float] Threshold limit for detection of extremely large 
%                  fluctuations (uV) {default: 1000}
%   'electrodes' - [integer] Use these channel indices for detection
%                  of improbable data {default: all channels}
%   'icacomps'   - [integer] Use these component activities (instead of 
%                  channel data) for detection of improbable data.
%                  {default: none}
%   'startprob'  - [float] Probability threshold (in std. dev.) for
%                  detection of improbable data {default: 5}
%   'maxrej'     - [float] Maximum % of total trials to reject per
%                  iteration {default: 5}
%   'eegplot'    - ['on'|'off'] Allow for visual inspection of rejected
%                  trials {default: 'off' if called from the command line, 
%                  'on' if called from the gui}
%   'nogui'      - ['on'|'off'] Do not pop up a gui window to ask for 
%                  input parameters {default: 'off'}
%
% Outputs:
%   EEG          - EEGLAB data structure
%   rmepochs     - [integer] rejected trial indices
%
% Function description: 
%     pop_autorej() first detects extremely large potential fluctuations;
%     this is mostly to detect artifacts from electrical surges or other 
%     unreasonably large amplitude events. Then it applies the following.
%     The function rejects data epochs containing data values outside a 
%     given standard deviation (s.d.) threshold entered by the user (e.g., 
%     3 s.d.'s). In each iteration, if the number of epochs that are
%     thus marked for rejection are fewer than 'maxrej' (by default, 5%), it  
%     then rejects the beyond-threshold data epochs and iterates. If the number 
%     of epochs marked  for rejection is more than 5% of the total number 
%     of data epochs, it does not reject them, but instead  increases the 
%     s.d. threshold by 0.5 s.d. and iterates. When no more data epochs are 
%     found to exceed the current s.d. threshold, it lowers the threshold 
%     by 0.5 s.d. and continues to iterate until either no more epochs are 
%     rejected or until 8 iterations have been performed.
%
% Authors: Julie Onton and Arnaud Delorme, SCCN/INC/UCSD, 2007-

% Copyright (C) Arnaud Delorme, SCCN/INC/UCSD, 2007
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

function [EEG, rmep, com ] = pop_autorej(EEG, varargin);
    
    DEFAULT_STARTPROB = 5;    % std devs
    DEFAULT_MAXREJ    = 5;    % std devs
    DEFAULT_THRESH    = 1000; % uV
    if nargin < 1
        help pop_autorej;
        return;
    end;
    rmep = [];
    com = '';
    
    opt = finputcheck(varargin, { 'startprob'    'real'    []     DEFAULT_STARTPROB; ...
                                  'electrodes'   'real'    []     [1:EEG.nbchan]; ...
                                  'icacomps'     'real'    []     []; ...
                                  'maxrej'       'real'    []     DEFAULT_MAXREJ; ...
                                  'eegplot'      'string'  { 'on';'off' }    'off'; ...
                                  'nogui'        'string'  { 'on';'off' }    'off'; ...
                                  'threshold'    'real'    []     DEFAULT_THRESH }, 'pop_autorej');
    if isstr(opt), error(opt); end;
    
    % pop-up GUI for rejecting artifacts
    % ----------------------------------
    if strcmpi(opt.nogui, 'off')
        
        promptstr    = { { 'style' 'text'       'string' 'Detection of extremelly large fluctuations (channels only)' 'fontweight' 'bold' } ...
                         { 'style' 'text'       'string' 'Threshold limit (microV)' } ...
                         { 'style' 'edit'       'string' '1000' }  ...
                         { 'style' 'text'       'string' ' ' }  ...
                         { 'style' 'text'       'string' 'Detection of unprobable activity (channels or ICA)' 'fontweight' 'bold' } ...
                         { 'style' 'text'       'string' 'Do not use these channel indices (default=all)' } ...
                         { 'style' 'edit'       'string' '' }  ...
                         { 'style' 'text'       'string' 'Use these ICA components instead of data channels' } ...
                         { 'style' 'edit'       'string' '' }  ...
                         { 'style' 'text'       'string' 'Probability threshold (std. dev.)' } ...
                         { 'style' 'edit'       'string' '5' } ...
                         { 'style' 'text'       'string' 'Maximum % of total trials to reject per iteration' } ...
                         { 'style' 'edit'       'string' '5' }  ...
                         { 'style' 'text'       'string' ' ' }  ...
                         { 'style' 'text'       'string' 'Check box for visual inspection of results' } ...
                         { 'style' 'checkbox'   'string' '' 'value' 1 } };
        geometry = { [1] [2 1] [1] [1] [2 1] [2 1] [2 1] [2 1] [1] [2 1]};
        result       = inputgui( 'geometry', geometry, 'uilist', promptstr, ...
                                 'helpcom', 'pophelp(''pop_autorej'')', ...
                                 'title', 'Automatic artifact rejection -- pop_autorej()');
        if isempty(result), return; end;
        
        options = { 'nogui' 'on' };
        if ~strcmpi(result{1}, '1000'), options = { options{:} 'threshold' str2num(result{1}) }; end;
        if ~isempty(result{2}),         options = { options{:} 'electrodes' setdiff([1:EEG.nbchan], str2num(result{2})) }; end;
        if ~isempty(result{3}),         options = { options{:} 'icacomps'   str2num(result{3}) }; end;
        if ~strcmpi(result{4}, '5'),    options = { options{:} 'startprob' str2num(result{4}) }; end;
        if ~strcmpi(result{5}, '5'),    options = { options{:} 'maxrej'    str2num(result{5}) }; end;
        if result{6}, options = { options{:} 'eegplot' 'on' }; end;
    
        [ EEG rmep com ] = pop_autorej(EEG, options{:});
        return;
    end;
  
    EEGIN = EEG; % backup EEG structure
    if ~isempty(opt.icacomps), 
        processdat = 0;
        complist   = opt.icacomps;
    else
        processdat = 1;
        complist   = opt.electrodes;
    end;
    
    % rejection of extremelly large fluctuations
    % ------------------------------------------
    fprintf('\nRunning auto-rejection protocol...\n');
    rmep = zeros(1,0);
    alleps = [1:EEG.trials];
    EEG = pop_eegthresh(EEG,1,[1:size(EEG.data,1)],-opt.threshold,opt.threshold,EEG.xmin,EEG.xmax,0,0);
    numrej = length(find(EEG.reject.rejthresh));  % count number of epochs marked
%     if numrej > 0
%         rmep(1,end+1:end+length(find(EEG.reject.rejthresh))) = alleps(find(EEG.reject.rejthresh));
%         alleps(find(EEG.reject.rejthresh)) = [];
%         EEG = pop_rejepoch( EEG,EEG.reject.rejthresh,0); % actually reject high prob epochs
%         fprintf('\nRe-baselining after large amplitude artifact removed (does not affect the data)...\n');
%         EEG = pop_rmbase( EEG, [EEG.xmin*1000 EEG.xmax*1000]);
%     end;
    
    %--------------------------------------------

    EEG = pop_jointprob(EEG, processdat, complist ,opt.startprob,opt.startprob,0,0);% calculate component probabilities
    if processdat % if rejection based on channels
        numrej = length(find(EEG.reject.rejjp));  % count number of epochs marked
        
    else% if rejection based on ICs
        numrej = length(find(EEG.reject.icarejjp));  % count number of epochs marked
    end;
    if (numrej/EEG.trials) < opt.maxrej/100
        if processdat
            rmep(1,end+1:end+length(find(EEG.reject.rejjp))) = alleps(EEG.reject.rejjp);
            alleps(EEG.reject.rejjp) = [];
            EEG = pop_rejepoch( EEG,EEG.reject.rejjp,0); % actually reject high prob epochs
        else
            rmep(1,end+1:end+length(find(EEG.reject.icarejjp))) = alleps(EEG.reject.icarejjp);
            alleps(EEG.reject.icarejjp) = [];
            EEG = pop_rejepoch( EEG,EEG.reject.icarejjp,0); % actually reject high prob epochs
        end;
    else
        fprintf('Re-adjusting probability limits and running again...*********\n');
        opt.startprob = opt.startprob + .5;
    end;
    repeat = 1; maxiter = 0;
    while repeat == 1 % keep running probability until there are no epochs above threshold
        if numrej > 0
            EEG = pop_jointprob(EEG,processdat,complist ,opt.startprob,opt.startprob,0,0);
            if processdat % if rejection based on channels
                numrej = length(find(EEG.reject.rejjp));  % count number of epochs marked
                
            else% if rejection based on ICs
                numrej = length(find(EEG.reject.icarejjp));  % count number of epochs marked
            end;
            if (numrej/EEG.trials) < opt.maxrej/100
                if processdat
                    rmep(1,end+1:end+length(find(EEG.reject.rejjp))) = alleps(EEG.reject.rejjp);
                    alleps(EEG.reject.rejjp) = [];
                    EEG = pop_rejepoch( EEG,EEG.reject.rejjp,0); % actually reject high prob epochs
                else
                    rmep(1,end+1:end+length(find(EEG.reject.icarejjp))) = alleps(EEG.reject.icarejjp);
                    alleps(EEG.reject.icarejjp) = [];
                    EEG = pop_rejepoch( EEG,EEG.reject.icarejjp,0);
                end;
            else
                opt.startprob = opt.startprob + 0.5; EEG.reject.icarejjp = [];EEG.reject.rejjpE = [];
                fprintf('Re-adjusting probability limits and running again...*********\n');
            end;
        else
            if opt.startprob > 5 & maxiter < 8 % don't decrease and startover more than 8 times
                fprintf('Decreasing probability limits for final pruning...######\n');
                opt.startprob = opt.startprob - 0.5; numrej = 1; maxiter = maxiter+1; % repeat process back to 5 stds
            else
                if maxiter > 8
                    opt.maxrej = 15; % go through last round with a high threshold
                else
                    repeat = 0;
                end;
            end;
        end;
    end;
    
    % run kurtosis check
    % ------------------
    disp('Final kurotsis reject...');
    EEG = pop_rejkurt(EEG,processdat,complist ,6,6,0,0);
    numrej = length(find(EEG.reject.icarejkurt));  % count number of epochs marked
    if numrej > 0
        rmep(1,end+1:end+length(find(EEG.reject.icarejkurt))) = alleps(EEG.reject.icarejkurt);
        alleps(EEG.reject.icarejjp) = [];
        EEG = pop_rejepoch( EEG,EEG.reject.icarejkurt,0);
    end;
    
    %--------------------------------------------
    
    EEG = EEGIN;
    if strcmpi(opt.eegplot, 'on')
        EEG = EEGIN;
        EEG.reject.rejauto  = zeros(1, length(EEG.trials));
        EEG.reject.rejauto(rmep) = 1;
        EEG.reject.rejautoE = zeros(EEG.nbchan, EEG.trials);
        
       	colrej = EEG.reject.rejmanualcol;
        rej    = EEG.reject.rejauto;
        rejE   = EEG.reject.rejautoE;
        elecrange = complist;
        superpose = 0;
        icacomp   = processdat;
        macrorej  = 'EEG.reject.rejauto';
        macrorejE = 'EEG.reject.rejautoE';     
        reject = 1;
        
        eeg_rejmacro; % script macro for generating command and old rejection arrays
        
        eegplot( EEG.data(elecrange,:,:), 'srate', EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}); 
    else
        EEG = pop_rejepoch( EEG, rmep, 0); % actually reject high prob epochs        
    end;
    
    com = sprintf('EEG = pop_autorej(EEG, %s);', vararg2str( varargin ));
    
