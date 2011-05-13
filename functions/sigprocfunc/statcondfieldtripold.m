% statcondfiledtrip()  - same as statcond except that it uses the fieldtrip
%                        statistical functions. This is useful to perform
%                        a wider variety of corrections for multiple 
%                        comparisons for instance.
% Usage:
%          >> [stats, df, pvals, surrog] = statcond( data, 'key','val'... );
% Inputs:
%   data       = same as for statcond()
%
% Optional inputs:
%   'paired'   = ['on'|'off'] pair the data array {default: 'on' unless 
%                the last dimension of data array is of different lengths}.
%   'mode'     = ['perm'|'param'] mode for computing the p-values:
%                 'param' = parametric testing (standard ANOVA or t-test); 
%                 'perm' = non-parametric testing using surrogate data 
%                  made by permuting the input data {default: 'perm'}
%   'naccu'    = this input is passed on as 'numrandomization' to Fieldtrip
%   'chanlocs' = EEGLAB channel location structure to perfom statistics
%                and cluster correction for multiple comparisons accross 
%                channels (native fieldtrip entries may also be used).
%
% Fieldtrip options:
%   Any option to the freqanalysis, the statistics_montecarlo, the
%   statistics_analysis, statistics_stat, statistics_glm may be used
%   using 'key', val argument pairs.
%
% Outputs:
%   stats      = F- or T-value array of the same size as input data without 
%                the last dimension. A T value is returned only when the data 
%                includes exactly two conditions.
%   df         = degrees of freedom, a (2,1) vector, when F-values are returned
%   pvals      = array of p-values. Same size as input data without the last
%                data dimension. All returned p-values are two-tailed.
%   surrog     = surrogate data array (same size as input data with the last 
%                dimension filled with a number ('naccu') of surrogate data sets.
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, 2005-
%         With thanks to Robert Oostenveld for fruitful discussions 
%         and advice on this function.
%
% See also: freqanalysis(), statistics_montecarlol()

% Copyright (C) Arnaud Delorme
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


function [ ori_vals, df, pvals, surrogval ] = statcondfieldtrip( data, varargin );
    
    if nargin < 1
        help statcondfieldtrip;
        return;
    end;
    
    [g cfgparams] = finputcheck( varargin, { 'naccu'    'integer'   [1 Inf]             200;
                                             'mode'     'string'    { 'param','perm','bootstrap' }  'param';
                                             'chanlocs' 'struct'    { }   struct([]);
                                             'method'   'string'    { 'montecarlo','analytic','stat','glm' } 'analytic';
                                             'paired'   'string'    { 'on','off' }      'on' }, 'statcond', 'ignore');
    if isstr(g), error(g); end;    
    
    if size(data,2) == 1, data  = transpose(data); end; % cell array transpose
    
    tmpsize   = size(data{1});
    
    % cfg configuration for Fieldtrip
    % -------------------------------
    cfg = struct(cfgparams{:});
    cfg.method      = g.method;
    if ~strcmpi(g.mode, 'param'), cfg.method = 'montecarlo'; end;
    if ~isempty(g.chanlocs)
        EEG = eeg_emptyset;
        EEG.chanlocs = g.chanlocs;
        EEG.nbchan   = length(g.chanlocs);
        EEG.data     = zeros(EEG.nbchan,100,1);
        EEG = eeg_checkset(EEG);
        tmpcfg = eeglab2fieldtrip(EEG, 'preprocessing', 'none');
        lay = prepare_layout(tmpcfg, tmpcfg);
        tmpcfg.layout        = lay;
        tmpcfg.neighbourdist = 30; %mean(abs([ EEG.chanlocs.X ]));
        cfg.neighbours    = neighbourselection(tmpcfg, []);
    else
        cfg.neighbours  = {};
    end;
    cfg.feedback    = 'no';
    cfg.ivar        = 1;
    cfg.numrandomization = g.naccu;
    
    % test if data can be paired
    % --------------------------
    if length(unique(cellfun('size', data, ndims(data{1}) ))) > 1
        g.paired = 'off'; 
    end;
    fprintf('%d x %d, ', size(data,1), size(data,2));
    if strcmpi(g.paired, 'on')
         fprintf('paired data, ');
    else fprintf('unpaired data, ');
    end;
    if size(data,1) == 1 & size(data,2) == 2
         fprintf('computing T values\n');
    else fprintf('computing F values\n');
    end;
    
    % output text
    % -----------
    if ~strcmpi(g.mode, 'param')
        fprintf('Accumulating (of %d):', g.naccu);
    end;

    if size(data,1) == 1, % only one row
        
        if size(data,2) == 2 & strcmpi(g.paired, 'on')
            
            % paired t-test (very fast)
            % -------------
            cfg.statistic   = 'depsamplesT';
            [newdata design1 design2 design3] = makefieldtripdata(data, 0);
            if ~isempty(g.chanlocs)
                for index = 1:length(newdata)
                    newdata{index}.powspctrm = squeeze(newdata{index}.powspctrm);
                    newdata{index}.label     = { g.chanlocs.labels };
                    newdata{index}.freq      = 1;
                end;
            end;
            cfg.design      = [ design1; design3 ];
            cfg.uvar        = 2;
            cfg.alpha       = 0.05;

            cfg
            stat            = freqstatistics(cfg, newdata{:});
            ori_vals        = stat.stat;
            df              = [];
            pvals           = stat.prob;
            return;
            
        elseif size(data,2) == 2 & strcmpi(g.paired, 'off')
            
            % paired t-test (very fast)
            % -------------
            cfg.statistic   = 'indepsamplesT';
            [newdata design1] = makefieldtripdata(data, 0);
            cfg.design      = design1;            
            stat            = freqstatistics(cfg, newdata{:});
            ori_vals        = stat.stat;
            df              = [];
            pvals           = stat.prob;
            return;
            
        elseif strcmpi(g.paired, 'on')
            
            % one-way ANOVA (paired) this is equivalent to unpaired t-test
            % -------------
            cfg.statistic   = 'depsamplesF';
            [newdata design1 design2 design3] = makefieldtripdata(data, 0);
            cfg.design      = [ design1; design3 ];
            cfg.uvar        = 2;
            stat            = freqstatistics(cfg, newdata{:});
            ori_vals        = stat.stat;
            df              = [];
            pvals           = stat.prob;
            return;
            
        else
            % one-way ANOVA (unpaired) 
            % -------------
            cfg.statistic   = 'indepsamplesF';
            [newdata design1] = makefieldtripdata(data, 0);
            cfg.design      = [ design1 ];
            stat            = freqstatistics(cfg, newdata{:});
            ori_vals        = stat.stat;
            df              = [];
            pvals           = stat.prob;
            return;
            
        end;
    else
        if strcmpi(g.paired, 'on')
            
            % two-way ANOVA (paired) 
            % -------------
            cfg.statistic   = 'anovan';
            [newdata design1 design2 design3] = makefieldtripdata(data, 0);
            cfg.design      = [ design1; design2; design3 ];
            cfg.effect      = 'X1*X2';
            cfg.ivar        = [1 2];
            cfg.uvar        = 3;
            stat            = freqstatistics(cfg, newdata{:});
            ori_vals        = stat.stat;
            df              = [];
            pvals           = stat.prob;
            return;
            
        else
            
            % two-way ANOVA (unpaired)
            % -------------
            cfg.statistic   = 'anovan';
            [newdata design1 design2] = makefieldtripdata(data, 0);
            cfg.design      = [ design1; design2 ];
            cfg.effect      = 'X1*X2';
            cfg.ivar        = [1 2];
            stat            = freqstatistics(cfg, newdata{:});
            ori_vals        = stat.stat;
            df              = [];
            pvals           = stat.prob;
            return;
            
        end;
    end;
                
function val = myndims(a)
    if ndims(a) > 2
        val = ndims(a);
    else
        if size(a,1) == 1,
            val = 2;
        elseif size(a,2) == 1,
            val = 1;
        else
            val = 2;
        end;
    end; 
  
function [newdata, design1, design2, design3] = makefieldtripdata(data, chanflag);
    
    newdata = {};
    for i = 1:length(data(:))
    
        switch myndims(data{1})
          case 1, 
            newdata{i}.powspctrm = data{i};
            
          case 2,
            if chanflag
                 newdata{i}.powspctrm = transpose(data{i});
            else newdata{i}.powspctrm = reshape(transpose(data{i}), size(data{i},2), 1, size(data{i},1));
            end;
            
          case 3,
            if chanflag
                 newdata{i}.powspctrm = permute(data{i}, [3 1 2]);
            else newdata{i}.powspctrm = permute(data{i}, [3 4 1 2]);
            end;
            
          case 4,
            if chanflag
                 newdata{i}.powspctrm = permute(data{i}, [4 1 2 3]);
            else newdata{i}.powspctrm = permute(data{i}, [3 5 1 2]);
            end;

        end;
        
        newdata{i}.dimord    = 'rpt_chan_freq_time';
        newdata{i}.label     = cell(1,size(newdata{i}.powspctrm,2));
        newdata{i}.label(:)  = { 'cz' };
        newdata{i}.freq      = [1:size(newdata{i}.powspctrm,3)];
        newdata{i}.time      = [1:size(newdata{i}.powspctrm,4)];

    end;
    
    design1 = [];
    design2 = [];
    design3 = [];
    for i = 1:size(data,2)
        for j = 1:size(data,1)
            ij = j+(i-1)*size(data,1);
            design1 = [ design1 ones(1, size(newdata{i}.powspctrm,1))*i ];
            design2 = [ design2 ones(1, size(newdata{i}.powspctrm,1))*j ];
            design3 = [ design3 [1:size(newdata{i}.powspctrm,1)] ];
        end;
    end;
