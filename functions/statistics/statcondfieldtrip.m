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
%   'method'   = ['permutation'|'parametric'] method for computing the p-values:
%                'parametric' = parametric testing (standard ANOVA or t-test); 
%                'permutation' = non-parametric testing using surrogate data 
%                made by permuting the input data. Note that if 'bootstrap'
%                is given as input, it is interpreted as 'permutation' 
%                Default is 'parametric'. Note that 'parametric'
%                corresponds to the 'analytic' method of Fieldtrip and
%                'permutation' correspond to the 'montecarlo' method.
%   'naccu'    = this input is passed on as 'numrandomization' to Fieldtrip
%   'neighbours' = Fieldtrip channel neighbour structure to perform statistics
%                and cluster correction for multiple comparisons across 
%                channels.
%   'alpha'    = [float] p-value threshold value. Allow returning
%                confidence intervals and mask (requires structoutput below).
%   'structoutput' = ['on'|'off'] return an output structure instead of 
%                the regular output. Allow to output mask and confidence
%                intervals.
%
% Fieldtrip options:
%   Any option to the freqanalysis, the statistics_montecarlo, the
%   statistics_analysis, statistics_stat, statistics_glm may be used
%   using 'key', val argument pairs. Note that although 'fieldtripmcorrect'
%   is used by std_stat, this function uses 'mcorrect'. See Fieldtrip
%   documentation for more information.
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


function [ ori_vals, df, pvals ] = statcondfieldtrip( data, varargin )
    
    if nargin < 1
        help statcondfieldtrip;
        return;
    end
    
    [g, cfgparams] = finputcheck( varargin, { 'naccu'      ''          []             [];
                                             'method'     'string'    { }            'param';
                                             'mode'       'string'    { }            ''; % deprecated (old method)
                                             'chanlocs'   'struct'    { }            struct([]);
                                             'chandim'    'integer'   []             0;
                                             'alpha'      'real'      []             NaN;
                                             'neighbours' 'struct'    { }            struct([]);
                                             'structoutput' 'string'  { 'on','off' }      'off';
%                                             'method'    'string'    {  } 'analytic'; % 'montecarlo','analytic','stat','glm'
                                             'paired'     'string'    { 'on','off' }      'on' }, ... % NOT USED IN THIS FUNCTION
                                             'statcond', 'ignore');
    if ischar(g), error(g); end;    
    if ~isempty(g.mode), g.method = g.mode; end
    if strcmpi(g.method, 'parametric'), g.method = 'param'; end
    if strcmpi(g.method, 'permutation'), g.method = 'montecarlo'; end
    if isfield(g, 'fieldtripmcorrect'), error('THIS FUNCTION DOES NOT RECOGNIZE "fieldtripmcorrect" USE "mcorrect" INSTEAD'); end
    if ~isempty(g.neighbours) && isempty(g.chanlocs)
        g.chanlocs = struct('labels', { g.neighbours(:).label });
    end
    if size(data,2) == 1, data  = transpose(data); end; % cell array transpose
    alphaset = fastif(isnan(g.alpha) || isempty(g.alpha), 0, 1);
    
    % remove first dimension for all input if necessary
    % necessary for scalp topographies which are given as 1 x nelec x subj
    % -------------------------------------------------
    ndim = size(data{1});
    if size(data{1},1) == 1
        for index = 1:length(data(:))
            data{index} = squeeze(data{index});
        end
    end
    tmpsize = size(data{1});
    
    % find the channel dimension if any
    % ---------------------------------
    if ~isempty(g.neighbours) && g.chandim == 0
        for index = 1:ndims(data{1})
            if size(data{1},index) == length(g.neighbours);
                if g.chandim == 0
                    g.chandim = index;
                else
                    error('Multiple possibilities for the channel dimension, please specify manually');
                end
            end
        end
    end
    
    % cfg configuration for Fieldtrip
    % -------------------------------
    cfg = struct(cfgparams{:});
    cfg.method = g.method;
    if strcmpi(g.method, 'param') || strcmpi(g.method, 'parametric')
         cfg.method = 'analytic';
    elseif strcmpi(g.method, 'perm') && strcmpi(g.method, 'permutation') || strcmpi(g.method, 'bootstrap')
         cfg.method = 'montecarlo';
    end
    if ~isempty(g.neighbours)
        cfg.neighbours = g.neighbours;
    end
    if isfield(cfg, 'mcorrect')
         if strcmpi(cfg.mcorrect, 'none')
             cfg.mcorrect = 'no';
         end
         cfg.correctm = cfg.mcorrect;
    else cfg.mcorrect = [];
    end
    cfg.feedback    = 'no';
    cfg.ivar        = 1;
    cfg.alpha       = fastif(alphaset, g.alpha, 0.05);
    cfg.numrandomization = g.naccu; 
    
    % test if data can be paired
    % --------------------------
    if length(unique(cellfun('size', data, ndims(data{1}) ))) > 1
        g.paired = 'off'; 
    end
    fprintf('%d x %d, ', size(data,1), size(data,2));
    if strcmpi(g.paired, 'on')
         fprintf('paired data, ');
    else fprintf('unpaired data, ');
    end
    if size(data,1) == 1 && size(data,2) == 2
         fprintf('computing T values\n');
    else fprintf('computing F values\n');
    end
    
    % set randomizations
    % ------------------
    if strcmpi(cfg.method, 'montecarlo') && isempty(cfg.numrandomization)
        cfg.numrandomization = 200;
        if ~strcmpi(cfg.mcorrect, 'no'), cfg.numrandomization = cfg.numrandomization*20; end
    end
    cfg.correcttail = 'alpha';
    
    if size(data,1) == 1 % only one row
        
        if size(data,2) == 2 && strcmpi(g.paired, 'on')
            
            % paired t-test (very fast)
            % -------------
            cfg.statistic   = 'depsamplesT';
            [newdata design1 design2 design3] = makefieldtripdata(data, g.chandim, g.chanlocs);
            cfg.design      = [ design1; design3 ];
            cfg.uvar        = 2;
            stat            = ft_freqstatistics(cfg, newdata{:});
            if isfield(stat, 'df')
                 df = stat.df;
            else df = [];
            end
            
        elseif size(data,2) == 2 && strcmpi(g.paired, 'off')
            
            % paired t-test (very fast)
            % -------------
            cfg.statistic   = 'indepsamplesT';
            [newdata design1] = makefieldtripdata(data, g.chandim, g.chanlocs);
            cfg.design      = design1;            
            stat            = ft_freqstatistics(cfg, newdata{:});
            if isfield(stat, 'df')
                 df = stat.df;
            else df = [];
            end
            
        elseif strcmpi(g.paired, 'on')
            
            % one-way ANOVA (paired) this is equivalent to unpaired t-test
            % -------------
            cfg.tail        = 1;
            cfg.correcttail = 'no';
            tmpP = fileparts(which('ft_freqstatistics'));
            if exist(fullfile(tmpP, 'statfun', 'ft_statfun_depsamplesFmultivariate.m'))
                 cfg.statistic   = 'depsamplesFunivariate';
            else cfg.statistic   = 'depsamplesF';
            end
            [newdata design1 design2 design3] = makefieldtripdata(data, g.chandim, g.chanlocs);
            cfg.design      = [ design1; design3 ];
            cfg.uvar        = 2;
            stat            = ft_freqstatistics(cfg, newdata{:});
            if isfield(stat, 'dfnum')
                 df = [stat.dfnum stat.dfdenom];
            else df = [];
            end
            
        else
            % one-way ANOVA (unpaired) 
            % -------------
            cfg.tail        = 1;
            cfg.correcttail = 'no';
            cfg.statistic   = 'indepsamplesF';
            [newdata design1] = makefieldtripdata(data, g.chandim, g.chanlocs);
            cfg.design      = [ design1 ];
            warning off;
            stat            = ft_freqstatistics(cfg, newdata{:});
            warning on;
            if isfield(stat, 'dfnum')
                 df = [stat.dfnum stat.dfdenom];
            else df = [];
            end
            
        end
        
    else
        if strcmpi(g.paired, 'on')
            
            % two-way ANOVA (paired) 
            % -------------
            cfg.tail        = 1;
            cfg.correcttail = 'no';
            cfg.statistic   = 'depsamplesFmultivariate';
            [newdata design1 design2 design3] = makefieldtripdata(data, g.chandim, g.chanlocs);
            cfg.design      = [ design1; design2; design3 ];
            cfg.ivar        = [1 2];
            cfg.uvar        = 3;
            stat            = ft_freqstatistics(cfg, newdata{:});
            ori_vals        = stat.stat;
            if isfield(stat, 'df')
                 df = stat.df;
            else df = [];
            end
            
        else
            
            % two-way ANOVA (unpaired)
            % -------------
            cfg.tail        = 1;
            cfg.correcttail = 'no';
            cfg.statistic   = 'indepsamplesF';
            cfg.clustercritval = 4.5416; % 95 percentile of n =10000; a = { rand(n,10) rand(n,10); rand(n,10) rand(n,10) }; [F df p ] = statcondfieldtrip(a, 'paired', 'off');
            [newdata design1 design2] = makefieldtripdata(data, g.chandim, g.chanlocs);
            if ~isempty(g.chanlocs)
                for index = 1:length(newdata)
                    newdata{index}.powspctrm = squeeze(newdata{index}.powspctrm);
                    newdata{index}.label     = { g.chanlocs.labels };
                    newdata{index}.freq      = 1;
                end
            end
            cfg
            newdata{1}
            cfg.design      = [ design1; design2 ];
            cfg.effect      = 'X1*X2';
            cfg.ivar        = [1 2];
            stat            = ft_freqstatistics(cfg, newdata{:});
            ori_vals        = stat.stat;
            df              = stat.df;
            
        end
    end
    
    ori_vals  = stat.stat;
    pvals     = stat.prob;
    
    if size(ori_vals,1) ~= size(data{1},1) && size(ori_vals,1) == 1
        ori_vals = reshape(ori_vals, size(ori_vals,2), size(ori_vals,3), size(ori_vals,4));
        pvals    = reshape(pvals   , size(pvals   ,2), size(pvals   ,3), size(pvals   ,4));
        if isfield(stat, 'mask')
            stat.mask = reshape(stat.mask   , size(stat.mask   ,2), size(stat.mask   ,3), size(stat.mask   ,4));
        end
    end
    
    if strcmpi(g.structoutput, 'on')
        outputstruct.mask = stat.mask;
        outputstruct.pval = pvals;
        if length(data(:)) == 2
             outputstruct.t = ori_vals;
        else outputstruct.f = ori_vals;
        end
        outputstruct.stat   = ori_vals;
%         outputstruct.method = g.method;
%         outputstruct.pval   = pvals;
%         outputstruct.df     = df;
%         outputstruct.surrog = surrogval;
%         if length(data(:)) == 2
%              outputstruct.t = ori_vals;
%         else outputstruct.f = ori_vals;
%         end
        ori_vals = outputstruct;
    end
    
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
        end
    end; 
  
function [newdata, design1, design2, design3] = makefieldtripdata(data, chandim, chanlocs);
    
    newdata = {};
    swapdim = [];
    for i = 1:length(data(:))
    
        newdata{i}.dimord    = 'rpt_chan_freq_time';
        newdata{i}.powspctrmdimord = 'rpt_chan_freq_time';
        switch myndims(data{1})
          case 1, 
            newdata{i}.powspctrm = data{i};
            
          case 2,
            if chandim
                 newdata{i}.powspctrm = transpose(data{i});
            else newdata{i}.powspctrm = reshape(transpose(data{i}), size(data{i},2), 1, size(data{i},1));
            end
            
          case 3,
            if chandim == 2 % chandim can be 1 or 2
                swapdim = [2 1];
            end
            if chandim
                 newdata{i}.powspctrm = permute(data{i}, [3 1 2]);
            else newdata{i}.powspctrm = permute(data{i}, [3 4 1 2]); % 4 is a singleton dimension
            end
            
          case 4,
            newdata{i}.powspctrm = permute(data{i}, [4 3 1 2 ]); % Fixed dimension from [4 1 2 3]
        end
        
        newdata{i}.label     = cell(1,size(newdata{i}.powspctrm,2));
        newdata{i}.label(:)  = { 'cz' };
        for ic = 1:length(newdata{i}.label)
            newdata{i}.label{ic} = [ 'c' num2str(ic) ];
        end
        newdata{i}.freq      = [1:size(newdata{i}.powspctrm,3)];
        newdata{i}.time      = [1:size(newdata{i}.powspctrm,4)];
        
        % below in case channels are specified
        % not that statistics are done on time x frequencies or channels
        % so time x frequency x channels do not work yet here
        if ~isempty(chanlocs)
            newdata{i}.powspctrm = squeeze(newdata{i}.powspctrm);
            newdata{i}.label     = { chanlocs.labels };
            newdata{i}.freq      = 1;
            newdata{i}.time      = 1;
        end
        if isempty(chanlocs) && size(newdata{i}.powspctrm,2) ~= 1
            newdata{i}.dimord    = 'rpt_freq_time';
            newdata{i}.powspctrmdimord = 'rpt_freq_time';
        end
    end
    
    design1 = [];
    design2 = [];
    design3 = [];
    for i = 1:size(data,2)
        for j = 1:size(data,1)
            nrepeat = size(data{i}, ndims(data{i}));
            ij = j+(i-1)*size(data,1);
            design1 = [ design1 ones(1, nrepeat)*i ];
            design2 = [ design2 ones(1, nrepeat)*j ];
            design3 = [ design3 [1:nrepeat] ];
        end
    end
        
