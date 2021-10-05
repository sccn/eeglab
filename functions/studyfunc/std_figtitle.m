% std_figtitle() - Generate plotting figure titles in a cell array
%
% Usage:    
%   >> celltxt = std_figtitle(key1, val1, key2, val2);  
%
% Inputs:
%  'subject'     - [string] Subject name
%  'datatype'    - [string] data type (For example 'erp')
%  'chanlabels'  - [cell array or string] channel names
%  'compnames'   - [cell array or string] component names
%  'factor1'     - [cell array] name of factor1
%  'factor2'     - [cell array] name of factor2
%  'factor1vals' - [cell array] values of factor1
%  'factor2vals' - [cell array] values of factor1
%  'vals'        - [cell array or real array] value for plot panel (for 
%                  example, 0 ms)
%  'valsunit'    - [cell array or string] unit value (i.e. ms)
%
% Optional statistical titles:
%  'effect'        - ['main'|'marginal'] use main or marginal (default) effect.
%  'factor1stat'   - ['on'|'off'] stat enabled for factor 1. default is 'off'
%  'factor2stat'   - ['on'|'off'] default is 'off'
%  'mcorrect'      - [string] correction for multiple comparisons. Default is
%                  empty.
%  'statistics'  - [string] type of statistics
%  'threshold'   - [real] threshold value
% 
% Optional grouping:
%  'factor1grouped'  - ['on'|'off'] group factor 1 together default is 'off'
%  'factor2grouped'  - ['on'|'off'] default is 'off'
% 
% Outputs:
%   celltxt      - cell array of text string, one for each set of condition
%
% Example:
%  std_figtitle('subject', 'toto', 'condnames', { 'test1' 'test2' }, ...
%               'vals' , { [100] [4] }, 'valsunit', { 'ms' 'Hz' })
%
% Authors: Arnaud Delorme, SCCN/UCSD, Feb 2010

% Copyright (C) Arnaud Delorme, arno@salk.edu
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

function [all_titles, alllegends ] = std_figtitle(varargin)

if nargin < 1
    help std_figtitle;
    return;
end

alllegends = {};
opt = finputcheck( varargin, { 'chanlabels'  {'cell','string'}   []     {};
                               'factor1'     'string'            []     '';
                               'factor2'     'string'            []     '';
                               'factor1vals' {'cell','string'}   []     '';
                               'factor2vals' {'cell','string'}   []     '';
                               'factor1stat' 'string'            {'on','off'}     'off';
                               'factor2stat' 'string'            {'on','off'}     'off';
                               'factor1grouped' 'string'            {'on','off','together','apart'}     'off';
                               'factor2grouped' 'string'            {'on','off','together','apart'}     'off';
                               'condnames'   {'cell','string'}   []     '';
                               'cond2names'  {'cell','string'}   []     '';
                               'condstat'    'string'            {'on','off'}     'off';
                               'cond2stat'   'string'            {'on','off'}     'off';
                               'condgroup'   'string'            {'on','off','together','apart'}     'off';
                               'cond2group'  'string'            {'on','off','together','apart'}     'off';
                               'plotmode'    'string'            {'normal','condensed'}              'normal';
                               'plotsubjects' 'string'           {'on','off'}                        'off';
                               'threshold'   'real'              []     NaN;
                               'statistics'  'string'            []     '';
                               'mcorrect'    'string'            []     '';
                               'effect'      'string'            {'main','marginal'}     'marginal';
                               'datatype'    'string'            []     '';
                               'clustname'   'string'            []     '';
                               'compnames'   {'cell','string'}   []     {};
                               'vals'        {'cell','real'}     []     {}; % just for titles
                               'valsunit'    {'cell','string'}   []     {}; % just for titles
                               'subject'     'string'            []              '' }, 'std_figtitle'); %, 'ignore');
if ischar(opt), error(opt); end

% convert
if ~isempty(opt.condnames),  opt.factor1vals = opt.condnames; end
if ~isempty(opt.cond2names), opt.factor2vals = opt.cond2names; end
if ~isempty(opt.condstat),   opt.factor1stat = opt.condstat; end
if ~isempty(opt.cond2stat),  opt.factor2stat = opt.cond2stat; end
if ~isempty(opt.condgroup),  opt.factor1grouped = opt.condgroup; end
if ~isempty(opt.cond2group), opt.factor2grouped = opt.cond2group; end

if ~iscell(opt.vals),       opt.vals       = { opt.vals }; end
ncori  = length(opt.factor1vals);
nc2ori = length(opt.factor2vals);
if ~isempty(opt.vals) && ~isempty(opt.vals{1}) && ~isnan(opt.vals{1}(1)), opt.factor1grouped = 'off'; opt.factor2grouped = 'off'; end
if strcmpi(opt.plotmode, 'condensed'), opt.factor1grouped = 'on';  opt.factor2grouped = 'on';  end
if strcmpi(opt.plotsubjects, 'on'),    opt.factor1grouped = 'off'; opt.factor2grouped = 'off'; end
if strcmpi(opt.plotsubjects, 'on')
    opt.factor1grouped  = 'apart';
    opt.factor2grouped = 'apart';
end
if strcmpi(opt.effect, 'main') && ncori > 1 && nc2ori > 1 && (strcmpi(opt.factor1stat, 'on') || strcmpi(opt.factor2stat, 'on'))
    opt.factor1grouped = 'apart';
    opt.factor2grouped = 'apart';
end

if strcmpi(opt.factor1grouped,  'on'), opt.factor1grouped =  'together'; end
if strcmpi(opt.factor2grouped, 'on'), opt.factor2grouped = 'together'; end

if strcmpi(opt.factor1grouped, 'together') &&  strcmpi(opt.factor2stat, 'on') && length(opt.factor2vals) > 1, opt.factor1grouped = 'apart'; end
if strcmpi(opt.factor2grouped, 'together') &&  strcmpi(opt.factor1stat, 'on') && length(opt.factor1vals) > 1, opt.factor2grouped = 'apart'; end
if ~( strcmpi(opt.factor1grouped,  'together') && strcmpi(opt.factor2grouped, 'together') )
    if strcmpi(opt.factor1grouped, 'together') && ~isempty(opt.factor1vals), alllegends = opt.factor1vals; opt.factor1vals = ''; end
    if strcmpi(opt.factor2grouped, 'together') && ~isempty(opt.factor2vals), alllegends = opt.factor2vals; opt.factor2vals = ''; end
end

if ~iscell(opt.valsunit),   opt.valsunit   = { opt.valsunit }; end
if ~iscell(opt.chanlabels), opt.chanlabels = { opt.chanlabels }; end
if ~iscell(opt.factor1vals), opt.factor1vals  = { opt.factor1vals }; end
if ~iscell(opt.factor2vals), opt.factor2vals = { opt.factor2vals }; end
if isempty(opt.factor1vals), opt.factor1vals{1}  = ''; end
if isempty(opt.factor2vals), opt.factor2vals{1} = ''; end

for c1 = 1:length(opt.factor1vals)
    for c2 = 1:length(opt.factor2vals)

        % value (ms or Hz)
        % ----------------
        fig_title1 = '';
        for i = 1:length(opt.vals)
            if ~isempty(opt.vals{i}) && ~isnan(opt.vals{i}(1)) 
                if opt.vals{i}(1) == opt.vals{i}(end), fig_title1 = [ num2str(opt.vals{i}(1)) '' opt.valsunit{i} fig_title1];
                else                                   fig_title1 = [ num2str(opt.vals{i}(1)) '-' num2str(opt.vals{i}(2)) '' opt.valsunit{i} fig_title1 ];
                end
                if length(opt.vals) > i, fig_title1 = [ ' & ' fig_title1 ]; end
            end
        end

        % conditions
        % ----------
        if ~isempty(opt.factor1vals{c1})
            fig_title1 = [ value2str(opt.factor1vals{c1}) ', ' fig_title1];
        end
        if ~isempty(opt.factor2vals{c2})
            fig_title1 = [ value2str(opt.factor2vals{c2}) ', ' fig_title1];
        end

        % channel labels, component name, subject name and datatype
        % ---------------------------------------------------------
        fig_title2 = '';
        if length( opt.chanlabels ) == 1 && ~isempty( opt.chanlabels{1} )
            fig_title2 = [ opt.chanlabels{1} ', ' fig_title2 ];
        end    
        
        % cluster and component name
        % --------------------------
        if ~isempty( opt.clustname )
            if ~isempty( opt.compnames )
                if iscell( opt.compnames )
                     compstr = [ 'C' num2str(opt.compnames{min(c1,size(opt.compnames,1)),min(c2,size(opt.compnames,2))}) ];
                else compstr = [ opt.compnames ];
                end
            else compstr = '';
            end
            if ~isempty( opt.subject )
                if ~isempty( compstr )
                     fig_title2 = [ opt.subject '/' compstr ', ' fig_title2 ];
                else fig_title2 = [ opt.subject ', ' fig_title2 ];
                end
            elseif ~isempty( compstr )
                fig_title2 = [ compstr ', ' fig_title2 ];
            end
            
            if ~isempty( opt.datatype )
                 fig_title2 = [ opt.clustname ' ' opt.datatype ', ' fig_title2 ];
            else fig_title2 = [ opt.clustname ', ' fig_title2 ];
            end
            
        else
            if ~isempty( opt.compnames )
                if iscell( opt.compnames )
                else fig_title2 = [ 'C' num2str(opt.compnames{c1,c2}) ', ' fig_title2 ];
                     fig_title2 = [ opt.compnames        ', ' fig_title2 ];
                end
            end
            
            % subject and data type
            % ---------------------
            if ~isempty( opt.subject )
                if ~isempty( opt.datatype )
                     fig_title2 = [ opt.subject ' ' opt.datatype ', ' fig_title2 ];
                else fig_title2 = [ opt.subject ', ' fig_title2 ];
                end
            elseif ~isempty( opt.datatype )
                fig_title2 = [ opt.datatype ' - ' fig_title2 ];
            end
        end    
        
        if strcmpi(opt.factor2grouped, 'together') && strcmpi(opt.factor1grouped, 'together')
            fig_title = fig_title2;
            if ~isempty(fig_title1) && length(fig_title1) > 1 && strcmpi(fig_title1(end-1:end), ', '), fig_title1(end-1:end) = []; end
            if ~isempty(fig_title1) && length(fig_title1) > 1 && strcmpi(fig_title1(end-1:end), '- '), fig_title1(end-1:end) = []; end
            alllegends{c1, c2} = fig_title1;
        else
            fig_title = [ fig_title2 fig_title1 ];
        end
        if ~isempty(fig_title) && length(fig_title) > 1 && strcmpi(fig_title(end-1:end), ', '), fig_title(end-1:end) = []; end
        if ~isempty(fig_title) && length(fig_title) > 1 && strcmpi(fig_title(end-1:end), '- '), fig_title(end-1:end) = []; end
        all_titles{c1,c2}  = fig_title;
        
    end
end
if ~isempty(alllegends)
    alllegends = alllegends';
    alllegends = alllegends(:)';

    % convert legends to string if necessary
    % --------------------------------------
    for ileg = 1:length(alllegends), 
        alllegends{ileg} = value2str(alllegends{ileg}); 
    end
end

% statistic titles
% ----------------
if isnan(opt.threshold), 
    basicstat = '(p-value)';
else
    if length(opt.threshold) >= 1
        basicstat = sprintf([ '(p<%.' thresh_pres(opt.threshold(1)) 'f)'], opt.threshold(1));
    end
    if length(opt.threshold) >= 2
        basicstat = [ basicstat(1:end-1) sprintf([' red;p<%.' thresh_pres(opt.threshold(2))  'f brown)'], opt.threshold(2)) ];
    end
    if length(opt.threshold) >= 3
        basicstat = [ basicstat(1:end-1) sprintf([';\np<%.' thresh_pres(opt.threshold(3)) 'f black)' ], opt.threshold(3)) ];  
    end
end    
if ~isempty(opt.statistics), basicstat = [ basicstat ' ' opt.statistics    ]; end
if ~isempty(opt.mcorrect) && ~strcmpi(opt.mcorrect, 'none'),   basicstat = [ basicstat ' with ' opt.mcorrect ]; end
rown = size(all_titles,1);
if strcmpi(opt.factor1stat, 'on')
    rown = rown+1;
    if strcmpi(opt.effect, 'marginal')
        for c2 = 1:length(opt.factor2vals)
            all_titles{rown, c2} = [ value2str(opt.factor2vals{c2}) ' ' basicstat ];
        end
    else
        all_titles{rown, 1} = [ opt.factor1 ' ' basicstat ]; 
    end
end
coln = size(all_titles,2);
if strcmpi(opt.factor2stat, 'on')
    coln = coln+1;
    if strcmpi(opt.effect, 'marginal')
        for c1 = 1:length(opt.factor1vals)
            all_titles{c1, coln} = [ value2str(opt.factor1vals{c1}) ' ' basicstat ];   
        end
    else
        all_titles{1, coln} = [ opt.factor2 ' ' basicstat ];   
    end
end
if ~strcmpi(opt.effect, 'marginal') && (strcmpi(opt.factor1stat, 'on') && strcmpi(opt.factor2stat, 'on') && rown > 1 && coln > 1)
   all_titles{rown, coln} = [ 'Interaction ' basicstat ];   
end

function pres = thresh_pres(thresh_pres);
    if (round(thresh_pres*100)-thresh_pres*100) == 0
        pres = 2;
    elseif (round(thresh_pres*1000)-thresh_pres*1000) == 0
        pres = 3;
    else
        pres = -round(log10(thresh_pres))+1;
    end
    pres = num2str(pres);

% convert 
% -------
function str = value2str(value)
    
    if ischar(value)
        str = value;
    elseif isnumeric(value)
        if length(value) == 1
            str = num2str(value);
        else
            str = num2str(value(1));
            if length(value) <= 5
                for ind = 2:length(value)
                    str = [ str ' & ' num2str(value(ind)) ];
                end
            else
                str = [ str ' & ' num2str(value(2)) ' & ...' ];
            end
        end
    else % cell array
        str = value{1};
        for ind = 2:length(value)
            str = [ str ' & ' value{ind} ];
        end
    end
    
            


