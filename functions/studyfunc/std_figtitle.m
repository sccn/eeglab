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
%  'condnames'   - [cell array] names of conditions
%  'cond2names'  - [cell array] names of conditions
%  'vals'        - [cell array or real array] value for plot panel (for 
%                  example, 0 ms)
%  'valsunit'    - [cell array or string] unit value (i.e. ms)
%
% Optional statistical titles:
%  'condstat'    - ['on'|'off'] default is 'off'
%  'cond2stat'   - ['on'|'off'] default is 'off'
%  'mcorrect'    - [string] correction for multiple comparisons. Default is
%                  empty.
%  'statistics'  - [string] type of statictics
%  'threshold'   - [real] treshold value
% 
% Optional grouping:
%  'condgroup'   - ['on'|'off'] group conditions together default is 'off'
%  'cond2group'  - ['on'|'off'] default is 'off'
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

function [all_titles alllegends ] = std_figtitle(varargin)

if nargin < 1
    help std_figtitle;
    return;
end;

alllegends = {};
opt = finputcheck( varargin, { 'chanlabels'  {'cell','string'}   []     {};
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
                               'datatype'    'string'            []     '';
                               'clustname'   'string'            []     '';
                               'compnames'   {'cell','string'}   []     {};
                               'vals'        {'cell','real'}     []     {}; % just for titles
                               'valsunit'    {'cell','string'}   []     {}; % just for titles
                               'subject'     'string'            []              '' }, 'std_figtitle'); %, 'ignore');
if isstr(opt), error(opt); end;
if ~iscell(opt.vals),       opt.vals       = { opt.vals }; end;
ncori  = length(opt.condnames);
nc2ori = length(opt.cond2names);
if ~isempty(opt.vals) && ~isempty(opt.vals{1}) && ~isnan(opt.vals{1}(1)), opt.condgroup = 'off'; opt.cond2group = 'off'; end;
if strcmpi(opt.plotmode, 'condensed'), opt.condgroup = 'on';  opt.cond2group = 'on';  end;
if strcmpi(opt.plotsubjects, 'on'),    opt.condgroup = 'off'; opt.cond2group = 'off'; end;
if strcmpi(opt.plotsubjects, 'on')
    opt.condgroup  = 'apart';
    opt.cond2group = 'apart';
end;
if strcmpi(opt.condgroup,  'on'), opt.condgroup =  'together'; end;
if strcmpi(opt.cond2group, 'on'), opt.cond2group = 'together'; end;
    
if strcmpi(opt.condgroup,  'together') &&  strcmpi(opt.cond2stat, 'on'),  opt.condgroup  = 'apart'; end;
if strcmpi(opt.cond2group, 'together') &&  strcmpi(opt.condstat , 'on') , opt.cond2group = 'apart'; end;
if ~( strcmpi(opt.condgroup,  'together') && strcmpi(opt.cond2group, 'together') )
    if strcmpi(opt.condgroup,  'together'), alllegends = opt.condnames ; opt.condnames  = ''; end;
    if strcmpi(opt.cond2group, 'together'), alllegends = opt.cond2names; opt.cond2names = ''; end;
end;
if ~iscell(opt.valsunit),   opt.valsunit   = { opt.valsunit }; end;
if ~iscell(opt.chanlabels), opt.chanlabels = { opt.chanlabels }; end;
if ~iscell(opt.condnames),  opt.condnames  = { opt.condnames }; end;
if ~iscell(opt.cond2names), opt.cond2names = { opt.cond2names }; end;
if isempty(opt.condnames),  opt.condnames{1}  = ''; end;
if isempty(opt.cond2names), opt.cond2names{1} = ''; end;

for c1 = 1:length(opt.condnames)
    for c2 = 1:length(opt.cond2names)

        % value (ms or Hz)
        % ----------------
        fig_title1 = '';
        for i = 1:length(opt.vals)
            if ~isempty(opt.vals{i}) && ~isnan(opt.vals{i}(1)) 
                if opt.vals{i}(1) == opt.vals{i}(end), fig_title1 = [ num2str(opt.vals{i}(1)) '' opt.valsunit{i} fig_title1];
                else                                   fig_title1 = [ num2str(opt.vals{i}(1)) '-' num2str(opt.vals{i}(2)) '' opt.valsunit{i} fig_title1 ];
                end;
                if length(opt.vals) > i, fig_title1 = [ ' & ' fig_title1 ]; end;
            end;
        end;

        % conditions
        % ----------
        if ~isempty(opt.condnames{c1})
            fig_title1 = [ value2str(opt.condnames{c1}) ', ' fig_title1];
        end;
        if ~isempty(opt.cond2names{c2})
            fig_title1 = [ value2str(opt.cond2names{c2}) ', ' fig_title1];
        end;

        % channel labels, component name, subject name and datatype
        % ---------------------------------------------------------
        fig_title2 = '';
        if length( opt.chanlabels ) == 1 && ~isempty( opt.chanlabels{1} )
            fig_title2 = [ opt.chanlabels{1} ', ' fig_title2 ];
        end;    
        
        % cluster and component name
        % --------------------------
        if ~isempty( opt.clustname )
            if ~isempty( opt.compnames )
                if iscell( opt.compnames )
                     compstr = [ 'C' num2str(opt.compnames{min(c1,size(opt.compnames,1)),min(c2,size(opt.compnames,2))}) ];
                else compstr = [ opt.compnames ];
                end;
            else compstr = '';
            end;
            if ~isempty( opt.subject )
                if ~isempty( compstr )
                     fig_title2 = [ opt.subject '/' compstr ', ' fig_title2 ];
                else fig_title2 = [ opt.subject ', ' fig_title2 ];
                end;
            elseif ~isempty( compstr )
                fig_title2 = [ compstr ', ' fig_title2 ];
            end;
            
            if ~isempty( opt.datatype )
                 fig_title2 = [ opt.clustname ' ' opt.datatype ', ' fig_title2 ];
            else fig_title2 = [ opt.clustname ', ' fig_title2 ];
            end;
            
        else
            if ~isempty( opt.compnames )
                if iscell( opt.compnames )
                else fig_title2 = [ 'C' num2str(opt.compnames{c1,c2}) ', ' fig_title2 ];
                     fig_title2 = [ opt.compnames        ', ' fig_title2 ];
                end;
            end;
            
            % subject and data type
            % ---------------------
            if ~isempty( opt.subject )
                if ~isempty( opt.datatype )
                     fig_title2 = [ opt.subject ' ' opt.datatype ', ' fig_title2 ];
                else fig_title2 = [ opt.subject ', ' fig_title2 ];
                end;
            elseif ~isempty( opt.datatype )
                fig_title2 = [ opt.datatype ' - ' fig_title2 ];
            end;
        end;    
        
        if strcmpi(opt.cond2group, 'together') && strcmpi(opt.condgroup, 'together')
            fig_title = fig_title2;
            if ~isempty(fig_title1) && length(fig_title1) > 1 && strcmpi(fig_title1(end-1:end), ', '), fig_title1(end-1:end) = []; end;
            if ~isempty(fig_title1) && length(fig_title1) > 1 && strcmpi(fig_title1(end-1:end), '- '), fig_title1(end-1:end) = []; end;
            alllegends{c1, c2} = fig_title1;
        else
            fig_title = [ fig_title2 fig_title1 ];
        end;
        if ~isempty(fig_title) && length(fig_title) > 1 && strcmpi(fig_title(end-1:end), ', '), fig_title(end-1:end) = []; end;
        if ~isempty(fig_title) && length(fig_title) > 1 && strcmpi(fig_title(end-1:end), '- '), fig_title(end-1:end) = []; end;
        all_titles{c1,c2}  = fig_title;
        
    end;
end;
if ~isempty(alllegends)
    alllegends = alllegends';
    alllegends = alllegends(:)';

    % convert legends to string if necessary
    % --------------------------------------
    for ileg = 1:length(alllegends), 
        alllegends{ileg} = value2str(alllegends{ileg}); 
    end;
end;

% statistic titles
% ----------------
if isnan(opt.threshold), 
    basicstat = '(p-value)';
else
    if length(opt.threshold) >= 1
        basicstat = sprintf([ '(p<%.' thresh_pres(opt.threshold(1)) 'f)'], opt.threshold(1));
    end;
    if length(opt.threshold) >= 2
        basicstat = [ basicstat(1:end-1) sprintf([' red;p<%.' thresh_pres(opt.threshold(2))  'f brown)'], opt.threshold(2)) ];
    end;
    if length(opt.threshold) >= 3
        basicstat = [ basicstat(1:end-1) sprintf([';\np<%.' thresh_pres(opt.threshold(3)) 'f black)' ], opt.threshold(3)) ];  
    end;
end;    
if ~isempty(opt.statistics), basicstat = [ basicstat ' ' opt.statistics    ]; end;
if ~isempty(opt.mcorrect) && ~strcmpi(opt.mcorrect, 'none'),   basicstat = [ basicstat ' with ' opt.mcorrect ]; end;
if strcmpi(opt.condstat, 'on')
    rown = size(all_titles,1)+1;
    for c2 = 1:length(opt.cond2names)
        all_titles{rown, c2} = [ value2str(opt.cond2names{c2}) ' ' basicstat ];   
    end;
end;
if strcmpi(opt.cond2stat, 'on')
    coln = size(all_titles,2)+1;
    for c1 = 1:length(opt.condnames)
        all_titles{c1, coln} = [ value2str(opt.condnames{c1}) ' ' basicstat ];   
    end;
end;
if strcmpi(opt.condstat, 'on') && strcmpi(opt.cond2stat, 'on')
   all_titles{rown, coln} = [ 'Interaction ' basicstat ];   
end;

function pres = thresh_pres(thresh_pres);
    if (round(thresh_pres*100)-thresh_pres*100) == 0
        pres = 2;
    elseif (round(thresh_pres*1000)-thresh_pres*1000) == 0
        pres = 3;
    else
        pres = -round(log10(thresh_pres))+1;
    end;
    pres = num2str(pres);

% convert 
% -------
function str = value2str(value)
    
    if isstr(value)
        str = value;
    elseif isnumeric(value)
        if length(value) == 1
            str = num2str(value);
        else
            str = num2str(value(1));
            if length(value) <= 5
                for ind = 2:length(value)
                    str = [ str ' & ' num2str(value(ind)) ];
                end;
            else
                str = [ str ' & ' num2str(value(2)) ' & ...' ];
            end;
        end;
    else % cell array
        str = value{1};
        for ind = 2:length(value)
            str = [ str ' & ' value{ind} ];
        end;
    end;
    
            


