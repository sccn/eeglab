% std_readspec() - load spectrum measures for data channels or 
%                  for all components of a specified cluster.
%                  Called by plotting functions
%                  std_envtopo(), std_erpplot(), std_erspplot(), ...
% Usage:
%         >> [STUDY, specdata, allfreqs, setinds, cinds] = ...
%                   std_readspec(STUDY, ALLEEG, varargin);
% Inputs:
%       STUDY - studyset structure containing some or all files in ALLEEG
%      ALLEEG - vector of loaded EEG datasets
%
% Optional inputs:
%  'channels'  - [cell] list of channels to import {default: all}
%  'clusters'  - [integer] list of clusters to import {[]|default: all but
%                the parent cluster (1) and any 'NotClust' clusters}
%  'freqrange' - [min max] frequency range {default: whole measure range}
%  'rmsubjmean' - ['on'|'off'] remove mean subject spectrum from every
%                 channel spectrum, making them easier to compare
%                 { default: 'off' }
%  'subject'    - [string] select a specific subject {default:all}
%  'component'  - [integer] select a specific component in a cluster
%                 {default:all}
%  'singletrials' - ['on'|'off'] load single trials spectral data (if available)
%
% Output:
%  STUDY    - updated studyset structure
%  specdata - [cell array] spectral data (the cell array size is 
%             condition x groups)
%  freqs    - [float array] array of frequencies
%  setinds  - [cell array] datasets indices
%  cinds    - [cell array] channel or component indices
%
% Author: Arnaud Delorme, CERCO, 2006-

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

% $Log: not supported by cvs2svn $
% Revision 1.18  2007/02/05 16:17:13  arno
% fix crash for old study subtracting average spectrum
%
% Revision 1.17  2006/11/23 01:12:50  arno
% implement mean spectrum subtraction
%
% Revision 1.16  2006/11/10 00:08:43  arno
% reprogram for channel
%
% Revision 1.15  2006/10/04 23:39:49  toby
% Bug fix courtesy Bas de Kruif
%
% Revision 1.14  2006/03/28 15:38:13  scott
% help msg
%
% Revision 1.13  2006/03/14 02:32:32  scott
% help msg
%
% Revision 1.12  2006/03/11 07:30:01  arno
% freqrange input
%
% Revision 1.11  2006/03/11 07:25:37  arno
% header
%
% Revision 1.10  2006/03/10 16:33:37  arno
% selecting frequency range for reading
%
% Revision 1.9  2006/03/10 00:37:45  arno
% error msg
%
% Revision 1.8  2006/03/09 18:10:38  arno
% *** empty log message ***
%
% Revision 1.7  2006/03/09 18:10:18  arno
% do not use etc field any more
%
% Revision 1.6  2006/03/09 00:42:09  arno
% fix reading file
%
% Revision 1.5  2006/03/09 00:37:31  arno
% now writing matlab fileend
%
% Revision 1.4  2006/03/09 00:03:57  arno
% read spectrum form matlab file
%
% Revision 1.3  2006/03/08 21:06:37  arno
% rename func
%
% Revision 1.2  2006/03/07 22:21:12  arno
% use fullfile
%

function [STUDY, specdata, allfreqs] = std_readspec(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_readspec;
    return;
end
if ~isstruct(ALLEEG) % old calling format
    [STUDY, specdata, allfreqs] = std_readspecsub(STUDY, ALLEEG, varargin{:});
    return;
end;

[STUDY, specdata, allfreqs] = std_readerp(STUDY, ALLEEG, 'datatype', 'spec', varargin{:});
return;

STUDY = pop_specparams(STUDY, 'default');
[opt moreopts] = finputcheck( varargin, { ...
    'type'          { 'string' 'cell' } { [] [] } '';
    'channels'      'cell'    []             {};
    'clusters'      'integer' []             [];
    'freqrange'     'real'    []             STUDY.etc.specparams.freqrange;
    'rmsubjmean'    'string'  { 'on' 'off' } 'off';
    'singletrials'  'string'  { 'on' 'off' } 'off';
    'component'     'integer' []             [];
    'subject'       'string'  []             '' }, ...
    'std_readspec', 'ignore');
if isstr(opt), error(opt); end;
nc = max(length(STUDY.condition),1);
ng = max(length(STUDY.group),1);

% find channel indices
% --------------------
if ~isempty(opt.channels)
     finalinds = std_chaninds(STUDY, opt.channels);
else finalinds = opt.clusters;
end;

for ind = 1:length(finalinds)

    % find indices
    % ------------
    if ~isempty(opt.channels)
        tmpstruct = STUDY.changrp(finalinds(ind));
        allinds       = tmpstruct.allinds;
        for i=1:length(allinds(:)), allinds{i} = -allinds{i}; end; % invert sign for reading
        setinds       = tmpstruct.setinds;
    else
        [ tmpstruct setinds allinds ] = std_setcomps2cell(STUDY, finalinds(ind));
    end;

    % check if data is already here
    % -----------------------------
    dataread = 0;
    if strcmpi(opt.singletrials,'off') && isfield(tmpstruct, 'specdata') && ...
            isequal( STUDY.etc.specparams.freqrange, opt.freqrange) && ~isempty(tmpstruct.specdata)
        dataread = 1;
    end;
    if strcmpi(opt.singletrials,'on') && isfield(tmpstruct, 'specdatatrials') && ...
            isequal( STUDY.etc.specparams.freqrange, opt.freqrange) && ~isempty(tmpstruct.specdatatrials)
        if ~isempty(opt.channels) && strcmpi(tmpstruct.spectrialinfo, opt.subject)
            dataread = 1; 
        elseif isequal(tmpstruct.spectrialinfo, opt.component) 
            dataread = 1; 
        end;
    end;
    
    if ~dataread
        % reserve arrays
        % --------------
        allspec  = cell( max(length(STUDY.condition),1), max(length(STUDY.group),1) );
        filetype = 'spec';
        tmpind = 1; while(isempty(setinds{tmpind})), tmpind = tmpind+1; end;
        try,
            [ tmp allfreqs datatrialpresent] = std_readspecsub( ALLEEG, setinds{tmpind}(1), allinds{tmpind}(1), opt.freqrange, 0, fastif(strcmpi(opt.singletrials, 'on'),1,0));
        catch
            filetype = 'ersp';
            datatrialpresent = 0;
            disp('Cannot find spectral file, trying ERSP baseline file instead');
            [ tmpersp allfreqs alltimes tmpparams tmpspec] = std_readersp( ALLEEG, setinds{tmpind}(1), allinds{tmpind}(1), [], opt.freqrange);
        end;
        for c = 1:nc
            for g = 1:ng
                allspec{c, g} = repmat(single(0), [length(allfreqs), length(allinds{c,g}) ]);
            end;
        end;

        % read the data and select channels
        % ---------------------------------
        fprintf('Reading Spectrum data...');
        if strcmpi(opt.singletrials, 'on')
            if ~datatrialpresent
                fprintf('\n');
                errordlg2('No single trial data - recompute data files');
                specdata = [];
                return;
            end;
            allsubjects = { STUDY.datasetinfo(:).subject };
            for c = 1:nc
                for g = 1:ng
                    if ~isempty(opt.channels)
                        if ~isempty(opt.subject) inds = strmatch( opt.subject, allsubjects(setinds{c,g}));
                        else inds = 1:length(allinds{c,g}); end;
                    else
                        if ~isempty(opt.component) inds = find( allinds{c,g} == opt.component);
                        else inds = 1:length(allinds{c,g}); end;
                    end;
                    allspec{c, g} = squeeze(std_readspecsub( ALLEEG, setinds{c,g}(inds), allinds{c,g}(inds), opt.freqrange, 0, 1));
                end;
            end;
        else
            if strcmpi(filetype, 'spec')
                for c = 1:nc
                    for g = 1:ng
                        allspec{c, g} = std_readspecsub( ALLEEG, setinds{c,g}(:), allinds{c,g}(:), opt.freqrange)';
                    end;
                end;
            else % std_readersp cannot be converted to read multiple datasets since it subtracts data between conditions
                for c = 1:nc
                    for g = 1:ng
                        for indtmp = 1:length(allinds{c,g})
                            [ tmpersp allfreqs alltimes tmpparams tmpspec] = std_readersp( ALLEEG, setinds{c,g}(indtmp), allinds{c,g}(indtmp), [], opt.freqrange);
                            allspec{c, g}(:,indtmp) = 10*log(tmpspec(:));
                            fprintf('.');
                        end;
                    end;
                end;
            end;
        end;
        fprintf('\n');
        
        % remove mean of each subject across groups and conditions
        if strcmpi(opt.rmsubjmean, 'on') && ~isempty(opt.channels) && strcmpi(opt.singletrials, 'off')
            disp('Removing mean spectrum accross subjects');
            for indtmp = 1:length(allinds{c,g}) % scan subjects
               meanspec =zeros(size( allspec{1, 1}(:,indtmp) ));
               for c = 1:nc
                    for g = 1:ng
                        meanspec = meanspec + allspec{c, g}(:,indtmp)/(nc*ng);
                    end;
               end;
               for c = 1:nc
                    for g = 1:ng
                        allspec{c, g}(:,indtmp) = allspec{c, g}(:,indtmp) - meanspec; % subtractive model
                        % allspec{c, g}(:,indtmp) = allspec{c, g}(:,indtmp)./meanspec; % divisive model
                    end;
               end;
            end;
        end;

        if strcmpi(opt.singletrials, 'on')
             tmpstruct.specdatatrials = allspec;
             if ~isempty(opt.channels)
                  tmpstruct.spectrialinfo  = opt.subject;
             else tmpstruct.spectrialinfo  = opt.component;
             end;
        else tmpstruct.specdata = allspec;
        end;
        tmpstruct.specfreqs = allfreqs;
        
        % copy results to structure
        % -------------------------
        fieldnames = { 'specdata' 'specfreqs' 'allinds' 'setinds' 'specdatatrials' 'spectrialinfo' };
        for f = 1:length(fieldnames)
            if isfield(tmpstruct, fieldnames{f}),
                tmpdata = getfield(tmpstruct, fieldnames{f});
                if ~isempty(opt.channels)
                     STUDY.changrp = setfield(STUDY.changrp, { finalinds(ind) }, fieldnames{f}, tmpdata);
                else STUDY.cluster = setfield(STUDY.cluster, { finalinds(ind) }, fieldnames{f}, tmpdata);
                end;
            end;
        end;
    end;
end;

% return structure
% ----------------
allinds   = finalinds;
if ~isempty(opt.channels)
    structdat = STUDY.changrp;
    specdata = cell(nc, ng);
    for ind =  1:length(specdata(:))
        if strcmpi(opt.singletrials, 'on')
             specdata{ind} = zeros([ size(structdat(allinds(1)).specdatatrials{ind}) length(allinds)]);
        else specdata{ind} = zeros([ size(structdat(allinds(1)).specdata{ind}) length(allinds)]);
        end;
        for index = 1:length(allinds)
            if strcmpi(opt.singletrials, 'on')
                 specdata{ind}(:,:,index) = structdat(allinds(index)).specdatatrials{ind};
            else specdata{ind}(:,:,index) = structdat(allinds(index)).specdata{ind};
            end;
            allfreqs                 = structdat(allinds(index)).specfreqs;
            compinds                 = structdat(allinds(index)).allinds;
            setinds                  = structdat(allinds(index)).setinds;
        end;
        specdata{ind} = squeeze(permute(specdata{ind}, [1 3 2])); % time elec subjects
    end;
    if ~isempty(opt.subject) && strcmpi(opt.singletrials,'off')
        specdata = std_selsubject(specdata, opt.subject, setinds, { STUDY.datasetinfo(:).subject }, 2); 
    end;
else
    if strcmpi(opt.singletrials, 'on')
         specdata = STUDY.cluster(allinds(1)).specdatatrials;
    else specdata = STUDY.cluster(allinds(1)).specdata;
    end;
    allfreqs = STUDY.cluster(allinds(1)).specfreqs;
    compinds = STUDY.cluster(allinds(1)).allinds;
    setinds  = STUDY.cluster(allinds(1)).setinds;
    if ~isempty(opt.component) && length(allinds) == 1 && strcmpi(opt.singletrials,'off')
        specdata = std_selcomp(STUDY, specdata, allinds, setinds, compinds, opt.component);
    end;
end;

% std_readspecsub() - returns the stored mean power spectrum for an ICA component 
%                  in a specified dataset.  The spectrum is assumed to have been 
%                  saved in a Matlab file, "[dataset_name].icaspec", in the same
%                  directory as the dataset file. If this file doesn't exist,
%                  use std_spec() to create it or a pre-clustering function
%                  (pop_preclust() or std_preclust()) that calls it. 
% Usage:    
%  >> [spec, freqs] = std_readspecsub(ALLEEG, setindx, component, freqrange, rmsubjmean);  
%
% Inputs:
%   ALLEEG     - a vector of dataset EEG structures (may also be one dataset). 
%                Must contain the dataset of interest (the 'setindx' below).
%   setindx    - [integer] an index of an EEG dataset in the ALLEEG
%                structure for which to read a component spectrum.
%   component  - [integer] index of the component in the selected EEG dataset 
%                for which to return the spectrum
%   freqrange  - [min max in Hz] frequency range to return
%   rmsubjmean - [0|1] remove subject mean spectrum (0 is no and is the default)
%
% Outputs:
%   spec      - the log-power spectrum of the requested ICA component in the
%               specified dataset (in dB)
%   freqs     - vector of spectral frequencies (in Hz)
%
%  See also  std_spec(), pop_preclust(), std_preclust()
%
% Authors:  Arnaud Delorme, Hilit Serby, SCCN, INC, UCSD, February, 2005

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

% $Log: not supported by cvs2svn $

function [X, f, singletrialdatapresent] = std_readspecsub(ALLEEG, abset, comp, freqrange, rmsubjmean, singletrial);

if nargin < 4
    freqrange = [];
end;
if nargin < 5
    rmsubjmean = 0;
end;
if nargin < 6
    singletrial = 0;
end;
    
% if spectrum is not available but ERSP is, use ERSP
% --------------------------------------------------
tmpfile = fullfile( ALLEEG(abset(1)).filepath,[ ALLEEG(abset(1)).filename(1:end-3) 'datspec']);
if ~exist(tmpfile) && ~singletrial
    disp('Cannot find spectral file, trying ERSP baseline file instead');
    for indtmp = 1:length(abset)
        [ tmpersp f alltimes tmpparams tmpspec] = std_readersp( ALLEEG, abset(indtmp), comp(indtmp), [], opt.freqrange);
        if indtmp == 1, X = repmat(single(0), [ size(tmpspec) length(comp) ]); end;
        X(:,indtmp) = 10*log(tmpspec(:));
        fprintf('.');
    end;
    return;
end;

X = [];
if length(abset) < length(comp)
    abset = ones(1,length(comp))*abset;
end;

% convert components or channel indices
% -------------------------------------
if iscell(comp)
    % find channel indices list
    % -------------------------
    chanind  = [];
    chanlabs = lower({ ALLEEG(abset(1)).chanlocs.labels });
    for index = 1:length(comp)
        tmp = strmatch(lower(comp{index}), chanlabs, 'exact');
        if isempty(tmp)
            error([ 'Channel ''' comp{index} ''' not found in dataset ' int2str(abset)]);
        else    
            chanind = [ chanind tmp ];
        end;
    end;
    prefix = 'chan';
    inds   = chanind;
elseif comp(1) < 0
    prefix = 'chan';
    inds   = -comp;
else
    prefix = 'comp';
    inds   = comp;
end;

singletrialdatapresent = 1;
for k = 1:length(abset)

    if strcmpi(prefix, 'chan')
         filename = fullfile( ALLEEG(abset(k)).filepath,[ ALLEEG(abset(k)).filename(1:end-3) 'datspec']);
    else filename = fullfile( ALLEEG(abset(k)).filepath,[ ALLEEG(abset(k)).filename(1:end-3) 'icaspec']);
    end;

    try,
        warning backtrace off;
        if rmsubjmean == 0
             erpstruct = load( '-mat', filename, [ prefix int2str(inds(k)) ], 'freqs' );
        else erpstruct = load( '-mat', filename, [ prefix int2str(inds(k)) ], 'freqs', 'average_spec' );
        end;
        warning backtrace on;
    catch
        error( [ 'Cannot read file ''' filename '''' ]);
    end;

    tmpdat    = getfield(erpstruct, [ prefix int2str(inds(k)) ]);
    if singletrial == 0,
        if size(tmpdat,2) > 1 && size(tmpdat,1) > 1, tmpdat = mean(tmpdat,2); end;
    else
        if size(tmpdat,1) == 1 || size(tmpdat,2) == 1
            singletrialdatapresent = 0;
        end;
    end;
    if rmsubjmean, 
        if isfield(erpstruct, 'average_spec')
            if singletrial == 0
                 tmpdat = tmpdat - erpstruct.average_spec; 
            else tmpdat = tmpdat - repmat(erpstruct.average_spec', [1 size(tmpdat,2)]); 
            end;
        end;
    end;
    if singletrial
        if k == 1
            if size(tmpdat,1) == 1, tmpdat = tmpdat'; end;
            X = zeros([ 1 size(tmpdat) ]);
            X(1,:,:) = tmpdat;
        else
            X(1,:,end+1:end+size(tmpdat,2)) = tmpdat;
        end;
    else
        if k == 1
            if size(tmpdat,1) == 1, tmpdat = tmpdat'; end;
            X = zeros([ length(comp) size(tmpdat) ]);
        end;
        X(k,:,:) = tmpdat;
    end;
    f = getfield(erpstruct, 'freqs');
end;

% select frequency range of interest
% ----------------------------------
if ~isempty(freqrange)
    maxind = max(find(f <= freqrange(end)));
    minind = min(find(f >= freqrange(1)));
else
    %if not, use whole spectrum
    maxind = length(f);
    minind = 1;
end

f = f(minind:maxind);
X = X(:,minind:maxind,:);
return;
