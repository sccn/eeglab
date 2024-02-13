% STD_READDATA - load measures for data channels or 
%                  for all components of a specified cluster.
%                  Called by plotting functions
%                  STD_ENVTOPO, STD_ERPPLOT, STD_ERSPPLOT, ...
% Usage:
%         >> [STUDY, datavals, xvals, yvals, events, params] = ...
%                   std_readdata(STUDY, ALLEEG, varargin);
% Inputs:
%       STUDY - studyset structure containing some or all files in ALLEEG
%      ALLEEG - vector of loaded EEG datasets
%
% Optional inputs:
%  'design'    - [integer] read files from a specific STUDY design. Default
%                is empty (use current design in STUDY.currentdesign). Use
%                NaN to create a design with with all the data.
%  'channels'  - [cell] list of channels to import {default: none}
%  'clusters'  - [integer] list of clusters to import {[]|default: all but
%                the parent cluster (1) and any 'NotClust' clusters}
%  'singletrials' - ['on'|'off'] load single trials spectral data (if 
%                available). Default is 'off'.
%  'subject'   - [string] select a specific subject {default:all}
%  'datatype'  - ['erp'|'spec'|'ersp'|'itc'|'erpim'] select measure to load
%                Default is 'erp'.
%  'component' - [integer] select a specific component in a cluster.
%                This is the index of the component in the cluster not the
%                component number {default:all}
%  'ndim'       - [integer] number of dimension of output array
%  'customread' - [func] custom function to read data. See example.
%  'customparams' - [cell] custom parameters for function above. See example.
%
% ERP specific optional inputs:
%  'timerange' - [min max] time range {default: whole measure range}
%  'componentpol' - ['on'|'off'] invert ERP component sign based on
%                   scalp map match with component scalp map centroid.
%                   {default:'on'}
%
% Output:
%  STUDY    - updated studyset structure
%  datavals  - [cell array] erp data (the cell array size is 
%             condition x groups)
%  xvals    - [float array] array of first dim values (for example time)
%  yvals    - [float array] array of second dim values (for example frequencies)
%  events   - [cell array] events (corresponding to the data)
%  params   - [struct] structure containing parameters
%  setinds  - [cell array] index of the dataset for each cell member in
%             datavals. Could be ambiguous if datasets are merged in the 
%             datavals array (the setinds will then contain the index of 
%             one of the datasets).
%  datainfo - [cell array] trial and participant info corresponding to the 
%             datavals array.
%
% Important note: This function does not do baseline correction for ERSP. 
%                 To get the baseline corrected data use the function STD_ERSPPLOT
%
% Example:
%  std_precomp(STUDY, ALLEEG, { ALLEEG(1).chanlocs.labels }, 'erp', 'on');
%  [STUDY,erp,times] = std_readdata(STUDY, ALLEEG, 'channels', { ALLEEG(1).chanlocs(1).labels });
%
%  % To read field 'EEG.etc.eegstats.alpha_asymmetry' for each dataset, enter
%  % (note that you need to use the eegstats plugin to compute this value first)
%  [STUDY,aa] = std_readdata(STUDY, ALLEEG, 'customread', 'std_readeegfield', ...
%               'customparams', {{ 'etc', 'eegstats', 'alpha_asymmetry' }}, 'ndim', 1);
% 
% Author: Arnaud Delorme, CERCO, 2006-

% Copyright (C) Arnaud Delorme, arno@salk.edu
%
% This file is part of EEGLAB, see https://urldefense.com/v3/__http://www.eeglab.org__;!!Mih3wA!Hn3hFRnhxnzS1csHGMH6q95jPsdmpLu36WaSxZZ0BSw4smbm8RghssmOC-feg7O1JHd0ncZ3ACn-r8229A$ 
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

function [STUDY, datavals, xvals, yvals, events, params, setinds, datainfo] = std_readdata(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_readdata;
    return;
end

STUDY = pop_erpparams(STUDY, 'default');
STUDY = pop_specparams(STUDY, 'default');
STUDY = pop_erspparams(STUDY, 'default');
[opt, moreopts] = finputcheck( varargin, { ...
    'design'        'integer' []             STUDY.currentdesign;
    'channels'      'cell'    []             {};
    'clusters'      'integer' []             [];
    'ndim'          'integer' []             [];
    'timerange'     'real'    []             [];
    'freqrange'     'real'    []             [];
    'datatype'      'string'  { 'erp','spec' 'ersp' 'itc' 'erpim' 'custom' } 'erp';
    'singletrials'  'string'  { 'on','off' } 'off';
    'componentpol'  'string'  { 'on','off' } 'on';
    'component'     'integer' []             [];
    'customread'    ''        []             [];
    'customparams'  'cell'    []             {};
    'subject'       'string'  []             '' }, ...
    'std_readdata', 'ignore');
if ~isempty(moreopts)
    fprintf(2, 'Warning: std_readdata ignored at least one parameters: "%s"\n', moreopts{1})
end
if ischar(opt), error(opt); end

dtype = opt.datatype;

% get the file extension
% ----------------------
tmpDataType = opt.datatype;
if ~isempty(opt.customread)
    opt.datatype = 'custom';
end
if strcmpi(opt.datatype, 'ersp') || strcmpi(opt.datatype, 'itc')
    tmpDataType = 'timef'; 
    if isempty(opt.timerange), opt.timerange = STUDY.etc.erspparams.timerange;  end
    if isempty(opt.freqrange), opt.freqrange = STUDY.etc.erspparams.freqrange; end
elseif strcmpi(opt.datatype, 'erpim')
    if isempty(opt.timerange), opt.timerange = STUDY.etc.erpimparams.timerange;  end
elseif strcmpi(opt.datatype, 'erp')    
    if isempty(opt.timerange), opt.timerange = STUDY.etc.erpparams.timerange;  end
elseif strcmpi(opt.datatype, 'spec')    
    if isempty(opt.freqrange), opt.freqrange = STUDY.etc.specparams.freqrange; end
end
if ~isempty(opt.channels), fileExt = [ '.dat' tmpDataType ];
else                       fileExt = [ '.ica' tmpDataType ];
end

% list of subjects
% ----------------
allSubjects = { STUDY.datasetinfo.subject };
uniqueSubjects = unique(allSubjects);
STUDY.subject = uniqueSubjects;
if ischar(opt.subject) && ~isempty(opt.subject), subjectList = {opt.subject}; else subjectList = opt.subject; end
if isempty(subjectList)
    if isnan(opt.design), subjectList = STUDY.subject;
    else subjectList = STUDY.design(opt.design).cases.value; 
    end
end

% options
% -------
opts = {};
if ~isempty(opt.timerange), opts = { opts{:} 'timelimits', opt.timerange }; end
if ~isempty(opt.freqrange), opts = { opts{:} 'freqlimits', opt.freqrange }; end
opts = { opts{:} 'singletrials' opt.singletrials };
fprintf('Reading subjects'' data or looking up measure values in EEGLAB cache\n');

% determining component polarity if necessary
% -------------------------------------------
if isempty(opt.customread) && isempty(opt.channels) && strcmpi(dtype, 'erp') && strcmpi(opt.componentpol, 'on')
    componentPol = ones(1, length(STUDY.cluster(opt.clusters).comps)); % default is all 1
    disp('Reading component scalp topo polarities - this is done to invert some ERP component polarities');
    STUDY = std_readtopoclust(STUDY, ALLEEG, opt.clusters);
    componentPol = STUDY.cluster(opt.clusters).topopol;
    if isempty(componentPol)
        disp('Cluster topographies absent - cannot adjust single component ERP polarities');
    end
end

% get all sessions (same code as std_readdat)
% -------------------------------------------
allSessions = { STUDY.datasetinfo.session };
allSessions(cellfun(@isempty, allSessions)) = { 1 };
allSessions = cellfun(@num2str, allSessions, 'uniformoutput', false);
uniqueSessions = unique(allSessions);

for iSubj = 1:length(subjectList)
    fprintf('.');
    
    % check cache
    bigstruct = [];
    if ~isempty(opt.channels), bigstruct.channel = opt.channels;
    else                       bigstruct.cluster = opt.clusters; % there can only be one cluster
    end
    bigstruct.datatype     = opt.datatype;
    bigstruct.singletrials = opt.singletrials;
    bigstruct.subject      = subjectList{iSubj};
    bigstruct.component    = opt.component;
    bigstruct.options      = opts;
    bigstruct.moreoptions  = opt;
    if isnan(opt.design)
         bigstruct.design.variable = struct([]);
    else bigstruct.design.variable = STUDY.design(opt.design).variable;
    end

    % find component indices
    % ----------------------
    if ~isempty(opt.clusters)
        datasetInds = strmatch(subjectList{iSubj}, { STUDY.datasetinfo.subject }, 'exact');
        compList    = [];
        polList     = [];
        if size(STUDY.cluster(opt.clusters).sets,1) ~= length(datasetInds)
            error('Cannot process components from different ICA decompositions of the same subjects'); % sometimes different sessions
        end            
        if isempty(opt.component)
            for iDat = datasetInds(:)'
                indSet   = find(STUDY.cluster(opt.clusters).sets(1,:) == iDat); % each column contain info about the same subject so we many only consider the first row
                if ~isempty(indSet)
                    compList = [ compList STUDY.cluster(opt.clusters).comps(indSet) ];
                    if strcmpi(dtype, 'erp') && strcmpi(opt.componentpol, 'on')
                        polList  = [ polList  componentPol(indSet) ];
                    end
                end
            end
        else
            if ~isempty(intersect(datasetInds, STUDY.cluster(opt.clusters).sets(:,opt.component)))
                compList = [ compList STUDY.cluster(opt.clusters).comps(opt.component) ];
                if strcmpi(dtype, 'erp') && strcmpi(opt.componentpol, 'on')
                    polList  = [ polList  componentPol(opt.component) ];
                end
            end
        end
    end
    
    % read all channels/components at once
    hashcode = gethashcode(std_serialize(bigstruct));
    [STUDY.cache, tmpstruct] = eeg_cache(STUDY.cache, hashcode);
    
    if ~isempty(tmpstruct)
        dataTmp{iSubj}     = tmpstruct{1};
        xvals              = tmpstruct{2};
        yvals              = tmpstruct{3};
        eventsTmp{iSubj}   = tmpstruct{4};
        params             = tmpstruct{5};
        dataTmpSubj{iSubj} = tmpstruct{6};
    else
        datInds = find(strncmp( subjectList{iSubj}, allSubjects, max(cellfun(@length, allSubjects))));

        if isempty(opt.customread)
            fileName = getfilename({STUDY.datasetinfo(datInds).filepath}, STUDY.datasetinfo(datInds(1)).subject, { STUDY.datasetinfo(datInds).session }, fileExt, length(uniqueSessions) == 1);
            if ~isempty(opt.channels)
                 [dataTmp{iSubj}, params, xvals, yvals, eventsTmp{iSubj}, dataTmpSubj{iSubj} ] = std_readfile( fileName, 'designvar', struct(bigstruct.design.variable), opts{:}, 'channels', opt.channels);
            else [dataTmp{iSubj}, params, xvals, yvals, eventsTmp{iSubj}, dataTmpSubj{iSubj} ] = std_readfile( fileName, 'designvar', struct(bigstruct.design.variable), opts{:}, 'components', compList);
            end
        else
            % read custom data
            [dataTmp{iSubj}, params, xvals, yvals, eventsTmp{iSubj}, dataTmpSubj{iSubj}] = feval(opt.customread, STUDY.datasetinfo(datInds), ALLEEG(datInds), struct(bigstruct.design.variable), opt.customparams{:}, opts{:});
        end
        if ~strcmpi(opt.datatype, 'ersp') && ~strcmpi(opt.datatype, 'itc') && ~strcmpi(opt.datatype, 'erpim') % ERP or spectrum or Custom single trials
            % inverting ERP polarity when relevant
            if strcmpi(opt.datatype, 'erp') && ~isempty(opt.clusters) && strcmpi(opt.componentpol, 'on')
                polList = reshape(polList,[1 1 length(polList)]); % components are in the 3rd dim
                dataTmp{iSubj} = cellfun(@(x)bsxfun(@times, x, fastif(isempty(x), [], polList)), dataTmp{iSubj}, 'uniformoutput', false);
            end
            if strcmpi(opt.singletrials, 'off')
                dataTmp{iSubj} = cellfun(@(x)squeeze(mean(x,2)), dataTmp{iSubj}, 'uniformoutput', false); % average
                for iCond = 1:length(dataTmpSubj{iSubj}(:)), if ~isempty(dataTmpSubj{iSubj}{iCond}), dataTmpSubj{iSubj}{iCond} = dataTmpSubj{iSubj}{iCond}(1); end; end
            end
            if strcmpi(opt.datatype, 'spec') && isfield(params, 'logtrials') && strcmpi(params.logtrials, 'off') % if log trial if off it means that single trials are raw power so we need to take the log of the mean
                dataTmp{iSubj} = cellfun(@(x)squeeze(10*log10(x)), dataTmp{iSubj}, 'uniformoutput', false); % average
            end
        elseif strcmpi(opt.datatype, 'erpim')
            %dataTmp{iSubj} = cellfun(@(x)processerpim(x, xvals, params), dataTmp{iSubj}, 'uniformoutput', false);
            for iCond = 1:length(dataTmp{iSubj}(:))
                if all(isnan(eventsTmp{iSubj}{iCond})), eventsTmp{iSubj}{iCond} = []; end
                [dataTmp{iSubj}{iCond}, eventsTmp{iSubj}{iCond}] = processerpim(dataTmp{iSubj}{iCond}, eventsTmp{iSubj}{iCond}, xvals, params);
            end
            nonEmptyCell = find( cellfun(@isempty, dataTmp{iSubj}) == 0);
            if ~isempty(nonEmptyCell)
                yvalsERPim = 1:size(dataTmp{iSubj}{nonEmptyCell(1)},1);
            end
        elseif strcmpi(opt.datatype, 'custom')
            % Nothing to do for custom data, this is done at reading time
        else
            dataTmp{iSubj} = cellfun(@(x)processtf(x, xvals, opt.datatype, opt.singletrials, params), dataTmp{iSubj}, 'uniformoutput', false);
        end
        STUDY.cache = eeg_cache(STUDY.cache, hashcode, { dataTmp{iSubj} xvals yvals eventsTmp{iSubj} params dataTmpSubj{iSubj} });
    end
end
if strcmpi(opt.datatype, 'erpim')
    yvals = yvalsERPim;
end
fprintf('\n');

% if single trials, swap the last 2 dim (put channels before trials)
if strcmpi(opt.singletrials, 'on') && length(opt.channels) > 1
    if ndims(dataTmp{1}{1}) == 3
        for iCase = 1:length(dataTmp)
            for iItem = 1:length(dataTmp{1}(:))
                dataTmp{iCase}{iItem} = permute(dataTmp{iCase}{iItem}, [1 3 2]);
            end
        end
    else
        for iCase = 1:length(dataTmp)
            for iItem = 1:length(dataTmp{1}(:))
                dataTmp{iCase}{iItem} = permute(dataTmp{iCase}{iItem}, [1 2 4 3]);
            end
        end
    end
end

% store data for all subjects
if isempty(opt.ndim)
    if strcmpi(opt.datatype, 'erp') || strcmpi(opt.datatype, 'spec')
        if length(opt.channels) > 1, dim = 3; else dim = 2; end
    else 
        if length(opt.channels) > 1, dim = 4; else dim = 3; end
    end
else
    dim = opt.ndim;
end

% check that all ERPimages have the same number of lines
if strcmpi(opt.datatype, 'erpim')
    [dataTmp,eventsTmp] = checkdataerpimage(dataTmp,eventsTmp);
    events = reorganizedata(eventsTmp, 2);
else
    events = {};
end

if ~isempty(opt.clusters)
    % Split ICA components from the same subjects need to be made 
    % as if coming from different subjects
    dataTmp2 = {};
    dataTmpSubj2 = {};
    correspInd = [];
    realDim  = dim;
    if strcmpi(opt.singletrials, 'on'), realDim = realDim+1; end
    for iDat1 = 1:length(dataTmp)
        for iDat2 = 1:length(dataTmp{iDat1}(:))
            if isempty(dataTmp{iDat1}{iDat2})
                dataTmp{iDat1}{iDat2} = double.empty(0,0,0); % sometimes empty but all dim not 0
            end
        end
        compNumbers = cellfun(@(x)size(x, realDim), dataTmp{iDat1});
        uniqComps = unique(compNumbers);
        if length(uniqComps) > 1 
            if ~(uniqComps(1) == 0 && length(uniqComps) == 2)
                error('Cannot handle conditions with different number of components for a given subject');
            end
        end
        
        if any(any(compNumbers))
            for iDat2 = 1:length(dataTmp{iDat1}(:))
                if compNumbers(iDat2)
                    for iComps = 1:compNumbers(iDat2)
                        dataTmp2{end+1}     = cell(size(dataTmp{iDat1}));
                        dataTmpSubj2{end+1} = cell(size(dataTmp{iDat1}));
                        correspInd(end+1) = iDat1;
                        % check dimensions of components
                        if ~isempty(dataTmp{iDat1}{iDat2})
                            if strcmpi(opt.singletrials, 'on') && ...
                                    (strcmpi(tmpDataType, 'timef') || strcmpi(tmpDataType, 'erpim'))
                                dataTmp2{end}{iDat2} = dataTmp{iDat1}{iDat2}(:,:,:,iComps);
                                dataTmpSubj2{end}{iDat2} = dataTmpSubj{iDat1}{iDat2}(:,:,:,1);
                            elseif strcmpi(opt.singletrials, 'on') || ...
                                    (strcmpi(tmpDataType, 'timef') || strcmpi(tmpDataType, 'erpim'))
                                dataTmp2{end}{iDat2} = dataTmp{iDat1}{iDat2}(:,:,iComps);
                                dataTmpSubj2{end}{iDat2} = dataTmpSubj{iDat1}{iDat2}(:,:,1);
                            else                                                                                                       
                                dataTmp2{end}{iDat2} = dataTmp{iDat1}{iDat2}(:,iComps);
                                dataTmpSubj2{end}{iDat2} = dataTmpSubj{iDat1}{iDat2}(:,1);
                            end
                        end
                    end
                end
            end
        end
    end
    dataTmp = dataTmp2;
    dataTmpSubj = dataTmpSubj2;
else
    correspInd = 1:length(dataTmp); % identity for channels 
end
datavals = reorganizedata(dataTmp, dim);

% fix setinds index
if nargout > 6
    datainfo = reorganizedatastruct(dataTmpSubj);
    try
        setinds = cellfun(@(x)[x.index], datainfo, 'UniformOutput', false);
    catch
        error('Error generating dataset indices, either remove the setinds output or recompute your measures')
    end
end

% reorganize data
% ---------------
function datavals = reorganizedata(dataTmp, dim)
    nonEmptyCell = find( cellfun(@isempty, dataTmp) == 0);
    datavals = cell(size(dataTmp{nonEmptyCell(1)}));
        
    % copy data
    for iItem=1:length(dataTmp{nonEmptyCell(1)}(:)')
        if dim > 1
            numItems    = sum(cellfun(@(x)size(x{iItem},dim)*(size(x{iItem},1) > 1), dataTmp)); % the size > 1 allows to detect empty array which have a non-null last dim
        else
            numItems    = sum(cellfun(@(x)size(x{iItem},dim), dataTmp)); % the size > 1 allows to detect empty array which have a non-null last dim
        end
        ind         = find(~cellfun(@(x)isempty(x{iItem}), dataTmp)); 
        if ~isempty(ind)
            ind = ind(1);
            switch dim
                case 1, datavals{iItem} = zeros([ size(dataTmp{ind}{iItem},1) numItems], 'single'); 
                case 2, datavals{iItem} = zeros([ size(dataTmp{ind}{iItem},1) numItems], 'single'); 
                case 3, datavals{iItem} = zeros([ size(dataTmp{ind}{iItem},1) size(dataTmp{ind}{iItem},2) numItems], 'single'); 
                case 4, datavals{iItem} = zeros([ size(dataTmp{ind}{iItem},1) size(dataTmp{ind}{iItem},2) size(dataTmp{ind}{iItem},3) numItems], 'single'); 
            end
        end
    end
    for iItem=1:length(dataTmp{nonEmptyCell(1)}(:)') % conditions * group
        count = 1;
        for iCase = 1:length(dataTmp) % subjects
            if ~isempty(dataTmp{iCase}{iItem})
                if dim > 1
                    numItems = size(dataTmp{iCase}{iItem},dim) * (size(dataTmp{iCase}{iItem},1) > 1); % the size > 1 allows to detect empty array which have a non-null last dim
                else
                    numItems = size(dataTmp{iCase}{iItem},dim);
                end
                switch dim
                    case 1, datavals{iItem}(:,count:count+numItems-1) = dataTmp{iCase}{iItem};
                    case 2, datavals{iItem}(:,count:count+numItems-1) = dataTmp{iCase}{iItem}; 
                    case 3, datavals{iItem}(:,:,count:count+numItems-1) = dataTmp{iCase}{iItem};
                    case 4, datavals{iItem}(:,:,:,count:count+numItems-1) = dataTmp{iCase}{iItem};
                end
                count = count+numItems;
            end
        end
    end

% reorganize data struct
% ---------------
function [datavals,setinds] = reorganizedatastruct(dataTmp)
    nonEmptyCell = find( cellfun(@isempty, dataTmp) == 0);
    datavals = cell(size(dataTmp{nonEmptyCell(1)}));
    setinds  = cell(size(dataTmp{nonEmptyCell(1)}));
        
    % copy data
    for iItem=1:length(dataTmp{nonEmptyCell(1)}(:)') % conditions * group
        count = 1;
        for iCase = 1:length(dataTmp) % subjects
            if ~isempty(dataTmp{iCase}{iItem})
                numItems = length(dataTmp{iCase}{iItem});
                setinds{iItem}(end+1) = iCase;
                datavals{iItem}(count:count+numItems-1) = dataTmp{iCase}{iItem};
                count = count+numItems;
            end
        end
    end    
    
% check data for ERPIMAGE
% -----------------------
function [dataTmp,eventTmp] = checkdataerpimage(dataTmp, eventTmp)
    
    % check second dim for ERPimage
    allsizes = [];
    for iItem=1:length(dataTmp(:))
        allsize2 = cellfun(@(x)size(x,2), dataTmp{iItem});
        allsize2( allsize2 == 0 ) = [];
        allsizes = [ allsizes allsize2 ];
    end
    if length(unique(allsizes(:))) > 1
        disp('********* Discrepancy between the number of lines in ERP-image');
    else
        return;
    end
    commonSize = min(allsizes(:));
    
    % copy data
    for iItem=1:length(dataTmp{1}(:)')
        for iCase = 1:length(dataTmp)
            if ~isempty(dataTmp{iCase}{iItem})
                % special case for ERPimage - one line missing or one line too many
                if size(dataTmp{iCase}{iItem},2)+1 == commonSize
                    dataTmp{iCase}{iItem}(:,end+1) = dataTmp{iCase}{iItem}(:,end); % duplicate last line
                    eventTmp{iCase}{iItem}(end+1) = eventTmp{iCase}{iItem}(end); % duplicate last line
                    disp('******** ERPimage discrepancy between the number of lines detected and corrected')
                elseif size(dataTmp{iCase}{iItem},2)-1 == commonSize
                    dataTmp{iCase}{iItem}(:,end) = [];
                    eventTmp{iCase}{iItem}( end) = [];
                    disp('******** ERPimage discrepancy between the number of lines detected and corrected')
                end
                if size(dataTmp{iCase}{iItem},2) ~= commonSize
                    error('ERPimage discrepancy between the number of lines');
                end
            end
        end
    end    
    
% reorganize data 2
% -----------------
function datavals = reorganizedata2(dataTmp, eventTmp)
    for iCase = 1:length(dataTmp)
        datavals(iCase).data  = dataTmp{iCase};
        datavals(iCase).event = eventTmp{iCase};
    end
    
% call newtimef (duplicate function in std_erspplot)
% --------------
function [dataout,erspbase] = processtf(dataSubject, xvals, datatype, singletrials, g)

    % compute ITC or ERSP
    if strcmpi(datatype, 'ersp')
        P = dataSubject .* conj(dataSubject);
        dataout = newtimeftrialbaseln(P, xvals, g);
        % common baseline is removed in std_erspplot
        if strcmpi(singletrials, 'off')
            dataout = squeeze(mean(dataout, 3));
        end
    else
        dataout = dataSubject;
        if strcmpi(singletrials, 'off')
            if ~isfield(g, 'itctype'), g.itctype = 'phasecoher'; end
            if ndims(dataSubject) == 4
                dataSubject = permute(dataSubject, [4 1 2 3]);
                dataout = newtimefitc(dataSubject, g.itctype);
                dataout = permute(dataout, [2 3 1]);
            else
                dataout = newtimefitc(dataSubject, g.itctype);
            end
            dataout = abs(dataout); % required for plotting scalp topo
        end
    end
  
% call erpimage
% -------------
function [dataout, eventout] = processerpim(dataSubject, events, xvals, g)

    if isempty(dataSubject), dataout = []; eventout = []; return; end
    if ~isfield(g, 'nlines'), finallines = 10; else finallines = g.nlines; end
    if ~isfield(g, 'smoothing'), smoothing = 10; else smoothing = g.smoothing; end
    
    % remove all fields and create new parameter list
    fieldList = { 'nlines' 'smoothing' 'sorttype' 'sortwin' 'sortfield' 'channels' ...
                  'interp' 'trialinfo' 'concatenate' 'savetrials' 'recompute' 'fileout' 'events' 'rmcomps'};
    params = {};
    fieldN = fieldnames(g);
    for iField = 1:length(fieldN)
        if ~ismember(fieldN{iField}, fieldList)
            params{end+1} = fieldN{iField};
            params{end+1} = g.(fieldN{iField});
        end
    end
    
    % reverse engeeneering the number of lines for ERPimage
    if ~isempty(events)
         if all(isnan(events))
             error('Cannot sort trials for one of the dataset');
         end
         lastx  = sum(~isnan(events));
    else lastx  = size(dataSubject,2);
    end
    if lastx < finallines + floor((g.smoothing-1)/2) + 3
        error('The default number of ERPimage lines is too large for one of the dataset');
    end
    firstx = 1;
    xwidth = g.smoothing;
    %xadv   = lastx/finallines;
    nout   = finallines; %floor(((lastx-firstx+xadv+1)-xwidth)/xadv);
    nlines = (lastx-xwidth)/(nout-0.5)*i; % make it imaginary
    %nlines = ceil(lastx/((lastx-firstx+1-xwidth)/(nout-1)));
           
    if ~isempty(params) && ischar(params{1}) && strcmpi(params{1}, 'components')
        params(1:2) = [];
    end
    for iChan = 1:size(dataSubject,3)
        [dataout(:,:,iChan), eventout] = erpimage(dataSubject(:,:,iChan), events, xvals, '', smoothing, nlines, 'noplot', 'on', params{:});
    end
    dataout = permute(dataout, [2 1 3]); % for tftopo
    if ~isempty(events)
        eventout = eventout'; % needs to be a column vector
    else
        eventout = [];
    end

% get file base name: filepath and sess are cell array (in case 2 files per subject)
% ----------------------------------------------------------------------------------
function filebase = getfilename(filepath, subj, sess, fileSuffix, onlyOneSession)
if onlyOneSession
    filebase = fullfile(filepath{1}, [ subj fileSuffix ] );
    if exist(filebase, 'file')
        return;
    end
    clear filebase;
end

if isempty(sess)
    sess = { '1' };
end
for iSess = 1:length(sess)
    if isnumeric(sess{iSess})
        sesStr   = [ '0' num2str(sess{iSess}) ];
    else
        sesStr   = [ '0' sess{iSess} ];
    end
    filebase{iSess} = fullfile(filepath{iSess}, [ subj '_ses-' sesStr(end-1:end) fileSuffix ] );
end
if length(unique(filebase)) < length(filebase)
    filebase = unique(filebase); % order is not important
end
    
