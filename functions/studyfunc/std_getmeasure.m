function [data labels] = std_getmeasure(STUDY, ALLEEG, datind, measure)

% INPUT datind the subject index
%       measure is daterp datersp datspec icaspec icaerp icaersp

if nargin < 3
    help std_getmeasure;
    return;
end;

typeOfMeas = fastif(strcmpi(measure(1:3),'dat'), 'chan', 'comp');
trials     = zeros(1,ALLEEG(datind).trials);
data   = [];
labels = {};
suffix = '';

if strcmpi(typeOfMeas, 'chan')
    maxChans = ALLEEG(datind).nbchan;
else
    maxChans = size(ALLEEG(datind).icaweights,1);
end;
dataDim = 2;
if strcmpi(measure(end-3:end), 'ersp'), dataDim = 3; end;
if strcmpi(measure(end-3:end), 'itc'),  dataDim = 3; end;

for des = 1:length(STUDY.design)
    for ind = 1:length(STUDY.design(des).cell)
        curDesign = STUDY.design(des).cell;
        
        for iCell = 1:length(curDesign)
            if length(curDesign(iCell).dataset) == 1 && ...
                    curDesign(iCell).dataset == datind
                
                % load data
                fileName = [ curDesign(iCell).filebase '.' measure ];
                if exist(fileName)
                    curData = load('-mat', fileName);
                    
                    % find the right suffix
                    if ~isfield(curData, [ typeOfMeas '1' suffix ])
                        suffix = '_ersp';
                        if ~isfield(curData, [ typeOfMeas '1' suffix ])
                            suffix = '_itc';
                            if ~isfield(curData, [ typeOfMeas '1' suffix ])
                                fprintf('Cound not find channels or components')
                            end;
                        end;
                    end;
                            
                    tmpData = getfield(curData, [ typeOfMeas '1' suffix ]);
                    if isfield(curData, 'labels'), labels  = curData.labels; end;

                    skipFile = 0;
                    if dataDim == 2 && (size(tmpData,1) == 1 || size(tmpData,2) == 1)
                        fprintf('Dimension issue for file %s\n', fileName);
                        skipFile = 1;
                    end;
                    if dataDim == 3 && ndims(tmpData) ~= 3
                        fprintf('Dimension issue for file %s\n', fileName);
                        skipFile = 1;
                    end;
                    
                    if ~skipFile
                        
                        % create output matrix
                        if isempty(data), data = zeros( [ maxChans size(tmpData) ])*NaN; end;
                    
                        % tag trials
                        trials(curDesign(iCell).trials{1}) = 1;

                        % scan channels
                        for iChan = 1:maxChans
                            tmpData = getfield(curData, [ typeOfMeas int2str(iChan) suffix ]);
                            if dataDim == 3, data(iChan, :, :, curDesign(iCell).trials{1}) = tmpData;
                            else data(iChan, :, curDesign(iCell).trials{1}) = tmpData;
                            end;
                        end;
                        
                    end;
                end;
            end;
            if all(trials), return; end;
        end;
    end;
end;

error('Could not find all trials for this dataset');
