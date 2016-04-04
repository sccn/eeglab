function [tmpdmat colLabels] = std_builddesignmat(design, trialinfo, expanding)

if nargin < 3, expanding = 0; end;
ntrials = length(trialinfo);
varindx = 1:length(design.variable);

% checking if 'group' var (temporal commit to detect group vars, right now only detecting variable 'group')
groupindx = find(strcmp({design.variable.label},'group'));
if ~isempty(groupindx)
    varindx(groupindx) = [];
end

tmpdmat = NaN(ntrials,length(varindx));
catflag = zeros(1,length(varindx));
for i = 1 : length(varindx)
    
    % case for continous variables
    catflag(i) = strcmpi(design.variable(varindx(i)).vartype, 'categorical');
    if ~catflag
         varlength = 1;
    else
        varvaluetmp = design.variable(varindx(i)).value;
        cellindx    = find(cellfun(@iscell,varvaluetmp));
        c = 1;
        for ivar = 1:length(varvaluetmp)
                if ~iscell(varvaluetmp{ivar})
                    varlist{c} = varvaluetmp{ivar};
                    c = c+1;
                else
                    for ival = 1: length(varvaluetmp)
                        varlist(c) =  varvaluetmp{ivar}(ival);
                        c = c+1;
                    end
                end
        end        
        varlength = length(varlist);

    end
    if varlength == 0, varlength = 1; end;
    colLabels{i} = design.variable(varindx(i)).label;
    
    for j = 1 : varlength
        if catflag(i)
            
            %facval = cell2mat(design.variable(varindx(i)).value(j));
            facval = varlist{j};
            if isnumeric(facval)
                if isempty(cellindx)
                    facval_indx = find(facval == cell2mat(design.variable(varindx(i)).value));
                else
                end
            else
                if isempty(cellindx)
                    facval_indx = find(strcmp(facval,design.variable(varindx(i)).value));
                else
                    for ivar = 1:length(varvaluetmp)
                        hittmp = find(strcmp(facval,varvaluetmp{ivar}));
                        if ~isempty(hittmp)
                            facval_indx = ivar;
                        end
                    end
                end
            end
           
            %
            if isnumeric( cell2mat(varlist(j)))
            %if isnumeric( cell2mat(design.variable(varindx(i)).value(j)))
                varval = cell2mat(design.variable(varindx(i)).value(j));
            else
                varval = varlist(j);
            end
        else
            varval = '';
        end
        [trialindsx, eventvals] = std_gettrialsind(trialinfo,design.variable(varindx(i)).label, varval);
        
        if ~isempty(trialindsx)
            % case for continous variables
            if ~catflag(i)
                facval_indx = eventvals;
            end
            tmpdmat(trialindsx,i) = facval_indx;
            
        end
    end
end

% expand categ var
if expanding == 1
    % count number of columns
    nCols = 0;
    for iCol = 1:size(tmpdmat,2)
        if catflag(iCol)
             nCols = nCols+length(design.variable(iCol).value);
        else nCols = nCols+1;
        end;
    end;
    
    tmpdmatExpanded = NaN(size(tmpdmat,1),nCols);
    countCol = 0;
    for iCol = 1:size(tmpdmat,2)
        if catflag(iCol)
            % get unique values for this given categ var
            uniqueVals = unique(tmpdmat(:,iCol));
            uniqueVals(isnan(uniqueVals)) = [];
            
            % scan unique values
            prevCol = countCol;
            for iUnique = 1:length(design.variable(iCol).value)
                countCol    = countCol+1;
                try tmpval = uniqueVals(iUnique); catch, tmpval = NaN; end
                trialSelect = tmpdmat(:,iCol) == tmpval;
                otherCols = prevCol+setdiff(1:length(uniqueVals), iUnique);
                % putting 0 in other columns ensures we process the NaNs 
                tmpdmatExpanded(trialSelect, countCol ) = 1;
                tmpdmatExpanded(trialSelect, otherCols) = 0;
    
                % get the label for that column
                if isstr(design.variable(iCol).value{iUnique})
                    colLabels{countCol} = [design.variable(iCol).label '-' design.variable(iCol).value{iUnique}];
                elseif iscell(design.variable(iCol).value{iUnique})
                    % Concat all vals
                    varnametmp = design.variable(iCol).value{iUnique}{1};
                    for ivar = 2: length(design.variable(iCol).value{iUnique})
                        varnametmp = [varnametmp '&' design.variable(iCol).value{iUnique}{ivar}];
                    end
                    colLabels{countCol} = [design.variable(iCol).label '-' varnametmp];
                elseif isnumeric(design.variable(iCol).value{iUnique})
                    colLabels{countCol} = [design.variable(iCol).label '-' int2str(design.variable(iCol).value{iUnique})];
                end
            end;
        else
            countCol = countCol+1;
            colLabels{countCol} = design.variable(iCol).label;
            tmpdmatExpanded(:,countCol) = (tmpdmat(:,iCol)-min(tmpdmat(:,iCol)))/(max(tmpdmat(:,iCol))-min(tmpdmat(:,iCol)));
        end;
    end;

    tmpdmat = tmpdmatExpanded;
end;

tmpdmat(:,end+1) = 1;
colLabels{numel(colLabels)+1}  = 'constant';
