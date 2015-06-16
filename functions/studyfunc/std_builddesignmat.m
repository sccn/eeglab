function [tmpdmat colLabels] = std_builddesignmat(design, trialinfo, expanding)

if nargin < 3, expanding = 0; end;
ntrials = length(trialinfo);

tmpdmat = NaN(ntrials,length(design.variable));
catflag = zeros(1,length(design.variable));
for i = 1 : length(design.variable)
    
    % case for continous variables
    catflag(i) = strcmpi(design.variable(i).vartype, 'categorical');
    if ~catflag
         varlength = 1;
    else varlength = length(design.variable(i).value);
    end
    if varlength == 0, varlength = 1; end;
    colLabels{i} = design.variable(i).label;
    
    for j = 1 : varlength
        if catflag(i)
            facval = cell2mat(design.variable(i).value(j));
            if isnumeric(facval)
                facval_indx = find(facval == cell2mat(design.variable(i).value));
            else
                facval_indx = find(strcmp(facval,design.variable(i).value));
            end
        end
        
        if catflag(i)
            if isnumeric( cell2mat(design.variable(i).value(j)))
                varval = cell2mat(design.variable(i).value(j));
            else
                varval = design.variable(i).value(j);
            end
        else
            varval = '';
        end
        [trialindsx, eventvals] = std_gettrialsind(trialinfo,design.variable(i).label, varval);
        
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
                trialSelect = tmpdmat(:,iCol) == uniqueVals(iUnique);
                otherCols = prevCol+setdiff(1:length(uniqueVals), iUnique);
                % putting 0 in other columns ensures we process the NaNs 
                tmpdmatExpanded(trialSelect, countCol ) = 1;
                tmpdmatExpanded(trialSelect, otherCols) = 0;
    
                % get the label for that column
                if isstr(design.variable(iCol).value{iUnique})
                     colLabels{countCol} = [design.variable(iCol).label '-' design.variable(iCol).value{iUnique}];
                else colLabels{countCol} = [design.variable(iCol).label '-' int2str(design.variable(iCol).value{iUnique})];
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
