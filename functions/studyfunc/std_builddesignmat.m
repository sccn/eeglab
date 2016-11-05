% std_builddesignmat() - Build the design matrix for an specific design
% specified in the structure ''design'' provided as input
%
% Usage:
%
% [tmpdmat,colLabels,catflag] = std_builddesignmat(design, trialinfo, 1)
% 
% Inputs:
%  design      - Design structure as in the STUDY
%  trialinfo   - Structure of trial information. Each field should be a
%                cell array with one element for each trial
% expanding    - Expand the design matrix
% Optional inputs:
%
% Outputs:
%   tmpdmat    - Design matrix
%   colLabels  - Labels for each column of the deisgn matrix
%   catflag    - Binary vector with dimension equal to the number of columns in the design matrix. 
%                [0] mean a continuous regressor, [1] means a categotical one.
%            
% See also: std_combtrialinfo , std_plodtmat
%   
% Authors: Ramon Martinez-Cancino
%                
% Copyright (C) 2015  Ramon Martinez-Cancino, UCSD, INC, SCCN
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

function [tmpdmat,colLabels,catflag] = std_builddesignmat(design, trialinfo, expanding)

if nargin < 3, expanding = 0; end;
ntrials = length(trialinfo);
varindx = 1:length(design.variable);

% checking if 'group' var (temporal commit to detect group vars, right now only detecting variable 'group')
groupindx = find(strcmp({design.variable.label},'group'));
if ~isempty(groupindx)
    varindx(groupindx) = [];
end

tmpdmat   = NaN(ntrials,length(varindx));
catflag   = strcmp({design.variable.vartype}, 'categorical');
colLabels = {design.variable.label};

for i = 1 : length(varindx)
     % case for cont variables
    if ~catflag(i)
         [trialindsx, eventvals] = std_gettrialsind(trialinfo,design.variable(varindx(i)).label, '');
         if ~isempty(trialindsx)
             tmpdmat(trialindsx,i) = eventvals;
         end
    else % case for cat variables
        varvaluetmp = design.variable(varindx(i)).value;

        % Expanding cells cells
        c = 1; varlist = {}; varvallength = []; facval_indx = []; dmatval = 1; jcount = 1;
        for ivar = 1:length(varvaluetmp)

            % Legth of varvaluetmp(i)
            if iscellstr(varvaluetmp(ivar))
                varvallength(ivar) = length(varvaluetmp(ivar)) ;
            elseif iscell(varvaluetmp(ivar))
                varvallength(ivar) = length(varvaluetmp{ivar}) ;
            else
                varvallength(ivar) = length(varvaluetmp{ivar}) ;
            end
            % Retreiving value and assigning index in design matrix
            if ~iscell(varvaluetmp{ivar})  && varvallength(ivar) == 1
                varlist{c} = varvaluetmp{ivar};
                varindxjoint{ivar} = c;
                facval_indx(c) = dmatval;
                c = c+1;
            else
                tmpcindx = [];
                for ival = 1: varvallength(ivar)
                    if iscellstr(varvaluetmp{ivar}(ival))
                        varlist{c} =  varvaluetmp{ivar}{ival};
                    else
                        varlist{c} =  varvaluetmp{ivar}(ival);
                    end
                    facval_indx(c) = dmatval;
                    tmpcindx       = cat(1,tmpcindx,c);
                    c = c+1;
                end
                varindxjoint{ivar} = tmpcindx;
            end
            dmatval = dmatval+1;
            %---
            for j = 1 :varvallength(ivar)
                %
                if iscellstr(varlist(varindxjoint{ivar}))
                    facval = varlist{varindxjoint{ivar}(j)};
                elseif iscell(varlist(varindxjoint{ivar}))
                    tmpval = varlist(varindxjoint{ivar});
                    if isnumeric(tmpval{j})
                        facval = tmpval{j};
                    else
                        facval = tmpval(j);
                    end
                elseif  isnumeric(varlist(varindxjoint{ivar}))
                    facval = varlist(varindxjoint{ivar}(j));
                elseif ischar(varlist(varindxjoint{ivar}))
                    %                     facval = varlist(varindxjoint{ivar}(j));
                end
                % Find indices of triasl for facval
                [trialindsx, eventvals] = std_gettrialsind(trialinfo,design.variable(varindx(i)).label, facval);

                % Populating the design matrix
                if ~isempty(trialindsx)
                    tmpdmat(trialindsx,i) = facval_indx(jcount);
                    jcount = jcount+1;
                end
            end
            %---
        end 
    end
end
% -------------------------------------------------------------------------
% expand categ var
if expanding == 1
    % count number of columns
    nCols = 0;
    for iCol = 1:size(tmpdmat,2)
        if catflag(iCol)
             nCols = nCols+length(design.variable(varindx(iCol)).value);
        else
            nCols = nCols+1;
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
            for iUnique = 1:length(design.variable(varindx(iCol)).value)
                countCol    = countCol+1;
                try tmpval = uniqueVals(iUnique); catch, tmpval = NaN; end
                trialSelect = tmpdmat(:,iCol) == tmpval;
                otherCols = prevCol+setdiff(1:length(uniqueVals), iUnique);
                % putting 0 in other columns ensures we process the NaNs 
                tmpdmatExpanded(trialSelect, countCol ) = 1;
                tmpdmatExpanded(trialSelect, otherCols) = 0;
    
                % get the label for that column
                if isstr(design.variable(varindx(iCol)).value{iUnique})
                    colLabels{countCol} = [design.variable(varindx(iCol)).label '-' design.variable(varindx(iCol)).value{iUnique}];
                elseif iscell(design.variable(varindx(iCol)).value{iUnique})
                    % Concat all vals
                    varnametmp = design.variable(varindx(iCol)).value{iUnique}{1};
                    for ivar = 2: length(design.variable(varindx(iCol)).value{iUnique})
                        varnametmp = [varnametmp '&' design.variable(varindx(iCol)).value{iUnique}{ivar}];
                    end
                    colLabels{countCol} = [design.variable(varindx(iCol)).label '-' varnametmp];
                elseif isnumeric(design.variable(varindx(iCol)).value{iUnique})
                    colLabels{countCol} = [design.variable(varindx(iCol)).label '-' int2str(design.variable(varindx(iCol)).value{iUnique})];
                end
            end;
        else
            countCol = countCol+1;
            colLabels{countCol} = design.variable(varindx(iCol)).label;
            tmpdmatExpanded(:,countCol) = tmpdmat(:,iCol);
        end;
    end;

    tmpdmat = tmpdmatExpanded;
end;

tmpdmat(:,end+1) = 1;
colLabels{numel(colLabels)+1}  = 'constant';