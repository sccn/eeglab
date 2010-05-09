% std_selsubject() - Helper function for std_erpplot(), std_specplot() 
%                    and std_erspplot() to select specific subject when
%                    plotting channel data.
% Usage:
%  >> data = std_selsubject( data, subject, setinds, allsubjects);
%
% Inputs:
%  data  -  [cell array] mean data for each subject group and/or data
%           condition. For example, to compute mean ERPs statistics from a  
%           STUDY for epochs of 800 frames in two conditions from three  
%           groups of 12 subjects,
%           >> data = { [800x12] [800x12] [800x12];... % 3 groups, cond 1
%                       [800x12] [800x12] [800x12] };  % 3 groups, cond 2
% subject - [string] subject name
% setinds - [cell array] set indices for each of the last dimension of the
%           data cell array.
%           >> setinds = { [12] [12] [12];... % 3 groups, cond 1
%                          [12] [12] [12] };  % 3 groups, cond 2
% allsubject - [cell array] all subjects (same order as in
%              STUDY.datasetinfo)
%
% Output:
%  data       - [cell array] data array with the subject or component selected
%
% Author: Arnaud Delorme, CERCO, CNRS, 2006-
% 
% See also: std_erpplot(), std_specplot() and std_erspplot()

function [data] = std_selsubject(data, subject, setinds, allsubjects, optndims);

if nargin < 2
    help std_selsubject;
    return;
end;

optndims = max(optndims, ndims(data{1}));
if isempty(strmatch(lower(subject), lower(allsubjects)))
    error(sprintf('Cannot select subject %s in list %s', subject, vararg2str({ allsubjects })));
end;

% plot specific subject
% ---------------------
for c = 1:size(data,1)
    for g = 1:size(data,2)
        for l=length(setinds{c,g}):-1:1
            if ~strcmpi(subject, allsubjects{setinds{c,g}(l)})
                if optndims == 2
                    data{c,g}(:,l) = []; %2-D
                elseif optndims == 3
                    data{c,g}(:,:,l) = []; %3-D
                else
                    data{c,g}(:,:,:,l) = []; %3-D
                end;
            end;
        end;
    end;
end;
