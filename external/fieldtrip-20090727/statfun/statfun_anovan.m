function [s] = statfun_xxx(cfg, dat, design);

% STATFUN_ANOVA is a function for computing ANOVA for the relation
% between biological data and a design vector containing trial
% classifications or another independent variable. This only works for
% repeated measures ANOVA (unpaired data repetition).
%
% The external interface of this function has to be
%   [s] = statfun_anova(cfg, dat, design);
% where
%   dat    contains the biological data, Nvoxels x Nreplications
%   design contains the independent variable,  1 x Nreplications
%
% Author: Arnaud Delorme, 2009 using the anova_rm function from Arash Salarian

cond1 = unique(design(1,:));
cond2 = unique(design(2,:));
alldat = {};
for ind1=1:length(cond1)
    for ind2=1:length(cond2)
        finds{ind1,ind2} = find(design(2,:) == cond2(ind2) &  design(1,:) == cond1(ind1));
    end;
end;

s.stat = zeros(size(dat,1),1);
s.prob = zeros(size(dat,1),1);
for ntest = 1:size(dat,1)
    
    % get the data
    alldat = {};
    finaldat = {};
    for ind1=1:length(cond1)
        for ind2=1:length(cond2)
            alldat{ind1, ind2} = dat(ntest, finds{ind1,ind2})';
        end;
        finaldat{ind1} = [ alldat{ind1, :} ];
    end;

    [p stats ] = anova_rm(finaldat, 'off');
    s.stat(ntest) = stats{4,5};
    s.prob(ntest) = stats{4,6};    
end;

%s.crival = finv(0.05, length(cond1)+length(cond2)-2);

function [p, table] = anova_rm(X, displayopt)
%   [p, table] = anova_rm(X, displayopt)
%   Single factor, tepeated measures ANOVA.
%
%   [p, table] = anova_rm(X, displayopt) performs a repeated measures ANOVA
%   for comparing the means of two or more columns (time) in one or more
%   samples(groups). Unbalanced samples (i.e. different number of subjects 
%   per group) is supported though the number of columns (followups)should 
%   be the same. 
%
%   DISPLAYOPT can be 'on' (the default) to display the ANOVA table, or 
%   'off' to skip the display. For a design with only one group and two or 
%   more follow-ups, X should be a matrix with one row for each subject. 
%   In a design with multiple groups, X should be a cell array of matrixes.
% 
%   Example: Gait-Cycle-times of a group of 7 PD patients have been
%   measured 3 times, in one baseline and two follow-ups:
%
%   patients = [
%    1.1015    1.0675    1.1264
%    0.9850    1.0061    1.0230
%    1.2253    1.2021    1.1248
%    1.0231    1.0573    1.0529
%    1.0612    1.0055    1.0600
%    1.0389    1.0219    1.0793
%    1.0869    1.1619    1.0827 ];
%
%   more over, a group of 8 controls has been measured in the same protocol:
%
%   controls = [
%     0.9646    0.9821    0.9709
%     0.9768    0.9735    0.9576
%     1.0140    0.9689    0.9328
%     0.9391    0.9532    0.9237
%     1.0207    1.0306    0.9482
%     0.9684    0.9398    0.9501
%     1.0692    1.0601    1.0766
%     1.0187    1.0534    1.0802 ];
%
%   We are interested to ssfdsgdee if the performance of the patients for the
%   followups were the same or not:
%  
%   p = anova_rm(patients);
%
%   By considering the both groups, we can also check to see if the 
%   follow-ups were significantly different and also check two see that the
%   two groups had a different performance:
%
%   p = anova_rm({patients controls});
%
%
%   ref: Statistical Methods for the Analysis of Repeated Measurements, 
%     C. S. Daivs, Springer, 2002
%
%   Copyright 2008, Arash Salarian
%   mailto://arash.salarian@ieee.org
%

if nargin < 2
    displayopt = 'on';
end

if ~iscell(X)
    X = {X};
end

%number of groups
s = size(X,2);  

%subjects per group 
n_h = zeros(s, 1);
for h=1:s
    n_h(h) = size(X{h}, 1);    
end
n = sum(n_h);

%number of follow-ups
t = size(X{1},2);   

% overall mean
y = 0;
for h=1:s
    y = y + sum(sum(X{h}));
end
y = y / (n * t);

% allocate means
y_h = zeros(s,1);
y_j = zeros(t,1);
y_hj = zeros(s,t);
y_hi = cell(s,1);
for h=1:s
    y_hi{h} = zeros(n_h(h),1);
end

% group means
for h=1:s
    y_h(h) = sum(sum(X{h})) / (n_h(h) * t);
end

% follow-up means
for j=1:t
    y_j(j) = 0;
    for h=1:s
        y_j(j) = y_j(j) + sum(X{h}(:,j));
    end
    y_j(j) = y_j(j) / n;
end

% group h and time j mean
for h=1:s
    for j=1:t
        y_hj(h,j) = sum(X{h}(:,j) / n_h(h));
    end
end

% subject i'th of group h mean
for h=1:s
    for i=1:n_h(h)
        y_hi{h}(i) = sum(X{h}(i,:)) / t;
    end
end

% calculate the sum of squares
ssG = 0;
ssSG = 0;
ssT = 0;
ssGT = 0;
ssR = 0;

for h=1:s
    for i=1:n_h(h)
        for j=1:t
            ssG  = ssG  + (y_h(h) - y)^2;
            ssSG = ssSG + (y_hi{h}(i) - y_h(h))^2;
            ssT  = ssT  + (y_j(j) - y)^2;
            ssGT = ssGT + (y_hj(h,j) - y_h(h) - y_j(j) + y)^2;
            ssR  = ssR  + (X{h}(i,j) - y_hj(h,j) - y_hi{h}(i) + y_h(h))^2;
        end
    end
end

% calculate means
if s > 1
    msG  = ssG  / (s-1);
    msGT = ssGT / ((s-1)*(t-1));
end
msSG = ssSG / (n-s);
msT  = ssT  / (t-1);
msR  = ssR  / ((n-s)*(t-1));


% calculate the F-statistics
if s > 1
    FG  = msG  / msSG;
    FGT = msGT / msR;
end
FT  = msT  / msR;
FSG = msSG / msR;


% single or multiple sample designs?
if s > 1
    % case for multiple samples
    pG  = 1 - fcdf(FG, s-1, n-s);
    pT  = 1 - fcdf(FT, t-1, (n-s)*(t-1));
    pGT = 1 - fcdf(FGT, (s-1)*(t-1), (n-s)*(t-1));
    pSG = 1 - fcdf(FSG, n-s, (n-s)*(t-1));

    p = [pT, pG, pSG, pGT];

    table = { 'Source' 'SS' 'df' 'MS' 'F' 'Prob>F'
        'Time'  ssT t-1 msT FT pT
        'Group' ssG s-1 msG FG pG
        'Ineratcion' ssGT (s-1)*(t-1) msGT FGT pGT
        'Subjects (matching)' ssSG n-s msSG FSG pSG
        'Error' ssR (n-s)*(t-1) msR  [] []
        'Total' [] [] [] [] []
        };
    table{end, 2} = sum([table{2:end-1,2}]);
    table{end, 3} = sum([table{2:end-1,3}]);

    if (isequal(displayopt, 'on'))
        digits = [-1 -1 0 -1 2 4];
        statdisptable(table, 'multi-sample repeated measures ANOVA', 'ANOVA Table', '', digits);
    end
else
    % case for only one sample
    pT  = 1 - fcdf(FT, t-1, (n-s)*(t-1));
    pSG = 1 - fcdf(FSG, n-s, (n-s)*(t-1));

    p = [pT, pSG];

    table = { 'Source' 'SS' 'df' 'MS' 'F' 'Prob>F'
        'Time'  ssT t-1 msT FT pT
        'Subjects (matching)' ssSG n-s msSG FSG pSG
        'Error' ssR (n-s)*(t-1) msR  [] []
        'Total' [] [] [] [] []
        };
    table{end, 2} = sum([table{2:end-1,2}]);
    table{end, 3} = sum([table{2:end-1,3}]);

    if (isequal(displayopt, 'on'))
        digits = [-1 -1 0 -1 2 4];
        statdisptable(table, 'repeated measures ANOVA', 'ANOVA Table', '', digits);
    end
end
