% teststat - EEGLAB statistical testing function 
%
% Statistics are critical for inference testing in Science. It is thus
% primordial to make sure than all the statistics implemented are 
% robust and at least bug free. Statistical function using complex
% formulas are inherently prone to bugs. EEGLAB functions are all the 
% more prone to bugs given that they only use complex Matlab code to
% avoid loops and speed up computation.
%
% This test function does not garantee that EEGLAB statistical functions
% are bug free. It does assure though that bugs are unlikely and minor
% if they are present.
%
% This function test 3 things. 
%
% * First, it checks that for vector inputs the EEGLAB functions return 
%   the same output as other reference functions from the Matlab statistical 
%   toolbox or from other packages tested against the SPSS software for 
%   repeated measure ANOVA (rm_anova2 function).
%
% * Second, it checks that array inputs with different number of dimensions
%   (from 1 to 3) the EEGLAB function return the same output.
%
% * Third, it checks that the permutation and bootstrap methods shuffle
%   the data properly by running multiple tests.

% Copyright (C) 2006 Arnaud Delorme
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

function teststat;

% testing paired t-test
% ---------------------
a = { rand(1,10) rand(1,10)+0.5 };
[t df pvals surog] = statcond(a, 'mode', 'param', 'verbose', 'off', 'paired', 'on');
[h p tmp stats] = ttest(a{1}, a{2});
fprintf('Statistics paired statcond    t-value %2.2f df=%d p=%0.4f\n', t, df, pvals); 
fprintf('Statistics paired ttest func. t-value %2.2f df=%d p=%0.4f\n', stats.tstat, stats.df, p); 
assertsame([t stats.tstat], [df stats.df], [pvals p]); 
disp('--------------------');

% testing unpaired t-test
% -----------------------
[t df pvals surog] = statcond(a, 'mode', 'param', 'verbose', 'off', 'paired', 'off');
[h p tmp stats] = ttest2(a{1}, a{2});
fprintf('Statistics paired statcond     t-value %2.2f df=%d p=%0.4f\n', t, df, pvals); 
fprintf('Statistics paired ttest2 func. t-value %2.2f df=%d p=%0.4f\n', stats.tstat, stats.df, p); 
assertsame([t stats.tstat], [df stats.df], [pvals p]); 
disp('--------------------');

% testing paired 1-way ANOVA
% --------------------------
a = { rand(1,10) rand(1,10) rand(1,10)+0.2; rand(1,10) rand(1,10)+0.2 rand(1,10) };
[F df pvals surog] = statcond(a(1,:), 'mode', 'param', 'verbose', 'off', 'paired', 'on');
z = zeros(10,1); o = ones(10,1); t = ones(10,1)*2;
stats = rm_anova2(  [ a{1,1}';a{1,2}';a{1,3}'], repmat([1:10]', [3 1]), [o;o;o], [z;o;t], {'a','b'});
fprintf('Statistics 1-way paired statcond        F-value %2.2f df1=%d df2=%d p=%0.4f\n', F, df(1), df(2), pvals); 
fprintf('Statistics 1-way paired rm_avova2 func. F-value %2.2f df1=%d df2=%d p=%0.4f\n', stats{3,5}, stats{3,3}, stats{6,3}, stats{3,6}); 
assertsame([F stats{3,5}], [df(1) stats{3,3}], [df(2) stats{6,3}], [pvals stats{3,6}]); 
disp('--------------------');

% testing paired 2-way ANOVA
% --------------------------
[F df pvals surog] = statcond(a, 'mode', 'param', 'verbose', 'off', 'paired', 'on');
z = zeros(10,1); o = ones(10,1); t = ones(10,1)*2;
stats = rm_anova2(  [ a{1,1}';a{1,2}';a{1,3}';a{2,1}';a{2,2}';a{2,3}' ], ...
            repmat([1:10]', [6 1]), [o;o;o;z;z;z], [z;o;t;z;o;t], {'a','b'});
fprintf('Statistics 2-way paired statcond        F-value %2.2f df1=%d df2=%d p=%0.4f\n', F{3}, df{3}(1), df{3}(2), pvals{3}); 
fprintf('Statistics 2-way paired rm_avova2 func. F-value %2.2f df1=%d df2=%d p=%0.4f\n', stats{4,5}, stats{4,3}, stats{7,3}, stats{4,6}); 
assertsame([F{3} stats{4,5}], [df{3}(1) stats{4,3}], [df{3}(2) stats{7,3}], [pvals{3} stats{4,6}]); 
disp('--------------------');

% testing 1-way unpaired ANOVA
% ----------------------------
[F df pvals surog] = statcond(a(1,:), 'mode', 'param', 'verbose', 'off', 'paired', 'off');
[p stats] = anova1( [ a{1,1}' a{1,2}' a{1,3}' ],{}, 'off');
fprintf('Statistics 1-way unpaired statcond     F-value %2.2f df1=%d df2=%d p=%0.4f\n', F, df(1), df(2), pvals); 
fprintf('Statistics 1-way unpaired anova1 func. F-value %2.2f df1=%d df2=%d p=%0.4f\n', stats{2,5}, stats{2,3}, stats{3,3}, stats{2,6}); 
assertsame([F stats{2,5}], [df(1) stats{2,3}], [df(2) stats{3,3}], [pvals stats{2,6}]); 
disp('--------------------');

% testing 2-way unpaired ANOVA
% ----------------------------
[F df pvals surog] = statcond(a, 'mode', 'param', 'verbose', 'off', 'paired', 'off');
[p stats] = anova2( [ a{1,1}' a{1,2}' a{1,3}'; a{2,1}' a{2,2}' a{2,3}' ], 10, 'off');
fprintf('Statistics 2-way paired statcond       F-value %2.2f df1=%d df2=%d p=%0.4f\n', F{3}, df{3}(1), df{3}(2), pvals{3}); 
fprintf('Statistics 1-way unpaired anova2 func. F-value %2.2f df1=%d df2=%d p=%0.4f\n', stats{4,5}, stats{4,3}, stats{5,3}, stats{4,6}); 
assertsame([F{3} stats{4,5}], [df{3}(1) stats{4,3}], [df{3}(2) stats{5,3}], [pvals{3} stats{4,6}]); 
disp('--------------------');

% testing different dimensions in statcond
% ----------------------------------------
a = { rand(1,10)      rand(1,10)+0.5      rand(1,10)};
b = { rand(10,10)     rand(10,10)+0.5     rand(10,10)};     b{1}(4,:)     = a{1}; b{2}(4,:)     = a{2}; b{3}(4,:)     = a{3};
c = { rand(5,10,10)   rand(5,10,10)+0.5   rand(5,10,10)};   c{1}(2,4,:)   = a{1}; c{2}(2,4,:)   = a{2}; c{3}(2,4,:)   = a{3};
d = { rand(2,5,10,10) rand(2,5,10,10)+0.5 rand(2,5,10,10)}; d{1}(1,2,4,:) = a{1}; d{2}(1,2,4,:) = a{2}; d{3}(1,2,4,:) = a{3};
[t1 df1 pvals1] = statcond(a(1:2), 'mode', 'param', 'verbose', 'off', 'paired', 'on');
[t2 df2 pvals2] = statcond(b(1:2), 'mode', 'param', 'verbose', 'off', 'paired', 'on');
[t3 df3 pvals3] = statcond(c(1:2), 'mode', 'param', 'verbose', 'off', 'paired', 'on');
[t4 df4 pvals4] = statcond(d(1:2), 'mode', 'param', 'verbose', 'off', 'paired', 'on');
fprintf('Statistics paired statcond t-test dim1 t-value %2.2f df=%d p=%0.4f\n', t1, df1, pvals1);
fprintf('Statistics paired statcond t-test dim2 t-value %2.2f df=%d p=%0.4f\n', t2(4), df2, pvals2(4));
fprintf('Statistics paired statcond t-test dim3 t-value %2.2f df=%d p=%0.4f\n', t3(2,4), df3, pvals3(2,4));
fprintf('Statistics paired statcond t-test dim4 t-value %2.2f df=%d p=%0.4f\n', t4(1,2,4), df4, pvals4(1,2,4));
assertsame([t1 t2(4) t3(2,4) t4(1,2,4)], [df1 df2 df3 df4], [pvals1 pvals2(4) pvals3(2,4) pvals4(1,2,4)]); 
disp('--------------------');
[t1 df1 pvals1] = statcond(a(1:2), 'mode', 'param', 'verbose', 'off', 'paired', 'off');
[t2 df2 pvals2] = statcond(b(1:2), 'mode', 'param', 'verbose', 'off', 'paired', 'off');
[t3 df3 pvals3] = statcond(c(1:2), 'mode', 'param', 'verbose', 'off', 'paired', 'off');
[t4 df4 pvals4] = statcond(d(1:2), 'mode', 'param', 'verbose', 'off', 'paired', 'off');
fprintf('Statistics unpaired statcond t-test dim1 t-value %2.2f df=%d p=%0.4f\n', t1, df1, pvals1);
fprintf('Statistics unpaired statcond t-test dim2 t-value %2.2f df=%d p=%0.4f\n', t2(4), df2, pvals2(4));
fprintf('Statistics unpaired statcond t-test dim3 t-value %2.2f df=%d p=%0.4f\n', t3(2,4), df3, pvals3(2,4));
fprintf('Statistics unpaired statcond t-test dim4 t-value %2.2f df=%d p=%0.4f\n', t4(1,2,4), df4, pvals4(1,2,4));
assertsame([t1 t2(4) t3(2,4) t4(1,2,4)], [df1 df2 df3 df4], [pvals1 pvals2(4) pvals3(2,4) pvals4(1,2,4)]); 
disp('--------------------');
[t1 df1 pvals1] = statcond(a, 'mode', 'param', 'verbose', 'off', 'paired', 'on');
[t2 df2 pvals2] = statcond(b, 'mode', 'param', 'verbose', 'off', 'paired', 'on');
[t3 df3 pvals3] = statcond(c, 'mode', 'param', 'verbose', 'off', 'paired', 'on');
[t4 df4 pvals4] = statcond(d, 'mode', 'param', 'verbose', 'off', 'paired', 'on');
fprintf('Statistics paired statcond anova 1-way dim1 t-value %2.2f df1=%d df2=%d p=%0.4f\n', t1, df1(1), df1(2), pvals1);
fprintf('Statistics paired statcond anova 1-way dim2 t-value %2.2f df1=%d df2=%d p=%0.4f\n', t2(4), df2(1), df2(2), pvals2(4));
fprintf('Statistics paired statcond anova 1-way dim3 t-value %2.2f df1=%d df2=%d p=%0.4f\n', t3(2,4), df3(1), df3(2), pvals3(2,4));
fprintf('Statistics paired statcond anova 1-way dim4 t-value %2.2f df1=%d df2=%d p=%0.4f\n', t4(1,2,4), df4(1), df4(2), pvals4(1,2,4));
assertsame([t1 t2(4) t3(2,4) t4(1,2,4)], [df1 df2 df3 df4], [pvals1 pvals2(4) pvals3(2,4) pvals4(1,2,4)]); 
disp('--------------------');
[t1 df1 pvals1] = statcond(a, 'mode', 'param', 'verbose', 'off', 'paired', 'off');
[t2 df2 pvals2] = statcond(b, 'mode', 'param', 'verbose', 'off', 'paired', 'off');
[t3 df3 pvals3] = statcond(c, 'mode', 'param', 'verbose', 'off', 'paired', 'off');
[t4 df4 pvals4] = statcond(d, 'mode', 'param', 'verbose', 'off', 'paired', 'off');
fprintf('Statistics unpaired statcond anova 1-way dim1 t-value %2.2f df1=%d df2=%d p=%0.4f\n', t1, df1(1), df1(2), pvals1);
fprintf('Statistics unpaired statcond anova 1-way dim2 t-value %2.2f df1=%d df2=%d p=%0.4f\n', t2(4), df2(1), df2(2), pvals2(4));
fprintf('Statistics unpaired statcond anova 1-way dim3 t-value %2.2f df1=%d df2=%d p=%0.4f\n', t3(2,4), df3(1), df3(2), pvals3(2,4));
fprintf('Statistics unpaired statcond anova 1-way dim4 t-value %2.2f df1=%d df2=%d p=%0.4f\n', t4(1,2,4), df4(1), df4(2), pvals4(1,2,4));
assertsame([t1 t2(4) t3(2,4) t4(1,2,4)], [df1 df2 df3 df4], [pvals1 pvals2(4) pvals3(2,4) pvals4(1,2,4)]); 
disp('--------------------');
a(2,:) = a; a{1} = a{1}/2;
b(2,:) = b; b{1} = b{1}/2;
c(2,:) = c; c{1} = c{1}/2;
d(2,:) = d; d{1} = d{1}/2;
[t1 df1 pvals1] = statcond(a, 'mode', 'param', 'verbose', 'off', 'paired', 'on');
[t2 df2 pvals2] = statcond(b, 'mode', 'param', 'verbose', 'off', 'paired', 'on');
[t3 df3 pvals3] = statcond(c, 'mode', 'param', 'verbose', 'off', 'paired', 'on');
[t4 df4 pvals4] = statcond(d, 'mode', 'param', 'verbose', 'off', 'paired', 'on');
fprintf('Statistics paired statcond anova 2-way dim1 t-value %2.2f df1=%d df2=%d p=%0.4f\n', t1{3}, df1{3}(1), df1{3}(2), pvals1{3});
fprintf('Statistics paired statcond anova 2-way dim2 t-value %2.2f df1=%d df2=%d p=%0.4f\n', t2{3}(4), df2{3}(1), df2{3}(2), pvals2{3}(4));
fprintf('Statistics paired statcond anova 2-way dim3 t-value %2.2f df1=%d df2=%d p=%0.4f\n', t3{3}(2,4), df3{3}(1), df3{3}(2), pvals3{3}(2,4));
fprintf('Statistics paired statcond anova 2-way dim4 t-value %2.2f df1=%d df2=%d p=%0.4f\n', t4{3}(1,2,4), df4{3}(1), df4{3}(2), pvals4{3}(1,2,4));
assertsame([t1{3} t2{3}(4) t3{3}(2,4) t4{3}(1,2,4)], [df1{3}(1) df2{3}(1) df3{3}(1) df4{3}(1)], [df1{3}(2) df2{3}(2) df3{3}(2) df4{3}(2)], [pvals1{3} pvals2{3}(4) pvals3{3}(2,4) pvals4{3}(1,2,4)]); 
disp('--------------------');
[t1 df1 pvals1] = statcond(a, 'mode', 'param', 'verbose', 'off', 'paired', 'off');
[t2 df2 pvals2] = statcond(b, 'mode', 'param', 'verbose', 'off', 'paired', 'off');
[t3 df3 pvals3] = statcond(c, 'mode', 'param', 'verbose', 'off', 'paired', 'off');
[t4 df4 pvals4] = statcond(d, 'mode', 'param', 'verbose', 'off', 'paired', 'off');
fprintf('Statistics unpaired statcond anova 2-way dim1 t-value %2.2f df1=%d df2=%d p=%0.4f\n', t1{3}, df1{3}(1), df1{3}(2), pvals1{3});
fprintf('Statistics unpaired statcond anova 2-way dim2 t-value %2.2f df1=%d df2=%d p=%0.4f\n', t2{3}(4), df2{3}(1), df2{3}(2), pvals2{3}(4));
fprintf('Statistics unpaired statcond anova 2-way dim3 t-value %2.2f df1=%d df2=%d p=%0.4f\n', t3{3}(2,4), df3{3}(1), df3{3}(2), pvals3{3}(2,4));
fprintf('Statistics unpaired statcond anova 2-way dim4 t-value %2.2f df1=%d df2=%d p=%0.4f\n', t4{3}(1,2,4), df4{3}(1), df4{3}(2), pvals4{3}(1,2,4));
assertsame([t1{3} t2{3}(4) t3{3}(2,4) t4{3}(1,2,4)], [df1{3}(1) df2{3}(1) df3{3}(1) df4{3}(1)], [df1{3}(2) df2{3}(2) df3{3}(2) df4{3}(2)], [pvals1{3} pvals2{3}(4) pvals3{3}(2,4) pvals4{3}(1,2,4)]); 
disp('--------------------');

% testing shuffling and permutation for bootstrap
% -----------------------------------------------
clear a;
m1 = [1:10];
m2 = [1:10]+100;
m3 = [1:10]+1000;
a{1} = { m1 m2 };
a{2} = { m1 m2 m3 };
a{3} = { [ zeros(9,10); m1] [ zeros(9,10); m2] };
a{4} = { [ zeros(9,10); m1] [ zeros(9,10); m2] [ zeros(9,10); m3] };
tmpa = zeros(9,8,10); tmpa(end,end,:) = m1;
tmpb = zeros(9,8,10); tmpb(end,end,:) = m2;
tmpc = zeros(9,8,10); tmpc(end,end,:) = m3;
a{5} = { tmpa tmpb };
a{6} = { tmpa tmpb tmpc };

for method = 1:2
    if method == 2, opt = {'arraycomp', 'off'}; else opt = {}; end
    for dim = 1:length(a)
        [sa1] = statcond(a{dim}, 'mode', 'bootstrap', 'verbose', 'off', 'paired', 'on',  'returnresamplingarray', 'on', opt{:}, 'naccu', 10);
        [sa2] = statcond(a{dim}, 'mode', 'perm'     , 'verbose', 'off', 'paired', 'on',  'returnresamplingarray', 'on', opt{:}, 'naccu', 10);
        [sa3] = statcond(a{dim}, 'mode', 'bootstrap', 'verbose', 'off', 'paired', 'off', 'returnresamplingarray', 'on', opt{:}, 'naccu', 10);
        [sa4] = statcond(a{dim}, 'mode', 'perm'     , 'verbose', 'off', 'paired', 'off', 'returnresamplingarray', 'on', opt{:}, 'naccu', 10);
        
        % select data
        nd = ndims(sa1{1});
        if nd == 2 && size(sa1{1},2) > 1
            for t=1:length(sa1), 
                sa1{t} = sa1{t}(end,:); 
                sa2{t} = sa2{t}(end,:); 
                sa3{t} = sa3{t}(end,:); 
                sa4{t} = sa4{t}(end,:); 
            end
        elseif nd == 3
            for t=1:length(sa1), 
                sa1{t} = squeeze(sa1{t}(end,end,:));
                sa2{t} = squeeze(sa2{t}(end,end,:));
                sa3{t} = squeeze(sa3{t}(end,end,:));
                sa4{t} = squeeze(sa4{t}(end,end,:));
            end
        elseif nd == 4
            for t=1:length(sa1), 
                sa1{t} = squeeze(sa1{t}(end,end,end,:)); 
                sa2{t} = squeeze(sa2{t}(end,end,end,:));
                sa3{t} = squeeze(sa3{t}(end,end,end,:));
                sa4{t} = squeeze(sa4{t}(end,end,end,:));
            end
        end
        
        % for paired bootstrap, we make sure that the resampling has only shuffled between conditions
        % for instance [101 2 1003 104 ...] is an acceptable sequence 
        if all(rem(sa1{1}(:)',10) == [1:9 0]) && all(rem(sa1{2}(:)',10) == [1:9 0])
            fprintf('Bootstrap paired dim%d resampling method %d Pass\n', dim, method);
        else  error('Bootstrap paired resampling Error');
        end
        % for paired permutation, in addition, we make sure that the sum accross condition is constant
        % which is not true for bootstrap
        msa = meansa(sa2); msa = msa(:)-msa(1);
        if all(rem(sa1{1}(:)',10) == [1:9 0]) && all(rem(sa1{2}(:)',10) == [1:9 0]) && ...
            all(round(msa) == [0:9]') && length(unique(sa2{1})) == 10 && length(unique(sa2{2})) == 10
            fprintf('Permutation paired dim%d resampling method %d Pass\n', dim, method);
        else  error('Permutation paired resampling Error');
        end
        % for unpaired bootstrap, only make sure there are enough unique
        % values
        if length(unique(sa3{1})) > 3 && length(unique(sa3{2})) > 3
             fprintf('Bootstrap unpaired dim%d reampling method %d Pass\n', dim, method);
        else   error('Bootstrap unpaired reampling Error');
        end
        % for unpaired permutation, the number of unique values must be 10
        % and the sum must be constant (not true for bootstrap)
        if length(unique(sa4{1})) == 10 && length(unique(sa4{2})) == 10 && ( floor(mean(meansa(sa4))) == 55 || floor(mean(meansa(sa4))) == 372 )
             fprintf('Permutation unpaired dim%d reampling method %d Pass\n', dim, method);
        else   error('Permutation unpaired reampling Error');
        end
       
        disp('------------------------');
    end
end

% function to check 
function assertsame(varargin)

for ind = 1:length(varargin)
    if length(varargin{1}) > 2
        for tmpi = 1:length(varargin{1})-1
            assertsame(varargin{1}(tmpi:tmpi+1));
        end
        return;
    else
        if (varargin{ind}(1)-varargin{ind}(2)) > abs(mean(varargin{ind}))*0.01
            error('Test failed');
        end
    end
end
disp('Test pass');

function [meanmat] = meansa(mat)

meanmat = zeros(size(mat{1}));
for index = 1:length(mat)
    meanmat = meanmat+mat{index}/length(mat);
end

function stats = rm_anova2(Y,S,F1,F2,FACTNAMES)
%
% function stats = rm_anova2(Y,S,F1,F2,FACTNAMES)
%
% Two-factor, within-subject repeated measures ANOVA.
% For designs with two within-subject factors.
%
% Parameters:
%    Y          dependent variable (numeric) in a column vector
%    S          grouping variable for SUBJECT
%    F1         grouping variable for factor #1
%    F2         grouping variable for factor #2
%    F1name     name (character array) of factor #1
%    F2name     name (character array) of factor #2
%
%    Y should be a 1-d column vector with all of your data (numeric).
%    The grouping variables should also be 1-d numeric, each with same
%    length as Y. Each entry in each of the grouping vectors indicates the
%    level # (or subject #) of the corresponding entry in Y.
%
% Returns:
%    stats is a cell array with the usual ANOVA table:
%      Source / ss / df / ms / F / p
%
% Notes:
%    Program does not do any input validation, so it is up to you to make
%    sure that you have passed in the parameters in the correct form:
%
%       Y, S, F1, and F2 must be numeric vectors all of the same length.
%
%       There must be at least one value in Y for each possible combination
%       of S, F1, and F2 (i.e. there must be at least one measurement per
%       subject per condition).
%
%       If there is more than one measurement per subject X condition, then
%       the program will take the mean of those measurements.
%
% Aaron Schurger (2005.02.04)
%   Derived from Keppel & Wickens (2004) "Design and Analysis" ch. 18
%

%
% Revision history...
%
% 11 December 2009 (Aaron Schurger)
% 
% Fixed error under "bracket terms"
% was: expY = sum(Y.^2);
% now: expY = sum(sum(sum(MEANS.^2)));
%

stats = cell(4,5);

F1_lvls = unique_bc(F1);
F2_lvls = unique_bc(F2);
Subjs = unique_bc(S);

a = length(F1_lvls); % # of levels in factor 1
b = length(F2_lvls); % # of levels in factor 2
n = length(Subjs); % # of subjects

INDS = cell(a,b,n); % this will hold arrays of indices
CELLS = cell(a,b,n); % this will hold the data for each subject X condition
MEANS = zeros(a,b,n); % this will hold the means for each subj X condition

% Calculate means for each subject X condition.
% Keep data in CELLS, because in future we may want to allow options for
% how to compute the means (e.g. leaving out outliers > 3stdev, etc...).
for i=1:a % F1
    for j=1:b % F2
        for k=1:n % Subjs
            INDS{i,j,k} = find(F1==F1_lvls(i) & F2==F2_lvls(j) & S==Subjs(k));
            CELLS{i,j,k} = Y(INDS{i,j,k});
            MEANS(i,j,k) = mean(CELLS{i,j,k});
        end
    end
end

% make tables (see table 18.1, p. 402)
AB = reshape(sum(MEANS,3),a,b); % across subjects
AS = reshape(sum(MEANS,2),a,n); % across factor 2
BS = reshape(sum(MEANS,1),b,n); % across factor 1

A = sum(AB,2); % sum across columns, so result is ax1 column vector
B = sum(AB,1); % sum across rows, so result is 1xb row vector
S = sum(AS,1); % sum across columns, so result is 1xs row vector
T = sum(sum(A)); % could sum either A or B or S, choice is arbitrary

% degrees of freedom
dfA = a-1;
dfB = b-1;
dfAB = (a-1)*(b-1);
dfS = n-1;
dfAS = (a-1)*(n-1);
dfBS = (b-1)*(n-1);
dfABS = (a-1)*(b-1)*(n-1);

% bracket terms (expected value)
expA = sum(A.^2)./(b*n);
expB = sum(B.^2)./(a*n);
expAB = sum(sum(AB.^2))./n;
expS = sum(S.^2)./(a*b);
expAS = sum(sum(AS.^2))./b;
expBS = sum(sum(BS.^2))./a;
expY = sum(sum(sum(MEANS.^2))); %sum(Y.^2);
expT = T^2 / (a*b*n);

% sums of squares
ssA = expA - expT;
ssB = expB - expT;
ssAB = expAB - expA - expB + expT;
ssS = expS - expT;
ssAS = expAS - expA - expS + expT;
ssBS = expBS - expB - expS + expT;
ssABS = expY - expAB - expAS - expBS + expA + expB + expS - expT;
ssTot = expY - expT;

% mean squares
msA = ssA / dfA;
msB = ssB / dfB;
msAB = ssAB / dfAB;
msS = ssS / dfS;
msAS = ssAS / dfAS;
msBS = ssBS / dfBS;
msABS = ssABS / dfABS;

% f statistic
fA = msA / msAS;
fB = msB / msBS;
fAB = msAB / msABS;

% p values
pA = 1-fcdf(fA,dfA,dfAS);
pB = 1-fcdf(fB,dfB,dfBS);
pAB = 1-fcdf(fAB,dfAB,dfABS);

% return values
stats = {'Source','SS','df','MS','F','p';...
         FACTNAMES{1}, ssA, dfA, msA, fA, pA;...
         FACTNAMES{2}, ssB, dfB, msB, fB, pB;...
         [FACTNAMES{1} ' x ' FACTNAMES{2}], ssAB, dfAB, msAB, fAB, pAB;...
         [FACTNAMES{1} ' x Subj'], ssAS, dfAS, msAS, [], [];...
         [FACTNAMES{1} ' x Subj'], ssBS, dfBS, msBS, [], [];...
         [FACTNAMES{1} ' x ' FACTNAMES{2} ' x Subj'], ssABS, dfABS, msABS, [], []};
 
 return

