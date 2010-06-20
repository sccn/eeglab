Resampling statistical toolkit
------------------------------
This package contains a set of functions for inferiential statistics
using resampling methods. Data are organized into arrays so multiple
tests can be run at once. This package is part of the EEGLAB
public software for electro-encephalography analysis and is distributed
under the GNU GPL license.

Requirements
------------
Matlab is required. For parametric statistics, the statistics Matlab
toolbox is required, although it is not required for resampling methods.

Quickstart
----------
Compute a paired t-test between 10 subjects in 2 conditions
for a time (100 points) x electrode (64 electrodes) design.

cond1 = rand(100,64,10); % the last dimension contains 10 subjects
cond2 = rand(100,64,10)+0.2; 

[t df p] = statcond({ cond1 cond2 }); % parametric test (uses statistic toolbox)
[t df p] = statcond({ cond1 cond2 }, 'mode', 'bootstrap', 'naccu', 1000); % bootstrap
[t df p] = statcond({ cond1 cond2 }, 'mode', 'perm', 'naccu', 1000); % permunation

Bootstrap and permutation take about 5.2 second on a Macbook pro at 2Gz
(that's 6.4 millions t-test and there is not a single loop). You then may 
use the included fdr.m function to correct for multiple comparaisons

p = fdr(p);

For more information see the help message of each m function.

Content
-------
statcond      - main function to compute paired or unpaired t-test, 1-way anova and 2-way anova
fdr           - implement False Detection Rate method to correct for multiple comparisons
teststat      - function testing the statcond function output
ttest_cell    - computes paired t-test
ttest2_cell   - computes unpaired t-test
anova1_cell   - Anova 1-way unpaired (same as anova1 from the statistics toolbox)
anova1rm_cell - Anova repeated measures 1-way (paired)
anova2_cell   - Anova 2-way unpaired (same as anova2 from the statistics toolbox)
anova1rm_cell - Anova repeated measures 2-way (paired)
combinedata   - combine different conditions, support function for statcond

The latest revision of these function may be checked out using SVN at

http://sccn.ucsd.edu/repos/software/eeglab/functions/statistics/

Arnaud Delorme, PhD
June 19th 2010
University of San Diego California, USA
CNRS, University of Toulouse Paul Sabatier, France
