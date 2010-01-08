
This dir contains original Matlab functions from the EEGLAB (formerly ICA/EEG)
Matlab toolbox, all released under the Gnu public license (see eeglablicence.txt). 
See the EEGLAB tutorial and reference paper (URLs given below) for more information.

Sub-directories:

 /functions - All distributed EEGLAB functions (admin, sigproc, pop, misc)
 /plugins   - Directory to place all downloaded EEGLAB plug-ins. dipfit (1.0)
              is present by default
 /sample_data -  Miscellaneous EEGLAB data using in tutorials and references
 /sample_locs -  Miscellaneous standard channel location files (10-10, 10-20)
              See the EEGLAB web site http://sccn.ucsd.edu/eeglab/ for more.

To use EEGLAB: 

1. Place the Matlab functions in a directory ($DIR) and add 
   $DIR/eeglab4.x to (Unix) your matlabpath environment variable. 
   Else, within Matlab 
                      >> addpath('full_path_here')

2. Optional: Edit file "icadefs.m" under the function directory to specify the location of
   the faster binary "ica" function (equivalent to Matlab runica() and called from within
   Matlab by binica()). This requires a (recommended) separate download from 
                       http://sccn.ucsd.edu/eeglab/binica/
   Also add the path to the EEGLAB tutorial. This requires another (recommended) download. See
                       http://sccn.ucsd.edu/eeglab/
   File "icadefs.m" also specifies various limits and constants used in EEGLAB functions.

3. Then start Matlab and type >> eeglab

4. Open the main EEGLAB tutorial page (your downloaded "eeglabtut.html",
   else browse http://sccn.ucsd.edu/eeglab/eeglabtut.html)

5. Please send feedback and suggestions to: eeglab@sccn.ucsd.edu

6. In publications, please reference:

Delorme, A., Makeig, S. (in press) EEGLAB: An open source toolbox for analysis of single-trial 
EEG dynamics including independent component analysis. Journal of Neuroscience Methods.
               http://sccn.ucsd.edu/eeglab/download/eeglab_jnm03.pdf
 
Enjoy using EEGLAB to explore and analyze your data. Consider contributing to your functions
and creativity to EEGLAB open source development (see http://sccn.ucsd.edu/eeglab for more details).

Arno Delorme & Scott Makeig
Fri Sep  6 11:44:29 PDT 2002
