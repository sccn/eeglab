

# What is EEGLAB?
EEGLAB is an open source signal processing environment for electrophysiological signals running on Matlab and Octave (command line only for Octave). This folder contains original Matlab functions from the EEGLAB (formerly ICA/EEG) Matlab toolbox, all released under the Gnu public license (see eeglablicence.txt). See the EEGLAB tutorial and reference paper (URLs given below) for more information.

# Instaling/cloning
Make sure you include sub-modules when you clone

git clone --recurse-submodules https://github.com/eeglabdevelopers/eeglab.git

If you forgot to clone the submodule, go to the eeglab folder and type

git submodule update --init --recursive

# Sub-directories:

 - /functions - All distributed EEGLAB functions (admin, sigproc, pop, misc)
 - /plugins   - Directory to place all downloaded EEGLAB plug-ins. dipfit (1.0) is present by default
 - /sample_data -  Miscellaneous EEGLAB data using in tutorials and references
 - /sample_locs -  Miscellaneous standard channel location files (10-10, 10-20). See the EEGLAB web site http://sccn.ucsd.edu/eeglab/ for more.

# To use EEGLAB: 

1. Place the Matlab functions in a directory ($DIR) and add $DIR/eeglabxx to (Unix) your matlabpath environment variable. 
   Else, within Matlab >> addpath('full_path_here')

2. Optional: Edit file "icadefs.m" which specifies various limits and constants used in EEGLAB functions.

3. Then start Matlab and type >> eeglab

4. Open the main EEGLAB tutorial page (your downloaded "eeglabtut.html",
   else browse http://sccn.ucsd.edu/wiki/EEGLAB_Wiki)

5. Please send feedback and suggestions to: eeglab@sccn.ucsd.edu

# In publications, please reference:

Delorme, A., Makeig, S. (in press) EEGLAB: An open source toolbox for analysis of single-trial EEG dynamics including independent component analysis. Journal of Neuroscience Methods. http://sccn.ucsd.edu/eeglab/download/eeglab_jnm03.pdf
 
Consider contributing to your functions and creativity to EEGLAB open source development (see http://sccn.ucsd.edu/wiki/EEGLAB_Wiki for more details).
