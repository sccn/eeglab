

# What is EEGLAB?
EEGLAB is an open source signal processing environment for electrophysiological signals running on Matlab and Octave (command line only for Octave). This folder contains original Matlab functions from the EEGLAB (formerly ICA/EEG) Matlab toolbox, all released under the Gnu public license (see eeglablicence.txt). See the EEGLAB tutorial and reference paper (URLs given below) for more information.

# Installing/cloning
**Recommended:** Download the official EEGLAB release from https://sccn.ucsd.edu/eeglab/download.php

**Do not download a ZIP file directly from GIT as it will not contain EEGLAB submodules**. Instead clone the reposity while pulling EEGLAB sub-modules.

```
git clone --recurse-submodules https://github.com/sccn/eeglab.git
```

If you forgot to clone the submodule, go to the eeglab folder and type

```
git submodule update --init --recursive --remote
git pull --recurse-submodules
```

# Sub-directories:

 - /functions - All distributed EEGLAB functions (admin, sigproc, pop, misc)
 - /plugins   - Directory to place all downloaded EEGLAB plug-ins. dipfit (1.0) is present by default
 - /sample_data -  Miscellaneous EEGLAB data using in tutorials and references
 - /sample_locs -  Miscellaneous standard channel location files (10-10, 10-20). See the EEGLAB web site http://sccn.ucsd.edu/eeglab/ for more.

# To use EEGLAB: 

1. Start Matlab

2. Use Matlab to navigate to the folder containing EEGLAB

3. Type "eeglab" at the Matlab command prompt ">>" and press enter

3. Open the main EEGLAB tutorial page (http://sccn.ucsd.edu/wiki/EEGLAB_Wiki)

4. Please send feedback and suggestions to: eeglab@sccn.ucsd.edu

# In publications, please reference:

Delorme, A., & Makeig, S. (2004). EEGLAB: an open source toolbox for analysis of single-trial EEG dynamics including independent component analysis. Journal of neuroscience methods, 134(1), 9-21. (See article [here](http://sccn.ucsd.edu/eeglab/download/eeglab_jnm03.pdf))
 
# Documentation:

EEGLAB documentation is available on the EEGLAB wiki (see http://sccn.ucsd.edu/wiki/EEGLAB_Wiki for more details).
