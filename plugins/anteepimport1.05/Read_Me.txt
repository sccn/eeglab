  ========================================
  September 2010
  MATLAB / EEGLAB import / plugin for ANT continuous data, trigger files and averages
  Version 1.07
  ========================================

  Update release adds support for 64-bit versions of Windows.
  For the latest version of our source-code, visit our project page here:
  http://sourceforge.net/projects/libeep/

  ========================================
  20 June 2008
  MATLAB / EEGLAB import / plugin for ANT continuous data, trigger files and averages
  Version 1.06
  ========================================

  Update release, fixed reading of triggers for epoch selection
  Deleted read_eep_trg.dll, restored default read_eep_trg.m
  (for general platform independency)

  ========================================
  15 January 2008
  MATLAB / EEGLAB import / plugin for ANT continuous data, trigger files and averages
  Version 1.05
  ========================================

  Update release, fixes a bug in read/write of the sampling rate
  (now always translates the "rate" variable to integer sampling frequency)
  Function write_eep_cnt.m still in beta version (not compatible with EEGLAB)
  Electrode positions file included ANT_WG_standard_346.ced, defines
  the electrode positions for the 10/5 percent system.

  ========================================
  26 September 2006
  MATLAB / EEGLAB import / plugin for ANT continuous data, trigger files and averages
  Version 1.04
  ========================================

  Update release, including Mac Mex files.
  Support for new format AVR (also still reads old-style)
  Support for Export of EEProbe CNT format (currently only in Windows)

  ========================================
  25 September 2006
  MATLAB / EEGLAB import / plugin for ANT continuous data, trigger files and averages
  Version 1.03
  ========================================

  Preliminary release, including Mac Mex files.
  Support for new format AVR (also still reads old-style)

  ========================================
  21 July 2005
  MATLAB / EEGLAB import / plugin for ANT continuous data, trigger files and averages
  Version 1.02
  ========================================

  Changed version number to the version of the EEGlab plugin for consistency.
  Update plugin menu for EEGlab -> File -> Import data -> From ANT EEProbe ...
  Fixed a bug in pop_loadeep.m, uilist variable, adding a space before the last {}

  ========================================
  1 February 2005
  MATLAB / EEGLAB import / plugin for ANT continuous data, trigger files and averages
  Version 4.4 beta
  ========================================

  This version fixes a bug that caused the Windows importer to crash.
  This error occurred when importing .cnt files where the <evt > chunk is missing.
  Further details can be obtained from ANT.

  ========================================
  17 October 2003
  IMPORTING ANT CNT DATA FILES INTO EEGLAB
  ========================================

  This compressed file contains Matlab scripts to import EEProbe CNT-RIFF raw data
  (compressed 32-bit file format) and EEProbe AVR average files (simple binary format)
  into EEGlab.

  Apart from these EEGlab import routines, additional Matlab import routines are provided to
  read the EEProbe TRG and REJ file formats (for triggers/events and rejection intervals).
  Also provided is a utility function (read_eep_trial.m) to read single-trial data.

  The importer tools use compiled Matlab MEX files. The enclosed version is compiled
  for both Microsoft Windows and Linux. When you are using a different platform,
  please request the corresponding MEX files from ANT at info@ant-software.nl

  For more information on using and writing plug-in see the EEGlab information.
  http://sccn.ucsd.edu/eeglab/

  EEGlab version 4.5 is necessary to run this import tool directly from EEGlab.
  The latest version of EEGlab can be downloaded as freeware from: http://sccn.ucsd.edu/eeglab/

  To use:

  1. Unpack this zip file into the plugin subdirectory in the main EEGlab directory on your computer.

  2. Start Matlab and type >> eeglab
     You will now see the following information in the command window:

     eeglab: adding plugin "eegplugin_eepimport"

     If this line does not show up, the package is not unpacked into the default EEGlab directory.
     The file "eegplugin_eepimport" has to be in the default EEG directory.

  4. The plug-in tool can be found in the "Import data" menu of EEGlab
     EEGlab>File>Import data>From ANT EEProbe .CNT file
     EEGlab>File>Import data>From ANT EEProbe .AVR file

  5. Please send feedback and suggestions to eeprobe@ant-neuro.com and/or info@ant-neuro.com

  Advanced Neuro Technology (ANT) BV, The Netherlands
  Internet: www.ant-neuro.com
  Email: info@ant-neuro.com

