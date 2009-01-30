To make these function work under the EEGLAB interface
you must uncompress the ctf folder in the plugin directory
of EEGLAB. Type

   cd eeglab4.311/plugins
   cvs -z3 -d:pserver:anonymous@cvs.sourceforge.net:/cvsroot/eeg checkout ctf

   [use "move" instead of "mv" under windows]

What these commands do is that they download all the files into 
a folder named ctf under the plugin sub-directory of your version of
EEGLAB. 

Arnaud Delorme, July 9, 2004
