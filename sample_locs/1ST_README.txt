This directory contains sample electrode location files for the standard 
10-20 and 10-10 Systems in EEGLAB format. Copy these files, 
then delete and/or add channels to match your electrode caps. You may
other standard electrode cap location files at the EEGLAB website.
".locs" files contain polar coordiantes. ".ced" files contain electrode 
positions in polar, 3-D cartesian and 3-D spherical coordinate frames.

These files can also be converted to other formats (e.g., BESA, etc.) 
using the EEGLAB channel editor: Under Matlab, move to this directory
(folder) and type 

        >> pop_chanedit(readlocs('Standard-10-20-Cap81.ced'));

Then, in the resulting pop-up window, press the button "Save others"
(Note: You may need to start 

        >> eeglab 

first to allow Matlab to find the channel editing function).
