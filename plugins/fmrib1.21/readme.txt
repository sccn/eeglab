README File for fmrib EEGLAB plugin

This plugin for EEGLAB adds a menu item under 'Tools' called
'FMRIB Tools' for removing artifacts form EEG data collected with FMRI.

Instructions:
	Place the folder fmribX.X (X.X depends on version) inside the 
'plugins' folder of EEGLAB.  When you run EEGLAB the plugin  will be 
detected and installed.  You should see the following message in Matlab 
'eeglab: adding plugin "fmrib1.0b" (see >> help eegplugin_fmrib)'
when you start EEGLAB.

Tools:
	1. FASTR (fmri artifact slice template removal):  This tool only 
requires that EEG data in EEGLAB to have timing events for each FMRI slice 
acquisition.  It uses this information to robustly subtract the gradient 
artifact.  More information about the algorithm used is available in the 
plugin documentation.

       2. QRS detection.  This tool allows the detection of heart beat /QRS 
complexes (to be subsequently used in removing hear-related artifact) from 
a single channel of ECG data and stores the results as events in the EEGLAB 
event structure.  It has been working quite robustly even with quite bad ECGs.  Again, documentation and references is included in the plugin.

      3. Pulse artifact removal.   Uses events from (2) to remove pulse 
artifacts using different methods of constructing an artifact template. 




 Copyright (C) 2004 Rami K. Niazy, FMRIB Centre, University of Oxford
 rami@fmrib.ox.ac.uk

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA