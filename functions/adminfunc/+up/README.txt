Information for when releasing a new EEGLAB version.
This message is for EEGLAB developpers only.

This plugins allows to perform automatic updates.
To determine the current version, type in

[dummy eeglabVersionNumber currentReleaseDateString] = eeg_getversion;
eeglabUpdater = up.updater(eeglabVersionNumber, 'http://sccn.ucsd.edu/eeglab/updater/latest_version.php', 'EEGLAB', currentReleaseDateString);
eeglabUpdater.checkForNewVersion({'eeglab_event' 'setup'});
eeglabUpdater

Copy the information about the new version to

/home/www/eeglab/updater/latestRelease.xml

Arnaud Delorme - April 23rd, 2012
