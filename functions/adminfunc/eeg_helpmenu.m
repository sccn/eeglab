% EEGLAB menus:
% File
%       Import data
%             Using EEGLAB functions and plugins
%                   From ASCII/float file or Matlab array - <a href="matlab:helpwin pop_importdata">pop_importdata</a>
%                   From Netstation binary simple file - <a href="matlab:helpwin pop_readegi">pop_readegi</a>
%                   From Netstation .MFF file - <a href="matlab:helpwin pop_readegimff">pop_readegimff</a>
%                   From Multiple seg. Netstation files - <a href="matlab:helpwin pop_readsegegi">pop_readsegegi</a>
%                   From Netstation Matlab files - <a href="matlab:helpwin pop_importegimat">pop_importegimat</a>
%                   From BCI2000 ASCII file - <a href="matlab:helpwin pop_loadbci">pop_loadbci</a>
%                   From Snapmaster .SMA file - <a href="matlab:helpwin pop_snapread">pop_snapread</a>
%                   From Neuroscan .CNT file - <a href="matlab:helpwin pop_loadcnt">pop_loadcnt</a>
%                   From Neuroscan .EEG file - <a href="matlab:helpwin pop_loadeeg">pop_loadeeg</a>
%                   From Biosemi BDF file (BIOSIG toolbox) - <a href="matlab:helpwin pop_biosig">pop_biosig</a>
%                   From EDF/EDF+/GDF files (BIOSIG toolbox) - <a href="matlab:helpwin pop_biosig">pop_biosig</a>
%                   From Biosemi BDF and EDF files (BDF plugin) - <a href="matlab:helpwin pop_readbdf">pop_readbdf</a>
%                   From BrainRT .SIG file - <a href="matlab:helpwin pop_readbrainrt">pop_readbrainrt</a>
%                   Import EKG file and adjust latencies - <a href="matlab:helpwin pop_importekgtxtfile">pop_importekgtxtfile</a>
%                   From ANT EEProbe .CNT file - <a href="matlab:helpwin pop_loadeep">pop_loadeep</a>
%                   From ANT EEProbe .AVR file - <a href="matlab:helpwin pop_loadeep_avg">pop_loadeep_avg</a>
%                   From BCI2000 .DAT file - <a href="matlab:helpwin pop_loadBCI2000">pop_loadBCI2000</a>
%                   From BIOPAC MATLAB files - <a href="matlab:helpwin pop_biopac">pop_biopac</a>
%                   From Brain Vis. Rec. .vhdr file - <a href="matlab:helpwin pop_loadbv">pop_loadbv</a>
%                   From Brain Vis. Anal. Matlab file - <a href="matlab:helpwin pop_loadbva">pop_loadbva</a>
%                   From CTF folder (MEG) - <a href="matlab:helpwin pop_ctf_read">pop_ctf_read</a>
%                   From ERPSS .RAW or .RDF file - <a href="matlab:helpwin pop_read_erpss">pop_read_erpss</a>
%                   From INStep .ASC file - <a href="matlab:helpwin pop_loadascinstep">pop_loadascinstep</a>
%                   From 4D .m4d pdf file - <a href="matlab:helpwin pop_read4d">pop_read4d</a>
%                   From nihonkohden data files - <a href="matlab:helpwin pop_nihonkohden">pop_nihonkohden</a>
%                   From Procom Infinity Text File - <a href="matlab:helpwin pop_importpi">pop_importpi</a>
%             Using the FILE-IO interface - <a href="matlab:helpwin pop_fileio">pop_fileio</a>
%             Using the BIOSIG interface - <a href="matlab:helpwin pop_biosig">pop_biosig</a>
%             Troubleshooting data formats...
%       Import epoch info
%             From Matlab array or ASCII file - <a href="matlab:helpwin pop_importepoch">pop_importepoch</a>
%             From Neuroscan .DAT file - <a href="matlab:helpwin pop_loaddat">pop_loaddat</a>
%       Import event info
%             From Matlab array or ASCII file - <a href="matlab:helpwin pop_importevent">pop_importevent</a>
%             From data channel - <a href="matlab:helpwin pop_chanevent">pop_chanevent</a>
%             From Presentation .LOG file - <a href="matlab:helpwin pop_importpres">pop_importpres</a>
%             From E-Prime ASCII (text) file - <a href="matlab:helpwin pop_importevent">pop_importevent</a>
%             From Neuroscan .ev2 file - <a href="matlab:helpwin pop_importev2">pop_importev2</a>
%       Export
%             Data and ICA activity to text file - <a href="matlab:helpwin pop_export">pop_export</a>
%             Weight matrix to text file - <a href="matlab:helpwin pop_expica">pop_expica</a>
%             Inverse weight matrix to text file - <a href="matlab:helpwin pop_expica">pop_expica</a>
%             Events to text file - <a href="matlab:helpwin pop_expevents">pop_expevents</a>
%             Data to EDF/BDF/GDF file - <a href="matlab:helpwin pop_writeeeg">pop_writeeeg</a>
%             Write Brain Vis. exchange format file - <a href="matlab:helpwin pop_writebva">pop_writebva</a>
%       Load existing dataset - <a href="matlab:helpwin pop_loadset">pop_loadset</a>
%       Save current dataset(s) - <a href="matlab:helpwin pop_saveset">pop_saveset</a>
%       Save current dataset as - <a href="matlab:helpwin pop_saveset">pop_saveset</a>
%       Clear dataset(s) - <a href="matlab:helpwin pop_delset">pop_delset</a>
%       Create study
%             Using all loaded datasets - <a href="matlab:helpwin pop_study">pop_study</a>
%             Browse for datasets - <a href="matlab:helpwin pop_study">pop_study</a>
%             Simple ERP STUDY - <a href="matlab:helpwin pop_studyerp">pop_studyerp</a>
%       Load existing study - <a href="matlab:helpwin pop_loadstudy">pop_loadstudy</a>
%       Save current study - <a href="matlab:helpwin pop_savestudy">pop_savestudy</a>
%       Save current study as - <a href="matlab:helpwin pop_savestudy">pop_savestudy</a>
%       Clear study
%       Memory and other options - <a href="matlab:helpwin pop_editoptions">pop_editoptions</a>
%       History scripts
%             Save dataset history script - <a href="matlab:helpwin pop_saveh">pop_saveh</a>
%             Save session history script - <a href="matlab:helpwin pop_saveh">pop_saveh</a>
%             Run script - <a href="matlab:helpwin pop_runscript">pop_runscript</a>
%       Quit - <a href="matlab:helpwin pop_saveh">pop_saveh</a>
% Edit
%       Dataset info - <a href="matlab:helpwin pop_editset">pop_editset</a>
%       Event fields - <a href="matlab:helpwin pop_editeventfield">pop_editeventfield</a>
%       Event values - <a href="matlab:helpwin pop_editeventvals">pop_editeventvals</a>
%       About this dataset - <a href="matlab:helpwin pop_comments">pop_comments</a>
%       Channel locations - <a href="matlab:helpwin pop_chanedit">pop_chanedit</a>
%       Select data - <a href="matlab:helpwin pop_select">pop_select</a>
%       Select data using events - <a href="matlab:helpwin pop_rmdat">pop_rmdat</a>
%       Select epochs or events - <a href="matlab:helpwin pop_selectevent">pop_selectevent</a>
%       Copy current dataset - <a href="matlab:helpwin pop_copyset">pop_copyset</a>
%       Append datasets - <a href="matlab:helpwin pop_mergeset">pop_mergeset</a>
%       Delete dataset(s) from memory - <a href="matlab:helpwin pop_delset">pop_delset</a>
%       Edit events & mark bad channels - <a href="matlab:helpwin pop_VisEd">pop_VisEd</a>
% Tools
%       Change sampling rate - <a href="matlab:helpwin pop_resample">pop_resample</a>
%       Filter the data
%             Basic FIR filter (new) - <a href="matlab:helpwin pop_eegfiltnew">pop_eegfiltnew</a>
%             Windowed sinc FIR filter - <a href="matlab:helpwin pop_firws">pop_firws</a>
%             Parks-McClellan (equiripple) FIR filter - <a href="matlab:helpwin pop_firpm">pop_firpm</a>
%             Moving average FIR filter - <a href="matlab:helpwin pop_firma">pop_firma</a>
%             Basic FIR filter (legacy) - <a href="matlab:helpwin pop_eegfilt">pop_eegfilt</a>
%             Basic FIR filter (new) - <a href="matlab:helpwin pop_eegfiltnew">pop_eegfiltnew</a>
%             Windowed sinc FIR filter - <a href="matlab:helpwin pop_firws">pop_firws</a>
%             Parks-McClellan (equiripple) FIR filter - <a href="matlab:helpwin pop_firpm">pop_firpm</a>
%             Moving average FIR filter - <a href="matlab:helpwin pop_firma">pop_firma</a>
%             Short non-linear IIR filter - <a href="matlab:helpwin pop_iirfilt">pop_iirfilt</a>
%       Re-reference - <a href="matlab:helpwin pop_reref">pop_reref</a>
%       Interpolate electrodes - <a href="matlab:helpwin pop_interp">pop_interp</a>
%       Reject continuous data by eye - <a href="matlab:helpwin pop_eegplot">pop_eegplot</a>
%       Extract epochs - <a href="matlab:helpwin pop_epoch">pop_epoch</a>
%       Remove baseline - <a href="matlab:helpwin pop_rmbase">pop_rmbase</a>
%       Run ICA - <a href="matlab:helpwin pop_runica">pop_runica</a>
%       Remove components - <a href="matlab:helpwin pop_subcomp">pop_subcomp</a>
%       Automatic channel rejection - <a href="matlab:helpwin pop_rejchan">pop_rejchan</a>
%       Automatic continuous rejection - <a href="matlab:helpwin pop_rejcont">pop_rejcont</a>
%       Automatic epoch rejection - <a href="matlab:helpwin pop_autorej">pop_autorej</a>
%       Reject data epochs
%             Reject data (all methods) - <a href="matlab:helpwin pop_rejmenu">pop_rejmenu</a>
%             Reject by inspection - <a href="matlab:helpwin pop_eegplot">pop_eegplot</a>
%             Reject extreme values - <a href="matlab:helpwin pop_eegthresh">pop_eegthresh</a>
%             Reject by linear trend/variance - <a href="matlab:helpwin pop_rejtrend">pop_rejtrend</a>
%             Reject by probability - <a href="matlab:helpwin pop_jointprob">pop_jointprob</a>
%             Reject by kurtosis - <a href="matlab:helpwin pop_rejkurt">pop_rejkurt</a>
%             Reject by spectra - <a href="matlab:helpwin pop_rejspec">pop_rejspec</a>
%             Export marks to ICA reject - <a href="matlab:helpwin eeg_checkset">eeg_checkset</a>
%             Reject marked epochs - <a href="matlab:helpwin pop_rejepoch">pop_rejepoch</a>
%       Reject data using ICA
%             Reject components by map - <a href="matlab:helpwin pop_selectcomps">pop_selectcomps</a>
%             Reject data (all methods) - <a href="matlab:helpwin pop_rejmenu">pop_rejmenu</a>
%             Reject by inspection - <a href="matlab:helpwin pop_eegplot">pop_eegplot</a>
%             Reject extreme values - <a href="matlab:helpwin pop_eegthresh">pop_eegthresh</a>
%             Reject by linear trend/variance - <a href="matlab:helpwin pop_rejtrend">pop_rejtrend</a>
%             Reject by probability - <a href="matlab:helpwin pop_jointprob">pop_jointprob</a>
%             Reject by kurtosis - <a href="matlab:helpwin pop_rejkurt">pop_rejkurt</a>
%             Reject by spectra - <a href="matlab:helpwin pop_rejspec">pop_rejspec</a>
%             Export marks to data reject - <a href="matlab:helpwin eeg_checkset">eeg_checkset</a>
%             Reject marked epochs - <a href="matlab:helpwin pop_rejepoch">pop_rejepoch</a>
%       Locate dipoles using DIPFIT 2.x
%             Head model and settings - <a href="matlab:helpwin pop_dipfit_settings">pop_dipfit_settings</a>
%             Coarse fit (grid scan) - <a href="matlab:helpwin pop_dipfit_gridsearch">pop_dipfit_gridsearch</a>
%             Fine fit (iterative) - <a href="matlab:helpwin pop_dipfit_nonlinear">pop_dipfit_nonlinear</a>
%             Autofit (coarse fit, fine fit & plot) - <a href="matlab:helpwin pop_multifit">pop_multifit</a>
%             Plot component dipoles - <a href="matlab:helpwin pop_dipplot">pop_dipplot</a>
%       Peak detection using EEG toolbox
%       FMRIB Tools
%             FASTR: Remove FMRI gradient artifacts - <a href="matlab:helpwin pop_fmrib_fastr">pop_fmrib_fastr</a>
%             Detect QRS events - <a href="matlab:helpwin pop_fmrib_qrsdetect">pop_fmrib_qrsdetect</a>
%             Remove pulse artifacts - <a href="matlab:helpwin pop_fmrib_pas">pop_fmrib_pas</a>
%       Locate dipoles using LORETA
%             Export components to LORETA - <a href="matlab:helpwin pop_eeglab2loreta">pop_eeglab2loreta</a>
% Plot
%       Channel locations
%             By name - <a href="matlab:helpwin topoplot">topoplot</a>
%             By number - <a href="matlab:helpwin topoplot">topoplot</a>
%       Channel data (scroll) - <a href="matlab:helpwin pop_eegplot">pop_eegplot</a>
%       Channel spectra and maps - <a href="matlab:helpwin pop_spectopo">pop_spectopo</a>
%       Channel properties - <a href="matlab:helpwin pop_prop">pop_prop</a>
%       Channel ERP image - <a href="matlab:helpwin pop_erpimage">pop_erpimage</a>
%       Channel ERPs
%             With scalp maps - <a href="matlab:helpwin pop_timtopo">pop_timtopo</a>
%             In scalp/rect. array - <a href="matlab:helpwin pop_plottopo">pop_plottopo</a>
%       ERP map series
%             In 2-D - <a href="matlab:helpwin pop_topoplot">pop_topoplot</a>
%             In 3-D - <a href="matlab:helpwin pop_headplot">pop_headplot</a>
%       Sum/Compare ERPs - <a href="matlab:helpwin pop_comperp">pop_comperp</a>
%       Component activations (scroll) - <a href="matlab:helpwin pop_eegplot">pop_eegplot</a>
%       Component spectra and maps - <a href="matlab:helpwin pop_spectopo">pop_spectopo</a>
%       Component maps
%             In 2-D - <a href="matlab:helpwin pop_topoplot">pop_topoplot</a>
%             In 3-D - <a href="matlab:helpwin pop_headplot">pop_headplot</a>
%       Component properties - <a href="matlab:helpwin pop_prop">pop_prop</a>
%       Component ERP image - <a href="matlab:helpwin pop_erpimage">pop_erpimage</a>
%       Component ERPs
%             With component maps - <a href="matlab:helpwin pop_envtopo">pop_envtopo</a>
%             With comp. maps (compare) - <a href="matlab:helpwin pop_envtopo">pop_envtopo</a>
%             In rectangular array - <a href="matlab:helpwin pop_plotdata">pop_plotdata</a>
%       Sum/Compare comp. ERPs - <a href="matlab:helpwin pop_comperp">pop_comperp</a>
%       Data statistics
%             Channel statistics - <a href="matlab:helpwin pop_signalstat">pop_signalstat</a>
%             Component statistics - <a href="matlab:helpwin pop_signalstat">pop_signalstat</a>
%             Event statistics - <a href="matlab:helpwin pop_eventstat">pop_eventstat</a>
%       Time-frequency transforms
%             Channel time-frequency - <a href="matlab:helpwin pop_newtimef">pop_newtimef</a>
%             Channel cross-coherence - <a href="matlab:helpwin pop_newcrossf">pop_newcrossf</a>
%             Component time-frequency - <a href="matlab:helpwin pop_newtimef">pop_newtimef</a>
%             Component cross-coherence - <a href="matlab:helpwin pop_newcrossf">pop_newcrossf</a>
%       Cluster dataset ICs - <a href="matlab:helpwin pop_miclust">pop_miclust</a>
% Study
%       Edit study info - <a href="matlab:helpwin pop_study">pop_study</a>
%       Select/Edit study design(s) - <a href="matlab:helpwin pop_studydesign">pop_studydesign</a>
%       Precompute channel measures - <a href="matlab:helpwin pop_precomp">pop_precomp</a>
%       Plot channel measures - <a href="matlab:helpwin pop_chanplot">pop_chanplot</a>
%       Precompute component measures - <a href="matlab:helpwin pop_precomp">pop_precomp</a>
%       PCA clustering (original)
%             Build preclustering array - <a href="matlab:helpwin pop_preclust">pop_preclust</a>
%             Cluster components - <a href="matlab:helpwin pop_clust">pop_clust</a>
%       Edit/plot clusters - <a href="matlab:helpwin pop_clustedit">pop_clustedit</a>
%       Cluster components by correlation (CORRMAP) - <a href="matlab:helpwin pop_corrmap">pop_corrmap</a>
% Datasets
%       Select multiple datasets - <a href="matlab:helpwin pop_chansel">pop_chansel</a>
% Help
%       Upgrade to the Latest Version
%       About EEGLAB
%       About EEGLAB help - <a href="matlab:helpwin eeg_helphelp">eeg_helphelp</a>
%       EEGLAB menus - <a href="matlab:helpwin eeg_helpmenu">eeg_helpmenu</a>
%       EEGLAB functions
%             Admin functions - <a href="matlab:helpwin eeg_helpadmin">eeg_helpadmin</a>
%             Interactive pop_ functions - <a href="matlab:helpwin eeg_helppop">eeg_helppop</a>
%             Signal processing functions - <a href="matlab:helpwin eeg_helpsigproc">eeg_helpsigproc</a>
%             Group processing (STUDY) functions - <a href="matlab:helpwin eeg_helpstudy">eeg_helpstudy</a>
%             Time-frequency functions - <a href="matlab:helpwin eeg_helptimefreq">eeg_helptimefreq</a>
%             Statistics functions - <a href="matlab:helpwin eeg_helpstatistics">eeg_helpstatistics</a>
%             Graphic interface builder functions - <a href="matlab:helpwin eeg_helpgui">eeg_helpgui</a>
%             Misceleanous functions (command line only) - <a href="matlab:helpwin eeg_helpmisc">eeg_helpmisc</a>
%       EEGLAB license
%       Web tutorial
%       Email the EEGLAB team
