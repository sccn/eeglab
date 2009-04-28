% script to plot an ERP and topographic map
clear
p = eeg_toolbox_defaults
p = eeg_open(p)
p = elec_open(p)
p = mesh_open(p)
p = gui_erp_plot(p)
p = eeg_contours_engine(p)
