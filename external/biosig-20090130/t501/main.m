
%General:
% To visualize e.g. a cross-spectrum or coherence in a 'head-in-head', 
% plot you  need some functional data (e.g. a cross-spectral matrx cs)
% and the locations of the channels. A location-variable for N channels 
% is either an Nx3 matrix, if the locations are given in 3D, or a 
% Nx2 matrix if they are given in 2D. 

%  $Id: main.m,v 1.1 2009-01-30 06:04:51 arno Exp $ 
%  Copyright (C) 2003,2004 Guido Nolte
%  Adapted by Alois Schloegl, 2007 
%  This function is part of the BioSig project
%  http://biosig.sf.net/	

% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the  License, or (at your option) any later version.
% 


%load eeg_sens; 
%load some eeg 3D-channels locations; 
testfile = '/home/neuro/data/BCI/bbciRaw/Friederike_06_06_06/arteFriederike.vhdr';  
HDR=sopen(testfile);HDR=sclose(HDR);

ix = all(~isnan(HDR.ELEC.XYZ),2);  % use only electrodes with known position
%HDR.Label(~ix);

clear para;para.rot=180; %setting this rotates locations by 180 degrees
locs_2D=mk_sensors_plane(HDR.ELEC.XYZ(ix,:),para); %make a variable locs_2D 
                                         % which contains informations 
                                         % how to plot coherence or covariance matrices,
                                         % the second argument para is optional 
                                         % If set as above, channel locations are rotated by 180 degrees
                                         % So, the channels close to the front are  up. 
% mk_sensors_plane does the following: 
% 1. it maps 3D locations to 2D locations; if the first argument (here eeg_sens) is an 
%     an Nx2 matrix (instead of Nx3 matrix) this step is omitted. 
% 2. it slightly shifts electrode locations to avoid overlapping spheres in these coherence plots
%     shown below.  This change is only effective for the locations of the circles and not for 
%     the locations of the electrodes within each circle. 
% mk_sensors_plane shows two figures, one with the 'original' 
% channel locations in 2D and one with the shifted locations where 
% the circles will be put.  
%
% The output locs_2D is an Nx5 matrix, the first column shows the indices of displayed 
% electrodes. If positive, it is used both as location for a circle and as electrode within 
% a circle. If negative, it is only used within a circle. (You need this when you have 
% very many channels). The second and third column indicate 'true' sensor locations 
% in 2 dimensions. These are the locations where (e.g.) coherence is plotted within 
% a small circle. The fourth and fifth column are the locations of the small circles. 
% They are similar to the 'true' locations but slightly shifted to avoid
% overlapping spheres. 


%load cs_eeg %load an example for a cross-spectrum (this is real data at 10Hz) 
cs = randn(size(locs_2D,1)); % use random data for demo 
figure;plot_coupling(imag(cs),locs_2D);  %show the imaginary part of the cross-spectrum

%If you want to show a potential or field you can use 
% the program showfield_general.
% E.g., to show the eigenvector of largest eigenvalue of the real 
% part of CS, calculated by
[u,s,v]=svd(real(cs));umax=u(:,1);
% then you can display it as:
figure;showfield_general(umax,locs_2D);
% 


% Second example
load locs_247 %load 3-D coordinates for 247 MEG channels.

clear para;
para.nin=50; % 247 channels are too many to display.
             % Setting pars.nin=50 selects a subset of 50 
             %  channels for small circles. The selection is 
             % done such that the minimum distance (in 3D) between selected channels 
             % is as large as possible. 
             % If you want 2D distances, first run this program and rerun it with 
             % locs_2D(2:3,:) as input. 

para.rot=90; % make a 90 degrees rotation 
locs_2D=mk_sensors_plane(locs_247,para); %make location file 
             
%load cs; % load example covariance matrix (This is a made up one)) 
cs = randn(size(locs_2D,1)); % use random data for demo 
figure;plot_coupling(cs,locs_2D,para); %plots 50 channels as 'head-in-head' with the above changes 


%to change parameters of the plotting you can write a structure, here named  'para',
% and include as a third argument. As examples, here's a list of things you can. 
% You can also set only a subset of parameters
para.resolution=10; %decreases the resolution; default is 25 
para.global_size=1.1; % increases the size of the whole set of circles by factor 1.1
para.relative_size=1.1;
para.head_right=-.02; %shifts the drawn head slightly to the right
para.head_up=-.03;  %shifts the drawn head slightly down 
para.scales=[-2,2];  % changes the colormap to range from -2 to 2
figure;plot_coupling(cs,locs_2D,para); %plots 30 channels as 'head-in-head' with the above changes 

% Finally, if you don't won't to shift circles you can use: 
clear para;para.rot=90;
para.circle_shift=0; %omits
locs_2D=mk_sensors_plane(locs_247,para); %make location file 

