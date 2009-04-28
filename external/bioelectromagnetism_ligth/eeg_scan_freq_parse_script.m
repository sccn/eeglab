% eeg_scan_freq_parse_script - generate Hz bands of Neuroscan AVG export file

clear; close all;

src = input('type working directory (eg, /ringo/eeg5/psmis/):\n','s');
file = input('\n\ntype fileprefix[.txt] to load (eg, text):\n','s');
ext = '.txt';

load (strcat(src, file, ext));
volt = eval(file); clear (file);


	% Generate fast fourier transform and power spectrum

ffx = fft(volt,256);         clear volt;
Pxx = ffx.* conj(ffx) / 256; clear ffx;


	% Plot power spectrum

%f = 4:4:128;
%plot(f,Pxx(1:32,3), '-+');


	% Mean power from 3-7 Hz, 8-13 Hz, 35-40 Hz, 40-45 Hz
	% across elec/trials

	% check that array coloumns correspond to correct f
	% especially given a lower resolution than 1 Hz.

theta  = Pxx( 1,:);
alpha  = Pxx( 2,:);
gamma1 = Pxx( 9,:);
gamma2 = Pxx(10,:);


	% create output data

th = fopen(strcat(src, file, '_th.dat'), 'wt');
al = fopen(strcat(src, file, '_al.dat'), 'wt');
g1 = fopen(strcat(src, file, '_g1.dat'), 'wt');
g2 = fopen(strcat(src, file, '_g2.dat'), 'wt');

fprintf(th, '%16.4f', theta);
fprintf(al, '%16.4f', alpha);
fprintf(g1, '%16.4f', gamma1);
fprintf(g2, '%16.4f', gamma2);

fclose('all');
