function plotHRSpectrum(r)
% Plots the spectrum of the heart rate.
%
% Plots the heart rate spectrum as calculated by calcHRSpectrum.
%
% Usage:
%   plotHRSpectrum(r);
%
% Input parameters:
%   r ... Input structure as obtained by calcHRSpectrum.m

% Copyright by Robert Leeb, Clemens Brunner
% $Revision: 1.1 $ $Date: 2009-01-30 06:04:51 $
% E-Mail: robert.leeb@tugraz.at

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

p = r{1}.p;
f = r{1}.f;

dispFreq = 0.5;  % Plot up to this frequency

idx = find(f<=dispFreq);
semilogy(f(idx), p(idx), 'k')
title('Heart Rate Spectrum', 'FontSize', 16, 'Interpreter', 'none');
xlabel('Frequency (Hz)');
ylabel('PSD (bpm^2/Hz)');
set(gca, 'FontSize', 16);
set(gca, 'XGrid', 'on');
v=axis;
axis([v(1) v(2) 0.001 10]);

set(gcf, 'PaperOrientation', 'portrait');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperType', 'A4');
set(gcf, 'PaperPosition', [1, 1, 19, 27.7]);