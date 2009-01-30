function [B,A]=paynter(RC,fs,mode)
% PAYNTER returns the filter coefficients for a Paynter Filter
%   Usually, the filter is applied to the rectified electromyogram (EMG).
%   Then, the output is amplitude-demodulated EMG 
%
%  The amplitude demodulated EMG can be obtained through
%      y = filter(B,A,abs(x));
%   with 
%       [B,A]=paynter(tau,fs)
%       [B,A]=paynter(tau,fs,'modified')
%       [B,A]=paynter(tau,fs,'bessel-modified')
%
%       tau	time constant (in [s])
%       fs	sampling rate (in [Hz])
%
%
% REFERENCE(S):
% [1] Platt, Ronald S., Eric A. Hajduk, Manuel Hulliger, and Paul A. Easton. 
%       A modified Bessel filter for amplitude demodulation of respiratory electromyograms. 
%       J. Appl. Physiol. 84(1): 378-388, 1998.
%       available online:  http://jap.physiology.org/cgi/content/full/84/1/378
% [2] Gottlieb, G.L. and Agarwal. 
%       Filtering of Electromyographic Signals. Am.J.Physical Medicine 49(3):142-146, 1970.
% [3] Hopp, F.A., J.L. Seagard, and J.P. Kampine. 
%     Comparison of four methods of averaging nerve activity. Am. J. Physiology 251:R700-R711, 1986.
% [4] Bruce, E. N., M. D. Goldman, and J. Mead. 
%       A digital computer technique for analyzing respiratory muscle EMGs. 
%       J. Appl. Physiol. 43: 551-556, 1977 
%       available online:  http://jap.physiology.org/cgi/content/abstract/43/3/551

%	$Id: paynter.m,v 1.1 2009-01-30 06:04:44 arno Exp $
%	Copyright (C) 2000, 2001, 2004, 2006 by Alois Schloegl <a.schloegl@ieee.org>	
%    	This is part of the BIOSIG-toolbox http://biosig.sf.net/

% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% Version 2 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 59 Temple Place - Suite 330,
% Boston, MA  02111-1307, USA.


if nargin<3,
        mode='pay';
end;

if length(mode)<3,
        fprintf(2,'Invalid mode in PAYNTER.M\n');
elseif mode(1:3)=='mod'
        mode='mod';
elseif mode(1:3)=='bes'
        mode='bes';
end;

if mode=='pay',
	%[B,A]=bilinear(1,rev([1, 3.2*RC, 4*RC*RC, 3.2*RC^3]),fs);
	[B,A]=bilinear(1,[1, 3.2*RC, 4*RC*RC, 3.2*RC^3],fs);
elseif mode=='mod'
	%[B,A]=bilinear(rev([1,0,RC*RC]),rev([1, 3.2*RC, 4*RC*RC, 3.2*RC^3]),fs);
	[B,A]=bilinear([1,0,RC*RC],[1, 3.2*RC, 4*RC*RC, 3.2*RC^3],fs);
elseif mode=='bes'
	B(1) = 1.6408979251842E-01; A(1) = 1.6408979251842E-01
	B(2) = 0;	A(2) = 1.1486285476290E+00
	B(3) = 1.3754616767430E-01;  A(3) = 3.7109537692628E+00
	B(4) = 0;	A(4) = 7.2157434402333E+00
	B(5) = 3.2736655946486E-02;  A(5) = 9.1836734693878E+00
	B(6) = 0;	A(6) = 7.7142857142857E+00
	B(7) = 1.9259308298786E-03;  A(7) = 4.0000000000000E+00
	B(8) = 0;	A(8) = 1.0000000000000E+00
end;

