function PLV = gaborPLV(x,y,fs,ntr,spt,f);
%
% GABOR_PLV   Computes the phase locking value (PLV) using the Gabor wavelet
%
% PLV = gabor_plv(x,y,fs,ntr,spt,f);
%
% Input parameters:
%   x, y ... Input signals
%   fs ... Sampling frequency
%   ntr ... Number of trials
%   spt ... Samples per trial
%   H ... Header containing information about x and y
%   f ... Frequency of interest
%
% Output parameters:
%   PLV ... Phase locking value between x and y

% Copyright by Clemens Brunner, original version by Joachim Reisinger
% $Revision: 1.1 $ $Date: 2009-01-30 06:04:44 $
% E-Mail: clemens.brunner@tugraz.at

% Create wavelet with frequency f
sig = 7/f;  % Standard deviation of Gaussian
T = 50/f;  % Length of wavelet
t1 = -T/2;
t2 =  T/2;
t = (t1:1/fs:t2-1/fs)';
if mod(length(t),2)
    t = (t1:1/fs:t2)';
end;
G = exp(- t.^2/(2*sig^2)) .* exp(j*2*pi*f*t);  % Gabor wavelet

SumTheta = 0;
Theta = [];

for tr = 1:ntr
    x_tr = x((tr-1)*spt+1:tr*spt);  % One trial
    y_tr = y((tr-1)*spt+1:tr*spt);
        
    Gx = conv(x_tr,G);  % Convolution, length(Gx)=length(x_tr) + length(G) - 1
    Gy = conv(y_tr,G);
    
    sGxy = size(Gx,1);  % Length of Gx or Gy
    % st = fix(size(t,1)/2);  % Länge der Wavelet Zeit/2 
    st = fix(size(t,1)/2);
    Gx = Gx(st:size(Gx,1)-st,:);   % Anfangs und Endbereichskorrektur
    Gy = Gy(st:size(Gy,1)-st,:);   % Anfangs und Endbereichskorrektur
      
    Phi_Gx = atan2(imag(Gx),real(Gx));       
    Phi_Gy = atan2(imag(Gy),real(Gy));
    Theta = Phi_Gx - Phi_Gy; % real
    Theta = unwrap(Theta);   % Phasensprungkorrektur
    SumTheta = SumTheta + exp(j*(Theta)); % komplex
end;

PLV = 1/ntr*(abs(SumTheta));