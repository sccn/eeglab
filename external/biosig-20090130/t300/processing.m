function [Y,Z] = processing(MODE,X,Z)
% Data processing in the BIOSIG toolbox, 
% This functions can be used as template as well as a wrapper 
% around different signal processing methods.  
%
% [Y,Zf] = processing(MODE,X,Zi)
%
%    X input signal 
%    Zi input status (optional), 
%	initial condition
%	if empty or not available, this initializes the method
%    Y  output signal	
%    Zf final condition 
%    MODE is a struct and determines the signal processing method
%
%    MODE =
%	{'ECG_envelope',Fs}       used for QRS-detection [1]
%	{'ECG_envelope',Fs,Threshold}  used for QRS-detection [1]
%    	'xyz'		-none-	
%
%
% Reference(s):
% [1] M.-E. Nygards, L. Sörnmo, Delineation of the QRS complex using the envelope of the e.c.g
%         Med. & Biol. Eng. & Comput., 1983, 21, 538-547.
%
%


%	$Revision: 1.1 $
%	$Id: processing.m,v 1.1 2009-01-30 06:04:44 arno Exp $
%	Copyright (C) 2000-2003 by Alois Schloegl <a.schloegl@ieee.org>	

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


if nargin<3, Z=[]; end;

if 0,~isstruct(MODE)
        D.method = MODE;
else
        D = MODE;
end;


if 0,
        
elseif strcmp(D{1},'AAR')
        if length(D)<2, D{2}=[2,3]; end;
        if length(D)<3, D{3}=[6]; end;
        if length(D)<4, D{4}=2^(-55/8); end;
        if length(D)<5, D{5}=eye(D{3}); end;
        if length(D)<6, D{6}=zeros(1,D{3}); end;
        
        [Y,e,REV] = aar(X, D{2}, D{3}, D{4}, D{5}, D{6});
        
elseif strcmp(D{1},'ECG_envelope')
        if isempty(Z),	% Initialization
                if length(D)<2,
                        fprintf(2,'ERROR: SampleRate not defined.\n');        
                        return;
                else 
                        Fs = D{2};   	%.SampleRate;
                        if length(D)>2,
                                Z.TH = D{3};	% Threshold
                                Z.datatype='QRS_events';    
                                Z.QRS_event = [];
                        else
                                Z.datatype='ECG_envelope';    
                                Z.ECGenvelope = [];
                        end;
                end;
                T = .005; % s
                K = floor(T*Fs+1);
                PAR.HI = zeros(1,4*(K+1)-1);
                PAR.HI(1:2:length(PAR.HI)) = 2/pi*1./[-2*K-1:2:2*K+1];
                PAR.HI(2*K+2) = -j;
                PAR.HI(:) = -j*PAR.HI(:);
                PAR.trgL = ceil(0.04*Fs);    % 40ms
                
                Z.B2 = [1:PAR.trgL+1 PAR.trgL:-1:1];
                Z.A2 = sum(abs(Z.B2));
                
                Z.B1 = conv(PAR.HI,[ones(1,ceil(T*Fs)),-ones(1,ceil(T*Fs))]*1000/Fs);       
                Z.A1 = 1;
                
                Z.tix = 0; 
                Z.RESULT =[];

                [tmp, Z.Z1] = filter(Z.B1,Z.A1,X);
                [Y,   Z.Z2] = filter(Z.B2,Z.A2,abs(tmp));
                
        else
                [tmp, Z.Z1] = filter(Z.B1,Z.A1,X,Z.Z1);
                [Y,   Z.Z2] = filter(Z.B2,Z.A2,abs(tmp),Z.Z2);
        end;
        if isfield(Z,'TH'),
                Y = gettrigger(Y,Z.TH);
                Z.QRS_event = [Z.QRS_event; (Y + Z.tix)/Fs];
                Z.tix    = Z.tix + size(Y,1);
        else
                Z.ECGenvelope = [Z.ECGenvelope; Y];
        end;
        
        
elseif strcmp(D{1},'CORRCOEF')
        Y.datatype='CORRCOEF';
        DIM  = size(X,2);
        % SUM, N and SSQ 
        %[Z.SUM,Z.N,Z.M2] = sumskipnan(X',2);
        [Y.SUM,Y.N] = sumskipnan(X',2);
        
        % Correlation, using all possible values 
        [Y.CC,Y.nN] = covm(X,'M');
        
        % Correlation of NaN's
        [Y.nanCC,z.nan.nN] = covm(X,isnan(X),'M');
        
        % all samples with any NaN rejected 
        tmp  = X(~any(isnan(X),2),:);
        N0 = size(tmp,1);
        Y.C0 = [ones(N0,1),tmp]' * tmp;
        Y.N0 = N0;
        
        if isempty(Z),	% Initialization
                Z = Y;         
        else
                % SUM, N and SSQ 
                Z.SUM = Z.SUM + Y.SUM;
                Z.N   = Z.N   + Y.N;
                %Z.M2  = Z.M2  + z.M2;
                
                % Correlation, using all possible values 
                Z.CC = Z.CC + Y.CC;
                Z.nN = Z.nN + Y.nN;
                
                % Correlation of NaN's
                Z.nanCC = Z.nanCC + Y.nanCC;
                %Z.nannN = Z.nannN + Y.nannN;
                
                % all samples with any NaN rejected 
                Z.C0 = Z.C0 + Y.C0;
                Z.N0 = Z.N0 + Y.N0;
        end;        
        
else 
        
end;
