% Benchmark for the BioSig toolbox 
%  The benchmark test performs a typical data analysis task.  
%  Except for data loading, it tests the performance of the computational speed 
%  and can be used to compare the performance of different platforms 
%  
%  Requirements: 
%  Octave 2.1 or higher, or Matlab 
%  BioSig4OctMat from http:/biosig.sf.net/
%
%  run bench_biosig and compare the results at other platforms  
%   http://bci.tugraz.at/~schloegl/biosig/bench/
% 
%  Send your benchmark result to <a.schloegl@ieee.org>
%


%	$Id: bench_biosig.m,v 1.1 2009-01-30 06:04:39 arno Exp $
%	Copyright (C) 2005,2006 by Alois Schloegl <a.schloegl@ieee.org>	
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



if ~exist('l1.gdf','file') % download test file if not available 
        system('wget http://hci.tugraz.at/schloegl/bci/bci7/l1.gdf');
end; 

tic; K=0; 

%[signal,HDR]=sload({[p,'x_train'],[p,'x_test']});
[signal,HDR]=sload('l1.gdf');

K=K+1;jo{K}='sload l1.gdf'; t(K)=toc 
z = zscore(signal);

try
bp = bandpower(z,HDR.SampleRate);
bp(bp<-10)=NaN;
K=K+1;jo{K}='bandpower'; t(K)=toc 
catch, end; 

[SIGMA1,PHI1,OMEGA1,m01,m11,m21] = wackermann(z(:,1:2),HDR.SampleRate);
[SIGMA2,PHI2,OMEGA2,m02,m12,m22] = wackermann(z(:,2:3),HDR.SampleRate);
K=K+1;jo{K}='wackermann'; t(K)=toc 

[a,f,s] = barlow(z,HDR.SampleRate);
BARLOW = [a,f,s]; 
K=K+1;jo{K}='barlow'; t(K)=toc 

[A,M,C] = hjorth(z,HDR.SampleRate);
HJORTH = [A,M,C];
K=K+1;jo{K}='hjorth'; t(K)=toc 

uc  = 30:5:80;
        a = [];
        for ch = 1:size(signal,2),
                fprintf(1,'.%i',ch);
                INI.MOP = [0,3,0];
                INI.UC  = 2^(-uc(7)/8);
                
                X = tvaar(signal(:,ch),INI);
                X = tvaar(signal(:,ch),X);
                
                a = [a,X.AAR];
                K = K+1; jo{K} = ['aar #',int2str(ch)]; t(K)=toc
        end


try,
CC1 = findclassifier(bp,HDR.EVENT.POS(HDR.EVENT.TYP==hex2dec('300'))-1,HDR.Classlabel,reshape(1:1152,16,72)',[1:72]'>24,'LD3');
CC1.TSD.T = CC1.TSD.T/HDR.SampleRate;
K=K+1; jo{K}='findclassifier bp'; t(K)=toc 
catch,end;


CC2 = findclassifier(BARLOW,HDR.EVENT.POS(HDR.EVENT.TYP==hex2dec('300'))-1,HDR.Classlabel,reshape(1:1152,16,72)',[1:72]'>24,'LD3');

CC2.TSD.T = CC2.TSD.T/HDR.SampleRate;
K=K+1; jo{K}='findclassifier barlow'; t(K)=toc 
CC3 = findclassifier(HJORTH,HDR.EVENT.POS(HDR.EVENT.TYP==hex2dec('300'))-1,HDR.Classlabel,reshape(1:1152,16,72)',[1:72]'>24,'LD3');
CC3.TSD.T = CC3.TSD.T/HDR.SampleRate;
K=K+1; jo{K}='findclassifier hjorth'; t(K)=toc 
CC4 = findclassifier(a,HDR.EVENT.POS(HDR.EVENT.TYP==hex2dec('300'))-1,HDR.Classlabel,reshape(1:1152,16,72)',[1:72]'>24,'LD3');
CC4.TSD.T = CC4.TSD.T/HDR.SampleRate;
K=K+1; jo{K}='findclassifier aar'; t(K)=toc 
CC5 = findclassifier([SIGMA1,PHI1,OMEGA1,SIGMA2,PHI2,OMEGA2],HDR.EVENT.POS(HDR.EVENT.TYP==hex2dec('300'))-1,HDR.Classlabel,reshape(1:1152,16,72)',[1:72]'>24,'LD3');
CC5.TSD.T = CC5.TSD.T/HDR.SampleRate;
K=K+1; jo{K}='findclassifier Wackermann'; t(K)=toc 

if exist('OCTAVE_VERSION','builtin');
	om = 'Octave'; 
elseif 1, 
	om = 'Matlab';
else
	om = 'FreeMat';
end;	
outfile = sprintf('bench_biosig1.75+_%s_%s_%s.mrk',computer,om,version); 

try
unix(['cat /proc/cpuinfo >"',outfile,'"'])
catch
end;
fid = fopen(outfile,'a'); 
fprintf(fid,'\n\nDate:\t%s\n',date);
fprintf(fid,'Revision:\t$Id: bench_biosig.m,v 1.1 2009-01-30 06:04:39 arno Exp $\n');
fprintf(fid,'Computer:\t%s\nSoftware:\t%s\nVersion:\t%s\n',computer,om,version);

tmp = [diff([0,t(:)']);t(:)']'; 
for k=1:K, 
	fprintf(fid,'%25s:\t%6.1f\t%6.1f\n',jo{k},tmp(k,:));
end;
fclose(fid); 

