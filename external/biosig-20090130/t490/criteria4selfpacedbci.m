function [X,O] = criteria4selfpacedbci(D,hdr,trig,cl,onset,duration,ref)
%    
%	X = criteria4selfpacedbci(D,trig,cl,onset,duration,ref)
%


%    $Id: criteria4selfpacedbci.m,v 1.1 2009-01-30 06:04:50 arno Exp $
%    Copyright (C) 2005 by Alois Schloegl <a.schloegl@ieee.org>	
%    This is part of the BIOSIG-toolbox http://biosig.sf.net/



%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


M = max(cl)+1;
[n,d] = size(D);
d=d+1;

[ix,iy] = meshgrid(trig,ref(1):ref(2));
t0 = ix(:)+iy(:);
CC.MD = repmat(NaN,[M,d,d]); 
CC.NN = repmat(NaN,[M,d,d]); 

[CC.MD(1,:,:),CC.NN(1,:,:)] = covm(D(t0,:),'E');
CC.datatype = 'classifier:statistical:ld3';
CC.Labels = 0:M-1;

kap = repmat(NaN,length(onset),length(duration)); 
for k1=1:length(onset)
for k2=1:length(duration)
[k1,k2],
	s = zeros(n,1); 
	for k=1:length(cl)
		s(trig(k)+onset:trig(k)+onset(k1)+duration(k2)) = cl(k); 
	end;
	for k3 = 2:M;
		[CC.MD(k3,:,:),CC.NN(k3,:,:)]=covm(D(find(s==CC.Labels(k3)),:),'E');
	end; 	
	R = test_sc(CC,D,[]);
	%R = test_sc(CC,D,'LD3',s+1);
	
	%H   = sparse(max(R.output,[],2),s+1,1,M,M);
	ix = ~any(isnan([R.classlabel',s]),2);
	H   = full(sparse(R.classlabel(ix)'+1,s(ix)+1,1,M,M));

	N   = sum(ix);
	p0  = sum(diag(H))/N;  %accuracy of observed agreement, overall agreement 

	p_i = sum(H); %sum(H,1);
	pi_ = sum(H'); %sum(H,2)';

	SA  = 2*diag(H)'./(p_i+pi_); % specific agreement 
	
	pe  = (p_i*pi_')/(N*N);  % estimate of change agreement

	px  = sum(p_i.*pi_.*(p_i+pi_))/(N*N*N);

	kap(k1,k2) = (p0-pe)/(1-pe);

end;
end;

kap,
H,
[k1,k2]=find(kap==max(kap(:)))
k1=k1(1);
k2=k2(1);

s = zeros(size(D,1),1); 
for k=1:length(cl)
	s(trig(k)+onset:trig(k)+onset(k1)+duration(k2)) = cl(k); 
end;
	for k3 = 2:M;
		[CC.MD(k3,:,:),CC.NN(k3,:,:)]=covm(D(find(s==CC.Labels(k3)),:),'E');
	end; 	
	R = test_sc(CC,D,[]);
	%H   = sparse(max(R.output,[],2),s+1,1,M,M);
	ix = ~any(isnan([R.classlabel(:),s(:)]),2);
	H   = full(sparse(R.classlabel(ix)'+1,s(ix)'+1,1,M,M))

%X.H   = sparse(ix,s,1,M,M));
X = CC; 
X.output = R.output;
[X.kappa, X.kapSD] = kappa(H);
X.tref   = ref; 
X.Onset  = onset(k1);
X.Offset = onset(k1)+duration(k2);
save 

if nargout>1,
	O = bci4eval(R.output,trig,cl,min([ref(:);min(onset);min(onset)+min(duration)]),max([ref(:);max(onset);max(onset)+max(duration)]),hdr.SampleRate)
end; 