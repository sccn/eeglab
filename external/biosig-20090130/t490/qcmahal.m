function [Q,D,R,ECM]=qcmahal(ECM,ECM2);
% Quality check of multiple discriminator classifier
% It can be used to optimize the LDA-based classifier
%
% [Q,D]=qcmahal(XC);
%
% Q	sum of total distance Q=sum(D(:))
% D	Distance matrix (mahalanobis distance of mean to each other)
%	D(k,l) gives the distance of the k-th mean to the class 'l'
%
% [Q,D,R,COV]=qcmahal(XC);
% optimization to reduce the number of classes
% 
% Q is list of distances.  
% R is a list how the classes should be summarized
% COV gives the suggested class definition
%
%	Copyright (c) 1999-2002 by Alois Schloegl
%	a.schloegl@ieee.org	
%	28.02.2001 Version 1.13
%	10.10.2001 Version 1.14
%	30.12.2002 Version 1.15

NC = size(ECM);
if length(NC)<3,
        if iscell(ECM(1)),
                NC = [max(NC(1:2)),size(ECM{1})];
		tmp = ECM;
		ECM = zeros([NC(1),size(tmp{1})]);
                for k = 1:NC(1),
                        ECM(k,:,:) = tmp{k};
                end;

	elseif isfield(ECM,'COV') & isfield(ECM,'NN')
    		ECM = ECM.COV./ECM.NN; 
    		NC  = size(ECM);
        
        elseif isstruct(ECM),
                x = ECM;
                NC=[length(x.IR),size(x.IR{1})];
        elseif NC(1)==NC(2)
                ECM{1}=ECM;
        end;

elseif (length(NC)==3) & (NC(2)==NC(3)),
        
elseif isfield(ECM,'COV') & isfield(ECM,'NN')
        ECM = ECM.COV./ECM.NN; 
        NC  = size(ECM);
        
elseif 0; 
        %ECM = num2cell(ECM,[2,3]);
        for k = 1:NC(1),
                IR{k} = squeeze(ECM(k,:,:));
        end;
        ECM = IR;
else
        
end


n=NC; %size(COV);
D=zeros(NC(1),NC(1));

if n(1)==1, Q=0; D=0; R=[];return; end;
        
for k=1:NC(1),%size(COV,1);
        %[M,SD1,XC01,xc01,N1] = decovm(squeeze(COV(k,:,:)));
        [M,SD1,XC01,xc01,N1(k)] = decovm(squeeze(ECM(k,:,:)));
        M1(:,k)=M';
        ri{k} = XC01\eye(n(2)-1,n(2)-1);
        XC{k}=XC01;

	d = NC(2)-1;
	x.logSF(k) = log(N1(k)) - d/2*log(2*pi) - det(ri{k})/2;
	x.logSF2(k)= log(N1(k)) - d/2*log(2*pi) - log(det(ri{k}))/2;
	x.SF(k) = N1(k)/sqrt((2*pi)^d * det(ri{k}));
end;        
if min(N1)<10*d,
        fprintf(1,'Warning QCMAHAL: ratio #samples/#features: %f\n',min(N1)/d);
end;

for k=1:n(1);
        for l=1:n(1),
                %E = ([1,M1(:,l)'] * [-M1(:,k)'; eye(n(2)-1,n(2)-1)]); % 
                E = [M1(:,l)-M1(:,k)]';
                D(k,l) = sum((E*ri{k}).*E,2);

                if 0;l~=k,
                        m = (M1(:,k)-M1(:,l))';
                        w{k,l} = squeeze(COV(k,2:n(2),1)-COV(l,2:n(2),1)) / (XC{k} + XC{l});
                        w{k,l} = [-sum(COV([k,l],2:n(2),1),1)*w{k,l}; w{k,l}];
                end;
    		%LogLik(k,l) = x.logSF2(k) - D(k,l)/2;
      	end;
end;

% for 2 classes this is equivalent to the likelihood ratio (i.e. differences of the log-likehood) 
Q = sum(D(:))/(n(1)*(n(1)-1)); 
%sum(log(D(~eye(n(1)))));

if nargout>2, % recursion for optimization 
	[tmp,i] = min(D+diag(diag(D)+inf));
	[tmp,k] = min(tmp);
        i = i(k); 
        
        %tmpCOV = COV;
	%tmpCOV(i,:,:) = tmpCOV(i,:,:)+tmpCOV(k,:,:);
	%tmpCOV(k,:,:) = [];
        tmpCOV = ECM;
	tmpCOV{i} = tmpCOV{i}+tmpCOV{k};
	tmpCOV(k) = [];
	[q,d,R,tmp1]  = qcmahal(tmpCOV);

	%if any(q>Q), COV = tmp1; end; % if actual value is not optimal, keep reduced matrix
        
        R = [[i,k];R];
	Q = [Q,q];
end;

